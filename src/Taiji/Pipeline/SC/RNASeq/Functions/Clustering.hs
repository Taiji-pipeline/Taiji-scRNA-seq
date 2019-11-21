{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.RNASeq.Functions.Clustering
    ( filterMatrix
    , spectral
    , mkKNNGraph
    , clustering
    , subMatrix
    , mkExprTable
    ) where 

import Data.Singletons.Prelude (Elem)
import qualified Data.ByteString.Char8 as B
import Data.List.Ordered (nubSort)
import           Data.CaseInsensitive                 (original)
import           Bio.RealWorld.GENCODE (readGenes, Gene(..))
import Data.ByteString.Lex.Integral (packDecimal)
import Bio.Utils.Functions (scale, filterFDR)
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import qualified Data.HashSet as S
import Data.Binary (decodeFile)
import Control.Arrow (first)
import qualified Data.HashMap.Strict as M
import qualified Data.Text as T
import Data.Binary (encodeFile)
import Shelly (shelly, run_)
import qualified Data.Matrix as Mat

import           Taiji.Pipeline.SC.RNASeq.Types
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Prelude
import Taiji.Utils
import Taiji.Utils.Plot
import Taiji.Utils.Plot.ECharts
import Taiji.Pipeline.SC.ATACSeq.Functions.Utils (computeSS, computeRAS, computeCDF)

filterMatrix :: (Elem 'Gzip tags ~ 'True, SCRNASeqConfig config)
             => FilePath
             -> RNASeq S (File tags 'Other)
             -> ReaderT config IO (RNASeq S (File '[] 'Tsv, File tags 'Other))
filterMatrix prefix input = do
    dir <- asks ((<> asDir prefix) . _scrnaseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_filt.mat.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
        rownames = printf "%s/%s_rep%d_rownames.txt" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ ( \fl -> liftIO $ do
        sp <- mkSpMatrix readInt $ fl^.location
        runResourceT $ runConduit $
            streamRows sp .| mapC f .| unlinesAsciiC .| sinkFile rownames
        counts <- do
            v <- UM.replicate (_num_col sp) 0
            runResourceT $ runConduit $ streamRows sp .| concatMapC snd .|
                mapC fst .| mapM_C (UM.unsafeModify v (+1))
            U.unsafeFreeze v
        let (zeros, nonzeros) = U.partition ((==0) . snd) $
                U.zip (U.enumFromN 0 (U.length counts)) counts
            (i, v) = U.unzip nonzeros
            idx = U.toList $ fst $ U.unzip $ U.filter ((>1.65) . snd) $ U.zip i $ scale v
        filterCols output (idx ++ U.toList (fst $ U.unzip zeros)) $ fl^.location
        return ( location .~ rownames $ emptyFile
               , location .~ output $ emptyFile ) )
  where
    f (nm, xs) = nm <> "\t" <> fromJust (packDecimal $ foldl1' (+) $ map snd xs)

-- | Reduce dimensionality using spectral clustering
spectral :: (Elem 'Gzip tags ~ 'True, SCRNASeqConfig config)
         => FilePath  -- ^ directory
         -> Maybe Int  -- ^ seed
         -> RNASeq S (a, File tags 'Other)
         -> ReaderT config IO (RNASeq S (a, File '[Gzip] 'Tsv))
spectral prefix seed input = do
    dir <- asks ((<> asDir prefix) . _scrnaseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_spectral.tsv.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \(rownames, fl) -> do
        shelly $ run_ "taiji-utils" $ ["reduce", "--method", "spectral"
            , "--distance", "cosine"
            , T.pack $ fl^.location, T.pack output] ++ maybe []
            (\x -> ["--seed", T.pack $ show x]) seed
        return (rownames, location .~ output $ emptyFile)
        )

mkKNNGraph :: SCRNASeqConfig config
           => FilePath
           -> RNASeq S [(File '[] 'Tsv, File '[Gzip] 'Tsv)]
           -> ReaderT config IO (RNASeq S (File '[] 'Tsv, File '[] 'Other, File '[] Tsv))
mkKNNGraph prefix input = do
    dir <- asks ((<> asDir ("/" ++ prefix)) . _scrnaseq_output_dir) >>= getPath
    let output_knn = printf "%s/%s_rep%d_knn.npz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        output_umap = printf "%s/%s_rep%d_umap.txt" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fls -> do
        shelly $ run_ "taiji-utils" $
            [ "knn"
            , T.pack $ intercalate "," $ map (^.location) $ snd $ unzip fls
            , T.pack output_knn
            , "-k", "50"
            , "--embed", T.pack output_umap
            , "--thread", "4" ]
        return ( head $ fst $ unzip fls
               , location .~ output_knn $ emptyFile
               , location .~ output_umap $ emptyFile )
        )

clustering :: SCRNASeqConfig config
           => FilePath
           -> RNASeq S (File '[] 'Tsv, File '[] 'Other, File '[] Tsv)
           -> ReaderT config IO (RNASeq S (File '[] 'Other))
clustering prefix input = do
    tmp <- asks _scrnaseq_tmp_dir
    dir <- asks ((<> asDir ("/" ++ prefix)) . _scrnaseq_output_dir) >>= getPath
    res <- asks _scrnaseq_cluster_resolution 
    optimizer <- asks _scrnaseq_cluster_optimizer 
    let output = printf "%s/%s_rep%d_clusters.bin" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        plotFl = printf "%s/%s_rep%d_clusters.html" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> do
        cls <- clustering' tmp optimizer res fl
        encodeFile output cls
        plotClusters plotFl cls
        return $ location .~ output $ emptyFile )

clustering' :: Maybe FilePath      -- ^ temp dir
            -> Optimizer 
            -> Double  -- ^ Clustering resolution
            -> (File '[] 'Tsv, File '[] 'Other, File '[] Tsv)
            -> IO [CellCluster]
clustering' dir method res (idx, knn, umap) = withTemp dir $ \tmpFl -> do
      let sourceCells = getZipSource $ (,,) <$>
              ZipSource (iterateC succ 0) <*>
              ZipSource seqDepthC <*>
              ZipSource ( sourceFile (umap^.location) .| linesUnboundedAsciiC .|
                mapC (map readDouble . B.split '\t') )
      shelly $ run_ "taiji-utils"
          [ "clust", T.pack $ knn^.location, T.pack tmpFl
          , "--res"
          , T.pack $ show res
          , "--optimizer"
          , case method of
              RBConfiguration -> "RB"
              CPM -> "CPM"
          ]
      cells <- runResourceT $ runConduit $ sourceCells .| mapC f .| sinkVector
      clusters <- map (map readInt . B.split ',') . B.lines <$> B.readFile tmpFl
      return $ zipWith (\i -> CellCluster $ B.pack $ "C" ++ show i) [1::Int ..] $
          map (map (cells V.!)) clusters
  where
    seqDepthC = sourceFile (idx^.location) .| linesUnboundedAsciiC .|
        mapC ((\[a,b] -> (a,b)) . B.split '\t')
    f (i, (bc, dep), [d1,d2,d3,d4,d5]) = Cell i (d1,d2) (d3,d4,d5) bc $ readInt dep
    f _ = error "formatting error"
{-# INLINE clustering' #-}

plotClusters :: FilePath
             -> [CellCluster]
             -> IO ()
plotClusters output input = do
    let (nms, num_cells) = unzip $ map (\(CellCluster nm cells) ->
            (T.pack $ B.unpack nm, fromIntegral $ length cells)) input
        plt = stackBar $ DF.mkDataFrame ["number of cells"] nms [num_cells]
    clusters <- sampleCells input
    case visualizeCluster clusters of
        [p] -> savePlots output [] $ [p, plt]
        [p1,p2] -> do
            let compos = composition input
            savePlots output [] $ [p1, p2, plt,
                clusterComposition compos, tissueComposition compos]

-- | Compute the normalized tissue composition for each cluster.
tissueComposition :: DF.DataFrame Int -> EChart
tissueComposition = stackBar . DF.map round' . DF.mapCols normalize .
    DF.transpose . DF.mapCols normalize . DF.map fromIntegral
  where
    round' x = fromIntegral (round $ x * 1000 :: Int) / 1000
    normalize :: V.Vector Double -> V.Vector Double
    normalize xs | V.all (==0) xs = xs
                 | otherwise = V.map (/V.sum xs) xs

-- | Compute the cluster composition for each tissue.
clusterComposition :: DF.DataFrame Int -> EChart
clusterComposition = stackBar . DF.map round' . DF.mapCols normalize .
    DF.map fromIntegral
  where
    round' x = fromIntegral (round $ x * 1000 :: Int) / 1000
    normalize :: V.Vector Double -> V.Vector Double
    normalize xs | V.all (==0) xs = xs
                 | otherwise = V.map (/V.sum xs) xs

-- | Compute the cluster x tissue composition table.
composition :: [CellCluster] -> DF.DataFrame Int
composition clusters = DF.mkDataFrame rownames colnames $
    flip map rows $ \x -> map (\i -> M.lookupDefault 0 i x) colnames
  where
    (rownames, rows) = unzip $ flip map clusters $ \CellCluster{..} ->
        ( T.pack $ B.unpack _cluster_name
        , M.fromListWith (+) $ map (\x -> (getName x, 1)) _cluster_member )
    colnames = nubSort $ concatMap M.keys rows
    getName Cell{..} =
        let prefix = fst $ B.breakEnd (=='+') _cell_barcode
        in if B.null prefix then "" else T.pack $ B.unpack $ B.init prefix

-- | Extract cluster submatrix
subMatrix :: SCRNASeqConfig config
          => FilePath   -- ^ Dir
          -> [RNASeq S (File tags 'Other)]   -- ^ Matrices
          -> File tag' 'Other   -- Clusters
          -> ReaderT config IO [RNASeq S (File tags 'Other)]
subMatrix prefix mats clFl = do
    dir <- asks _scrnaseq_output_dir >>= getPath . (<> (asDir prefix))
    liftIO $ do
        cls <- decodeFile $ clFl^.location
        mat <- mkSpMatrix id $ head mats ^. replicates._2.files.location
        let mkSink CellCluster{..} = filterC ((`S.member` ids) . fst) .|
                (sinkRows (S.size ids) (_num_col mat) id output >> return res)
              where
                ids = S.fromList $ map _cell_barcode _cluster_member
                output = dir <> B.unpack _cluster_name <> ".mat.gz"
                res = head mats & eid .~ T.pack (B.unpack _cluster_name)
                                & replicates._2.files.location .~ output
        runResourceT $ runConduit $ streamMatrices id mats .|
            sequenceSinks (map mkSink cls)

-- | Stream rows and add sample id to barcodes.
streamMatrices :: (B.ByteString -> a)   -- ^ Element decoder
               -> [RNASeq S (File tags 'Other)]
               -> ConduitT () (Row a) (ResourceT IO) ()
streamMatrices decoder inputs = forM_ inputs $ \input -> do
    mat <- liftIO $ mkSpMatrix decoder $ input^.replicates._2.files.location
    let f x = B.pack (T.unpack $ input^.eid) <> "+" <> x
    streamRows mat .| mapC (first f)

-- | Combine expression data into a table and output
mkExprTable :: SCRNASeqConfig config
            => FilePath
            -> [RNASeq S (File '[Gzip] 'Other)]
            -> ReaderT config IO (Maybe (FilePath, FilePath, FilePath))
mkExprTable _ [] = return Nothing
mkExprTable prefix inputs = do
    dir <- asks _scrnaseq_output_dir >>= getPath . (<> asDir prefix)
    geneNames <- asks _scrnaseq_annotation >>=
        liftIO . fmap (map (original . geneName)) . readGenes
    liftIO $ do
        let output1 = dir ++ "/gene_expression.tsv"
            output2 = dir ++ "/gene_specificity.tsv"
            output3 = dir ++ "/gene_specificity_pvalue.tsv"
        mat <- fmap transpose $ forM inputs $ \input -> fmap U.toList $
            computeRAS $ input^.replicates._2.files.location 
        let (genes, vals) = unzip $ map combine $ groupBy ((==) `on` fst) $
                sortBy (comparing fst) $ zip geneNames mat
            expr = DF.mkDataFrame (map (T.pack . B.unpack) genes)
                (map (^.eid) inputs) vals
            ss = computeSS expr

        DF.writeTable output1 (T.pack . show) expr
        DF.writeTable output2 (T.pack . show) ss
        cdf <- computeCDF expr
        DF.writeTable output3 (T.pack . show) $ DF.map (lookupP cdf) ss
        return $ Just (output1, output2, output3)
  where
    combine xs = (head gene, foldl1' (zipWith max) vals)
      where
        (gene, vals) = unzip xs
    lookupP (vec, res, n) x | p == 0 = 1 / n
                            | otherwise = p
      where
        p = vec V.! i
        i = min (V.length vec - 1) $ truncate $ x / res 

{-
-- | Get cell-specific genes
specificGene :: SCRNASeqConfig config
             => FilePath
             -> (FilePath, FilePath, FilePath)
             -> ReaderT config IO [(T.Text, File '[] 'Tsv)]
specificGene prefix (_, scoreFl, pvalueFl) = do
    dir <- asks ((<> asDir prefix) . _scrnaseq_output_dir) >>= getPath
    liftIO $ do
        scores <- DF.readTable scoreFl
        pvalues <- DF.readTable pvalueFl
        let table = DF.zip scores pvalues
            genes = V.fromList $ DF.rowNames table
            cls = zip (DF.colNames table) $ Mat.toColumns $
                DF._dataframe_data table
        forM cls $ \(nm, vec) -> do
            let output = dir <> T.unpack nm <> "_markers.tsv"
                (gs, sc, pval) = unzip3 $ fst $ unzip $ V.toList $
                    filterFDR fdr $ V.zipWith (\g (a,b) -> ((g,a,b), b)) genes vec 
            DF.writeTable output (T.pack . show) $
                DF.mkDataFrame gs ["specificity", "pvalue"] $
                transpose [sc, pval]
            return (nm, location .~ output $ emptyFile)
  where
    fdr = 0.001
    -}