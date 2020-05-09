{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.RNASeq.Functions.Clustering
    ( filterMatrix
    , selectFeature
    , selectFeature'
    , combineMatrices
    , spectral
    , mkKNNGraph
    , clustering
    , computeStability
    , pickResolution
    , segregateCells
    , mkExprTable
    ) where 

import Data.Singletons.Prelude (Elem)
import Data.Conduit.Zlib (multiple, ungzip, gzip)
import qualified Data.ByteString.Char8 as B
import Data.List.Ordered (nubSort)
import Data.Conduit.Internal (zipSinks)
import           Data.CaseInsensitive                 (original)
import           Bio.RealWorld.GENCODE (readGenes, Gene(..))
import Data.ByteString.Lex.Integral (packDecimal)
import Bio.Utils.Functions (scale)
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import qualified Data.HashSet as S
import Data.Binary (decodeFile)
import Control.Arrow (first, second)
import qualified Data.HashMap.Strict as M
import qualified Data.Text as T
import Data.Binary (encodeFile)
import Shelly (shelly, run_)
import Data.Char (toUpper)

import           Taiji.Pipeline.SC.RNASeq.Types (SCRNASeqConfig(..))
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Prelude
import Taiji.Utils
import Taiji.Utils.Plot
import Taiji.Utils.Plot.ECharts

combineMatrices :: SCRNASeqConfig config
                => [RNASeq S (Either
                        ( File '[RowName, Gzip] 'Tsv
                        , File '[ColumnName, Gzip] 'Tsv
                        , File '[Gzip] 'MatrixMarket ) 
                        ( File '[ColumnName, Gzip] 'Tsv
                        , File '[Gzip] 'Other ) 
                   )]
                -> ReaderT config IO
                    (Maybe (RNASeq S (File '[ColumnName, Gzip] 'Tsv, File '[Gzip] 'Other )))
combineMatrices [] = return Nothing
combineMatrices inputs = do
    dir <- asks _scrnaseq_output_dir >>= getPath . (<> (asDir "/Cluster/"))
    let output = dir <> "Merged.mat.gz"
        idxFl = dir <> "Merged.colnames.txt.gz"
    liftIO $ do
        (idx, mat) <- mergeMatrices <$> getMatrices inputs
        runResourceT $ runConduit $ yieldMany idx .| unlinesAsciiC .|
            gzip .| sinkFile idxFl
        saveMatrix output (fromJust . packDecimal) mat
        return $ Just $ (head inputs & eid .~ "Merged") &
            replicates._2.files .~
                ( location .~ idxFl $ emptyFile 
                , location .~ output $ emptyFile )
  where
    getMatrices fls = forM fls $ \x -> do
        let e = B.pack $ T.unpack $ x^.eid
        case x^.replicates._2.files of
            Left (idxFl, rownameFl, matFl) -> do
                rownames <- readIndex rownameFl
                idx <- V.map (B.map toUpper) <$> readIndex idxFl
                mat <- transformation (addPrefix e) <$>
                    mkSpMatrixMM (matFl^.location) rownames
                return (idx, mat)
            Right (idxFl, matFl) -> do
                idx <- V.map (B.map toUpper) <$> readIndex idxFl
                mat <- transformation (addPrefix e) <$>
                    mkSpMatrix readInt (matFl^.location)
                return (idx, mat)
    addPrefix x = mapC $ first (\n -> x <> "+" <> n) 
    readIndex x = runResourceT $ runConduit $ sourceFile (x^.location) .|
        multiple ungzip .| linesUnboundedAsciiC .|
        mapC (head . B.words) .| sinkVector

selectFeature :: SCRNASeqConfig config
              => FilePath
              -> RNASeq S (File '[ColumnName, Gzip] 'Tsv, File '[Gzip] 'Other)
              -> ReaderT config IO [Int]
selectFeature prefix input = do
    dir <- asks ((<> asDir prefix) . _scrnaseq_output_dir) >>= getPath
    let output = dir <> "/feature_selection.html"
        fl = input^.replicates.traversed.files._2.location
    liftIO $ do
        sp <- mkSpMatrix readInt fl
        mv <- meanVariance $ transformation (mapC $ second $ map (second fromIntegral)) sp
        let bins = getBins 20 $ U.toList $ U.map fst mv
            logVMR = U.map (\(m,v) -> if m == 0 then 0 else log' $ v / m) mv
            (hi, lo) = partition ((>0.5) . snd) $ scaleDispersion bins logVMR
            hi' = map (\(i, _) -> (log $ fst $ mv U.! i, logVMR U.! i)) hi
            lo' = map (\(i, _) -> (log $ fst $ mv U.! i, logVMR U.! i)) lo
            p = scatter' [("High", hi'), ("Low", lo')]
        savePlots output [] [p]
        return $ map fst lo
  where
    log' x = if x == 0 then -5 else log x

selectFeature' :: SCRNASeqConfig config
               => FilePath
               -> RNASeq S (File '[ColumnName, Gzip] 'Tsv, File '[Gzip] 'Other)
               -> ReaderT config IO [Int]
selectFeature' prefix input = do
    dir <- asks ((<> asDir prefix) . _scrnaseq_output_dir) >>= getPath
    let fl = input^.replicates.traversed.files._2.location
    liftIO $ do
        sp <- mkSpMatrix readInt fl
        counts <- do
            v <- UM.replicate (_num_col sp) 0
            runResourceT $ runConduit $ streamRows sp .| concatMapC snd .|
                mapC fst .| mapM_C (UM.unsafeModify v (+1))
            U.unsafeFreeze v
        let (zeros, nonzeros) = U.partition ((==0) . snd) $
                U.zip (U.enumFromN 0 (U.length counts)) counts
            (i, v) = U.unzip nonzeros
            idx = U.toList $ fst $ U.unzip $ U.filter ((>2) . snd) $ U.zip i $ scale v
        return $ idx ++ U.toList (fst $ U.unzip zeros)

scaleDispersion :: [[Int]] -> U.Vector Double -> [(Int, Double)]
scaleDispersion grp xs = concatMap f grp
  where
    f idx = zip idx $ U.toList $ scale $ U.fromList $ map (xs U.!) idx

getBins :: Int -> [Double] -> [[Int]]
getBins nBin xs = filter (not . null) $ go [] (lo + step) $
    sortBy (comparing fst) $ zip xs [0..]
  where
    go acc x' ((x, i):xs) | x <= x' = go (i:acc) x' xs
                          | otherwise = acc : go [] (x'+step) ((x,i):xs)
    go acc _ _ = [acc]
    (lo, hi) = (minimum xs, maximum xs)
    step = (hi - lo) / fromIntegral nBin


filterMatrix :: SCRNASeqConfig config
             => FilePath
             -> RNASeq S (File '[ColumnName, Gzip] 'Tsv, File '[Gzip] 'Other)
             -> ReaderT config IO (RNASeq S (File '[] 'Tsv, File '[Gzip] 'Other))
filterMatrix prefix input = do
    dir <- asks ((<> asDir prefix) . _scrnaseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_norm_filt.mat.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    let rownames = printf "%s/%s_rep%d_rownames.txt" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ ( \(_, fl) -> liftIO $ do
        sp <- mkSpMatrix readInt $ fl^.location
        let (r, c) = (_num_row sp, _num_col sp)
        _ <- runResourceT $ runConduit $ streamRows sp .| mapC f .|
            zipSinks (mapC fst .| unlinesAsciiC .| sinkFile rownames)
                (mapC snd .| sinkRows r c toShortest output)
        --filterCols output removed $ fl^.location
        return ( location .~ rownames $ emptyFile
               , location .~ output $ emptyFile )
        )
  where
    f (nm, xs) = ( nm <> "\t" <> fromJust (packDecimal totalReads)
                 , (nm, map normalize xs) )
      where
        totalReads = foldl1' (+) $ map snd xs
        normalize (i, x) = (i, log1p $ fromIntegral x / (fromIntegral totalReads / 10000))
        log1p x = log $ 1 + x :: Double

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
        shelly $ run_ "taiji-utils" $ [ "reduce"
            , "--distance", "cosine"
            , T.pack $ fl^.location, T.pack output ] ++ maybe []
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
           -> (Double, RNASeq S (File '[] 'Tsv, File '[] 'Other, File '[] Tsv))
           -> ReaderT config IO (RNASeq S (File '[] 'Other))
clustering prefix (res, input) = do
    tmp <- asks _scrnaseq_tmp_dir
    dir <- asks ((<> asDir ("/" ++ prefix)) . _scrnaseq_output_dir) >>= getPath
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
    f (i, (bc, dep), [d1,d2]) = Cell i (d1,d2) bc $ readInt dep
    f _ = error "formatting error"
{-# INLINE clustering' #-}

pickResolution :: [(Double, Double, Double)]
               -> Double
pickResolution xs = case filter (\x -> x^._3 > 0.8) (take 5 xs') of
    [] -> head xs' ^. _1
    x -> maximumBy (comparing (^._2)) x ^. _1
  where
    xs' = sortBy (flip (comparing (^._3))) xs

computeStability :: SCRNASeqConfig config
                 => (Double, RNASeq S (File '[] 'Tsv, File '[] 'Other, File '[] Tsv))
                 -> ReaderT config IO (Double, Double, Double)
computeStability (res, input) = do
    tmp <- asks _scrnaseq_tmp_dir
    optimizer <- asks _scrnaseq_cluster_optimizer 
    let knn = input^.replicates.traversed.files._2.location
    liftIO $ withTemp tmp $ \tmpFl -> do
        shelly $ run_ "taiji-utils"
            [ "clust", T.pack knn, T.pack tmpFl
            , "--stability"
            , "--res", T.pack $ show res
            , "--optimizer"
            , case optimizer of
                RBConfiguration -> "RB"
                CPM -> "CPM"
            ]
        [num, st] <- words . head . lines <$> readFile tmpFl
        return (res, read num, read st)

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
        _ -> undefined

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
segregateCells :: SCRNASeqConfig config
               => FilePath   -- ^ Dir
               -> RNASeq S (File '[ColumnName, Gzip] 'Tsv, File '[Gzip] 'Other)
               -> File tag' 'Other   -- Clusters
               -> ReaderT config IO 
                  ( File '[ColumnName, Gzip] 'Tsv
                  , [RNASeq S (File '[Gzip] 'Other)] )
segregateCells prefix matFl clFl = do
    dir <- asks _scrnaseq_output_dir >>= getPath . (<> (asDir prefix))
    liftIO $ do
        cls <- decodeFile $ clFl^.location
        mat <- mkSpMatrix id $ matFl ^. replicates._2.files._2.location
        let mkSink CellCluster{..} = filterC ((`S.member` ids) . fst) .|
                (sinkRows (S.size ids) (_num_col mat) id output >> return res)
              where
                ids = S.fromList $ map _cell_barcode _cluster_member
                output = dir <> B.unpack _cluster_name <> ".mat.gz"
                res = matFl & eid .~ T.pack (B.unpack _cluster_name)
                            & replicates._2.files .~ (location .~ output $ emptyFile)
        res <- runResourceT $ runConduit $
            streamRows mat .| sequenceSinks (map mkSink cls)
        return (matFl^.replicates._2.files._1, res)

-- | Combine expression data into a table and output
mkExprTable :: SCRNASeqConfig config
            => FilePath
            -> Maybe ( File '[ColumnName, Gzip] 'Tsv
               , [RNASeq S (File '[Gzip] 'Other)] )
            -> ReaderT config IO (Maybe (FilePath, FilePath, FilePath))
mkExprTable _ Nothing = return Nothing
mkExprTable prefix (Just (geneFl, inputs)) = do
    dir <- asks _scrnaseq_output_dir >>= getPath . (<> asDir prefix)
    liftIO $ do
        let output1 = dir ++ "/gene_expression.tsv"
            output2 = dir ++ "/gene_specificity.tsv"
            output3 = dir ++ "/gene_specificity_pvalue.tsv"
        genes <- runResourceT $ runConduit $
            sourceFile (geneFl^.location) .| multiple ungzip .|
            linesUnboundedAsciiC .| sinkList
        mat <- fmap transpose $ forM inputs $ \input -> fmap U.toList $
            computeRAS $ input^.replicates._2.files.location 
        let expr = DF.mkDataFrame (map (T.pack . B.unpack) genes)
                (map (^.eid) inputs) mat
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
