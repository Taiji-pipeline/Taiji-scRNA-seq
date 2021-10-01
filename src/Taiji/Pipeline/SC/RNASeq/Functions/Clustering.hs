{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.RNASeq.Functions.Clustering
    ( combineMatrices
    , batchCorrection
    , mkKNNGraph
    , plotClusters
    , segregateCells
    , mkExprTable
    , old_spectral
    ) where 

import Data.Singletons.Prelude (Elem)
import Data.Conduit.Zlib (multiple, ungzip, gzip)
import qualified Data.ByteString.Char8 as B
import Data.List.Ordered (nubSort)
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import qualified Data.HashSet as S
import Data.Binary (decodeFile)
import Control.Arrow (first)
import qualified Data.HashMap.Strict as M
import qualified Data.Text as T
import Data.Binary (encodeFile)
import Shelly (shelly, run_)
import Data.Char (toUpper)

import           Taiji.Pipeline.SC.RNASeq.Types (SCRNASeqConfig(..), SCRNASeq)
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Prelude
import Taiji.Utils
import Taiji.Utils.Plot
import Taiji.Utils.Plot.ECharts

combineMatrices :: SCRNASeqConfig config
                => [SCRNASeq S (Either
                        ( File '[RowName, Gzip] 'Tsv
                        , File '[ColumnName, Gzip] 'Tsv
                        , File '[Gzip] 'MatrixMarket ) 
                        ( File '[ColumnName, Gzip] 'Tsv
                        , File '[Gzip] 'Matrix) 
                   )]
                -> ReaderT config IO
                    (Maybe (SCRNASeq S (File '[ColumnName, Gzip] 'Tsv, File '[Gzip] 'Matrix )))
combineMatrices [] = return Nothing
combineMatrices inputs = do
    dir <- asks _scrnaseq_output_dir >>= getPath . (<> (asDir "/Cluster/"))
    let output = dir <> "Merged.mat.gz"
        idxFl = dir <> "Merged.colnames.txt.gz"
    liftIO $ do
        (idx, mat) <- mergeMatrices <$> getMatrices inputs
        zeros <- U.toList . U.imapMaybe (\i x -> if x < minCells then Just i else Nothing) <$>
            colSum (fmap (const (1 :: Int)) mat)
        let mat' = deleteCols zeros mat
            idx' = V.ifilter (\i _ -> not $ i `elem` zeros) idx
        runResourceT $ runConduit $ yieldMany idx' .| unlinesAsciiC .|
            gzip .| sinkFile idxFl
        saveMatrix output (fromJust . packDecimal) mat'
        return $ Just $ (head inputs & eid .~ "Merged") &
            replicates._2.files .~
                ( location .~ idxFl $ emptyFile 
                , location .~ output $ emptyFile )
  where
    minCells = 5
    getMatrices fls = forM fls $ \x -> do
        let e = B.pack $ printf "%s_%d" (T.unpack $ x^.eid) (x^.replicates._1)
        case x^.replicates._2.files of
            Left (idxFl, rownameFl, matFl) -> do
                rownames <- readIndex rownameFl
                idx <- V.map (B.map toUpper) <$> readIndex idxFl
                mat <- mapRows (addPrefix e) <$>
                    mkSpMatrixMM (matFl^.location) rownames
                return (idx, mat)
            Right (idxFl, matFl) -> do
                idx <- V.map (B.map toUpper) <$> readIndex idxFl
                mat <- mapRows (addPrefix e) <$>
                    mkSpMatrix readInt (matFl^.location)
                return (idx, mat)
    addPrefix x = first (\n -> x <> "+" <> n) 
    readIndex x = runResourceT $ runConduit $ sourceFile (x^.location) .|
        multiple ungzip .| linesUnboundedAsciiC .|
        mapC (head . B.words) .| sinkVector

-- | Reduce dimensionality using spectral clustering
old_spectral :: (Elem 'Gzip tags ~ 'True, SCRNASeqConfig config)
         => FilePath  -- ^ directory
         -> Maybe Int  -- ^ seed
         -> SCRNASeq S (File '[ColumnName, Gzip] 'Tsv, File tags 'Matrix )
         -> ReaderT config IO (SCRNASeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
old_spectral prefix seed input = do
    dir <- asks ((<> asDir prefix) . _scrnaseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_spectral.tsv.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
        rownames = printf "%s/%s_rep%d_rownames.txt" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \(_, fl) -> do
        sp <- mkSpMatrix readInt $ fl^.location
        _ <- runResourceT $ runConduit $ streamRows sp .| mapC f .| unlinesAsciiC .| sinkFile rownames
        shelly $ run_ "taiji-utils" $ [ "reduce"
            , "--distance", "pca"
            , T.pack $ fl^.location, T.pack output ] ++ maybe []
            (\x -> ["--seed", T.pack $ show x]) seed
        return ( location .~ rownames $ emptyFile, location .~ output $ emptyFile)
        )
  where
    f (nm, xs) = nm <> "\t" <> fromJust (packDecimal totalReads)
      where
        totalReads = foldl1' (+) $ map snd xs

batchCorrection :: SCRNASeqConfig config
                => String -> SCRNASeq S (File '[] 'Tsv, File '[Gzip] 'Tsv)
                -> ReaderT config IO (SCRNASeq S (File '[] 'Tsv, File '[Gzip] 'Tsv))
batchCorrection prefix input = do
    dir <- asks ((<> asDir prefix) . _scrnaseq_output_dir) >>= getPath
    asks _scrnaseq_batch_info >>= \case
        Nothing -> return input
        Just batchFl -> do
            let output = printf "%s/%s_rep%d_spectral_corrected.tsv.gz" dir
                    (T.unpack $ input^.eid) (input^.replicates._1)
            input & replicates.traversed.files %%~ liftIO . ( \(rownames, fl) -> do
                idToBatchMap <- M.fromListWith undefined <$> readBatchInfo batchFl
                let f x = let (i, r) = B.breakEnd (=='_') x
                          in case M.lookup (B.init i) idToBatchMap of
                              Nothing -> Nothing
                              Just (l, g) -> Just (l <> r, g)
                labels <- map (f . B.init . fst . B.breakEnd (=='+') . head . B.split '\t') . B.lines <$>
                    B.readFile (rownames^.location)
                if (all isNothing labels)
                    then return (rownames, fl)
                    else do
                        readData (fl^.location) >>= batchCorrect labels >>= writeData output
                        return (rownames, location .~ output $ emptyFile)
                )
  where
    readData fl = runResourceT $ runConduit $
        sourceFile fl .| multiple ungzip .| linesUnboundedAsciiC .|
        mapC (U.fromList . map readDouble . B.split '\t') .| sinkVector
    writeData output vec = runResourceT $ runConduit $ yieldMany vec .|
        mapC (B.intercalate "\t" . map toShortest . U.toList) .|
        unlinesAsciiC .| gzip .| sinkFile output

mkKNNGraph :: SCRNASeqConfig config
           => FilePath
           -> SCRNASeq S [(a, File '[Gzip] 'Tsv)]
           -> ReaderT config IO (SCRNASeq S (a, File '[] 'Other, File '[] Tsv))
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

plotClusters :: SCRNASeqConfig config
             => ( [(Double, (Int, Double, Double))]
                , [(Double, ([FilePath], FilePath))]
                , Maybe (SCRNASeq S (File '[] 'Tsv, File '[] 'Other, File '[] Tsv)) )
             -> ReaderT config IO (Maybe (SCRNASeq S (File '[] 'Other)))
plotClusters (_, _, Nothing) = return Nothing
plotClusters (params, clFl, Just knn) = do
    dir <- asks _scrnaseq_output_dir >>= getPath . (<> "/Cluster/")
    resAuto <- liftIO $ optimalParam (dir <> "/Clustering_parameters.html") params
    res <- fromMaybe resAuto <$> asks _scrnaseq_cluster_resolution 
    let clOutput = printf "%s/%s_rep%d_clusters.bin" dir (T.unpack $ knn^.eid)
            (knn^.replicates._1)
    fmap Just $ knn & replicates.traversed.files %%~ liftIO . ( \(idx, _, umap) -> do
        let sourceCells = getZipSource $ (,,) <$>
                ZipSource (iterateC succ 0) <*>
                ZipSource (sourceFile (idx^.location) .| linesUnboundedAsciiC .| mapC g) <*>
                ZipSource ( sourceFile (umap^.location) .|
                    linesUnboundedAsciiC .|
                    mapC (map readDouble . B.split '\t') )
            (perturbed, cl) = fromJust $ lookup res clFl
        clusters <- decodeFile cl
        perturbedCls <- mapM decodeFile perturbed
        let (clNames, repro) = unzip $ zipWith (\i x -> (T.pack $ "C" <> show i, x)) [1::Int ..] $
                computeReproducibility clusters perturbedCls
            figRepro = line $ DF.mkDataFrame ["reproducibility"] clNames [repro]
        cells <- runResourceT $ runConduit $ sourceCells .| mapC f .| sinkVector
        let cellCluster = zipWith3 (\i x y -> CellCluster (B.pack $ "C" ++ show i) x $ Just y)
                [1::Int ..] (map (map (cells V.!)) clusters) repro
        encodeFile clOutput cellCluster
        let clViz = printf "%s/%s_rep%d_cluster.html" dir (T.unpack $ knn^.eid)
                (knn^.replicates._1)
            (nms, num_cells) = unzip $ map (\(CellCluster nm cs _) ->
                (T.pack $ B.unpack nm, fromIntegral $ length cs)) cellCluster
            plt = stackBar $ DF.mkDataFrame ["number of cells"] nms [num_cells]
            compos = composition cellCluster
        clusters' <- sampleCells cellCluster
        savePlots clViz [] $ visualizeCluster clusters' ++
            figRepro : clusterComposition compos : tissueComposition compos : [plt]
        return $ location .~ clOutput $ emptyFile )
  where
    f (i, (bc, dep), [d1,d2]) = Cell i (d1,d2) bc $ readInt dep
    f _ = error "formatting error"
    g x = case B.split '\t' x of
        [a, b] -> (a, b)
        [a] -> (a, "0")
        _ -> error "formatting error"


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
               -> SCRNASeq S (File '[ColumnName, Gzip] 'Tsv, File '[Gzip] 'Matrix)
               -> File tag' 'Other   -- Clusters
               -> ReaderT config IO 
                  ( File '[ColumnName, Gzip] 'Tsv
                  , [SCRNASeq S (File '[Gzip] 'Matrix)] )
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
               , [SCRNASeq S (File '[Gzip] 'Matrix)] )
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
