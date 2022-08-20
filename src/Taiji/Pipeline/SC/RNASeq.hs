{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.SC.RNASeq ( builder) where

import           Control.Workflow
import qualified Data.Text as T
import Data.List.Ordered (nubSort)

import Taiji.Prelude
import           Taiji.Utils
import           Taiji.Pipeline.SC.RNASeq.Functions

builder :: Builder ()
builder = do
    node "Read_Input" [| readInput |] $ return ()

    uNode "Get_Fastq" [| return . getFastq |]
    nodePar "Barcode_Stat" [| getValidBarcode |] $ return ()
    nodePar "Demultiplex" [| demultiplex |] $ nCore .= 3
    path ["Read_Input", "Get_Fastq", "Barcode_Stat", "Demultiplex"]

    uNode "Get_Demulti_Fastq" [| \(input, fq) -> return $
        getDemultiplexedFastq input ++ fq |]
    ["Read_Input", "Demultiplex"] ~> "Get_Demulti_Fastq"

    node "Make_Index" [| mkIndex |] $ return ()
    nodePar "Align" [| tagAlign |] $ do
        nCore .= 8
        memory .= 50
    nodePar "Filter_Bam" [| filterNameSortBam |] $ nCore .= 2
    nodePar "Quantification" [| quantification |] $ do
        memory .= 8 
        nCore .= 3
    nodePar "Filter_Cell" [| filterCells |] $ return ()
    nodePar "Remove_Doublet" [| removeDoublet |] $ return ()
    path ["Get_Demulti_Fastq", "Make_Index", "Align", "Filter_Bam"
        , "Quantification", "Filter_Cell", "Remove_Doublet"]

    node "QC" [| plotQC |] $ return ()
    ["Quantification"] ~> "QC"

    node "Merge_Matrix" [| \(inputs, mats) -> do
        let a = getMatrix inputs & mapped.replicates._2.files %~ Left
            b = mats & mapped.replicates._2.files %~ Right . fst
        combineMatrices $ a ++ b
        |] $ return ()
    ["Read_Input", "Remove_Doublet"] ~> "Merge_Matrix"


--------------------------------------------------------------------------------
-- Normalization
--------------------------------------------------------------------------------
    node "Merged_Normalization" [| \case
        Nothing -> return Nothing
        Just input -> Just <$> logNormalize input
        |] $ return ()
    path ["Merge_Matrix", "Merged_Normalization"]
{-
    node "Merged_Fit_NB" [| \case
        Nothing -> return Nothing
        Just input -> Just <$> fitNB input
        |] $ return ()
    uNode "Merged_Normalization_Prep" [| \case
        (Just mat, Just param) -> fmap (map (\x ->
            let (i, j) = x^.replicates._2.files._2
            in x & eid .~ (T.pack $ printf "%d_%d" i j) ) . split) $
            param & replicates.traversed.files %%~ liftIO . ( \fl -> do
                let matFl = mat^.replicates._2.files._2
                nCells <- fmap _num_row $ mkSpMatrix id $ matFl^.location
                let go i | j <= nCells = (i, j) : go j
                         | i < nCells = [(i, nCells)]
                         | otherwise = []
                      where j = i + batchSize
                    batchSize = 50000
                return $ zip3 (repeat fl) (go 0) (repeat matFl) )
        _ -> return []
        |]
    nodePar "Merged_Normalization" 'normalization $ return ()
    node "Merged_Feature_Selection" [| selectFeatures 3000 |] $ return ()
    path ["Merge_Matrix", "Merged_Fit_NB"]
    ["Merge_Matrix", "Merged_Fit_NB"] ~> "Merged_Normalization_Prep"
    path ["Merged_Normalization_Prep", "Merged_Normalization"]
    ["Merge_Matrix", "Merged_Normalization", "Merged_Fit_NB"] ~> "Merged_Feature_Selection"
-}

--------------------------------------------------------------------------------
-- TEST
--------------------------------------------------------------------------------
    node "Merged_Reduce_Dims" [| \case
        Just input -> do
            let prefix = "/TEST/"
            fmap Just $ old_spectral prefix Nothing input
        _ -> return Nothing
        |] $ return ()
    ["Merged_Normalization"] ~> "Merged_Reduce_Dims"


--------------------------------------------------------------------------------
-- 
--------------------------------------------------------------------------------
{-
    node "Merged_Spectral_Sample" [| getSpectral 20000 |] $ return ()
    uNode "Merged_Spectral_Nystrom_Prep" [| \(a,b,c) -> 
        return $ zip3 a (repeat b) $ repeat c |]
    nodePar "Merged_Spectral_Nystrom" 'nystromExtend $ return ()
    ["Merged_Normalization", "Merged_Feature_Selection"] ~> "Merged_Spectral_Sample"
    ["Merged_Normalization", "Merged_Feature_Selection", "Merged_Spectral_Sample"] ~> "Merged_Spectral_Nystrom_Prep"
    path ["Merged_Spectral_Nystrom_Prep", "Merged_Spectral_Nystrom"]

    node "Merged_Reduce_Dims" 'outputReduced $ return ()
    ["Merge_Matrix", "Merged_Spectral_Nystrom"] ~> "Merged_Reduce_Dims"
    -}

    node "Merged_Batch_Correction" [| \case
        Nothing -> return Nothing
        Just input -> Just <$> batchCorrection "/Cluster/" input
        |] $ return ()
    node "Merged_Make_KNN" [| \case
        Nothing -> return Nothing
        Just input -> fmap Just $ mkKNNGraph "/Cluster/" $
            input & replicates.traverse.files %~ return
        |] $ nCore .= 4
    path ["Merged_Reduce_Dims", "Merged_Batch_Correction", "Merged_Make_KNN"]

--------------------------------------------------------------------------------
-- Selecting parameter
--------------------------------------------------------------------------------
    uNode "Merged_Param_Search_Prep" [| \case
        Just knn -> do
            ress <- asks _scrnaseq_cluster_resolution_list
            res <- maybe [] return <$> asks _scrnaseq_cluster_resolution
            optimizer <- asks _scrnaseq_cluster_optimizer 
            return $ flip map (nubSort $ res <> ress) $ \r ->
                ( optimizer, r
                , knn^.replicates._2.files._2.location )
        _ -> return []
        |]
    nodePar "Merged_Param_Search" [| \(optimizer, r, knn) -> do
        dir <- asks _scrnaseq_output_dir >>= getPath . (<> asDir ("/Cluster/Params/" <> show r))
        res <- liftIO $ evalClusters dir optimizer r knn
        return (r, res)
        |] $ return ()
    path ["Merged_Make_KNN", "Merged_Param_Search_Prep", "Merged_Param_Search"]


    uNode "Merged_Cluster_Metric_Prep" [| \(spec, x) -> case spec of
        Nothing -> return []
        Just fl -> return $ flip map x $ \(r, y) -> (r, y, fl^.replicates._2.files._2.location)
        |]
    nodePar "Merged_Cluster_Metric" [| \(res, cl, spec) -> do
        r <- liftIO $ computeClusterMetrics cl spec
        return (res, r)
        |] $ return ()
    node "Merged_Cluster" [| plotClusters |] $ return ()
    ["Merged_Reduce_Dims", "Merged_Param_Search"] ~> "Merged_Cluster_Metric_Prep"
    path ["Merged_Cluster_Metric_Prep", "Merged_Cluster_Metric"]
    ["Merged_Cluster_Metric", "Merged_Param_Search", "Merged_Make_KNN"] ~> "Merged_Cluster"

    node "Make_Cluster_Matrix" [| \case
        (Just mat, Just cl) -> fmap Just $
            segregateCells "/Quantification/Cluster/" mat $ cl^.replicates._2.files
        _ -> return Nothing
        |] $ return ()
    node "Make_Expr_Table" [| mkExprTable "/Quantification/" |] $ return ()
    ["Merge_Matrix", "Merged_Cluster"] ~> "Make_Cluster_Matrix"
    path ["Make_Cluster_Matrix", "Make_Expr_Table"]