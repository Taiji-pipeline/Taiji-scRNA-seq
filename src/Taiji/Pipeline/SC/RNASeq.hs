{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.SC.RNASeq ( builder) where

import           Control.Workflow

import Taiji.Prelude
import           Taiji.Pipeline.SC.RNASeq.Functions

builder :: Builder ()
builder = do
    node "Read_Input" 'readInput $ return ()

    uNode "Get_Fastq" 'getFastq
    nodePar "Demultiplex" 'extractBarcode $ return ()
    path ["Read_Input", "Get_Fastq", "Demultiplex"]

    uNode "Get_Demulti_Fastq" [| \(input, fq) -> getDemultiplexedFastq input ++ fq |]
    ["Read_Input", "Demultiplex"] ~> "Get_Demulti_Fastq"

    node "Make_Index" 'mkIndex $ return ()
    nodePar "Align" 'tagAlign $ do
        nCore .= 8
        memory .= 50
    nodePar "Filter_Bam" 'filterNameSortBam $ return ()
    nodePar "Quantification" 'quantification $ memory .= 8
    nodePar "Filter_Cell" 'filterCells $ return ()
    nodePar "Remove_Doublet" 'removeDoublet $ return ()
    path ["Get_Demulti_Fastq", "Make_Index", "Align", "Filter_Bam"
        , "Quantification", "Filter_Cell", "Remove_Doublet"]

    node "QC" 'plotQC $ return ()
    ["Quantification"] ~> "QC"

    node "Merge_Matrix" [| \(inputs, mats) -> do
        let a = getMatrix inputs & mapped.replicates._2.files %~ Left
            b = mats & mapped.replicates._2.files %~ Right . fst
        combineMatrices $ a ++ b
        |] $ return ()
    ["Read_Input", "Remove_Doublet"] ~> "Merge_Matrix"

    {-
    node "Feature_Selection" [| \case
        Nothing -> return Nothing
        Just input -> do
            let prefix = "/Cluster/"
            idx <- selectFeature' prefix input
            return $ Just (idx, input)
            --return $ Just ([], input)
        |] $ return ()
    -}
    node "Merged_Reduce_Dimension" [| \case
        Nothing -> return Nothing
        Just input -> do
            let prefix = "/Cluster/"
            fmap Just $ filterMatrix prefix input >>= spectral prefix Nothing
        |] $ return ()
    node "Merged_Make_KNN" [| \case
        Nothing -> return Nothing
        Just input -> fmap Just $ mkKNNGraph "/Cluster/" $
            input & replicates.traverse.files %~ return
        |] $ return ()
    node "Merged_Cluster" [| \case
        Nothing -> return Nothing
        Just input -> Just <$> clustering "/Cluster/" input
        |] $ return ()
    path ["Merge_Matrix", "Merged_Reduce_Dimension", "Merged_Make_KNN", "Merged_Cluster"]

    node "Make_Cluster_Matrix" [| \case
        (Just mat, Just cl) -> segregateCells "/Quantification/Cluster/" mat $
            cl^.replicates._2.files
        _ -> return []
        |] $ return ()
    node "Make_Expr_Table" [| mkExprTable "/Quantification/" |] $ return ()
    ["Merge_Matrix", "Merged_Cluster"] ~> "Make_Cluster_Matrix"
    path ["Make_Cluster_Matrix", "Make_Expr_Table"]