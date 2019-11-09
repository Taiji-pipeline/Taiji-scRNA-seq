{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.SC.RNASeq ( builder) where

import           Control.Workflow

import Taiji.Prelude
import           Taiji.Pipeline.SC.RNASeq.Functions

builder :: Builder ()
builder = do
    node "Read_Input" 'readInput $ return ()

    node "Get_Fastq" [| return . getFastq |] $ return ()
    nodePar "Extract_Barcode" 'extractBarcode $ return ()
    path ["Read_Input", "Get_Fastq", "Extract_Barcode"]

    node "Make_Index" 'mkIndex $ return ()
    nodePar "Align" 'tagAlign $ do
        nCore .= 8
        memory .= 50
    nodePar "Filter_Bam" 'filterNameSortBam $ return ()
    nodePar "Quantification" 'quantification $ return ()
    nodePar "Remove_Doublet" 'removeDoublet $ return ()
    path ["Extract_Barcode", "Make_Index", "Align", "Filter_Bam"
        , "Quantification", "Remove_Doublet"]

    nodePar "Reduce_Dimension" [| \input -> do
        let input' = input & replicates.traverse.files %~ snd . fst
            prefix = "/Cluster/"
        filterMatrix prefix input' >>= spectral "/Cluster/" Nothing 
        |] $ return ()
    path ["Remove_Doublet", "Reduce_Dimension"] 

    node "QC" 'plotQC $ return ()
    ["Quantification"] ~> "QC"
