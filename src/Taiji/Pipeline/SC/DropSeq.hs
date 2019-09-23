{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.SC.DropSeq ( builder) where

import           Control.Workflow

import Taiji.Prelude
import           Taiji.Pipeline.SC.DropSeq.Functions

builder :: Builder ()
builder = do
    node "Read_Input" 'readInput $ return ()

    node "Get_Fastq" [| return . getFastq |] $ return ()
    nodePar "Extract_Barcode" 'extractBarcode $ return ()
    path ["Read_Input", "Get_Fastq", "Extract_Barcode"]

    node "Make_Index" 'mkIndex $ return ()
    nodePar "Align" 'tagAlign $ nCore .= 4
    nodePar "Filter_Bam" 'filterBamSort $ return ()
    nodePar "Quantification" 'quantification $ return ()
    path ["Extract_Barcode", "Make_Index", "Align", "Filter_Bam", "Quantification"]

    node "QC_Anno" 'annoQC $ return ()
    ["Quantification"] ~> "QC_Anno"
