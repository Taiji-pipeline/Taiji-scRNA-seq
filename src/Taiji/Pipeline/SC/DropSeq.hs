{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.SC.DropSeq ( builder) where

import           Scientific.Workflow
import Control.Lens

import           Taiji.Pipeline.SC.DropSeq.Functions

builder :: Builder ()
builder = do
    nodeS "Read_Input" 'readInput $ submitToRemote .= Just False

    node' "Get_Fastq" 'getFastq $ submitToRemote .= Just False
    nodePS 1 "Extract_Barcode" 'extractBarcode $ remoteParam .= "--mem=20000 -p gpu"
    path ["Read_Input", "Get_Fastq", "Extract_Barcode"]

    nodeS "Make_Index" 'mkIndex $ remoteParam .= "--mem=40000 -p gpu"
    nodePS 1 "Align" 'tagAlign $ remoteParam .= "--ntasks-per-node=4 --mem=40000 -p gpu"
    nodePS 1 "Filter_Bam" 'filterBamSort $ return ()
    nodePS 1 "Quantification" 'quantification $ return ()
    path ["Extract_Barcode", "Make_Index", "Align", "Filter_Bam", "Quantification"]

    nodeS "Run_QC" 'reportQC $ submitToRemote .= Just False
    ["Quantification"] ~> "Run_QC"
