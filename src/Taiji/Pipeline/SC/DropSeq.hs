{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.SC.DropSeq
    ( builder
    , DropSeqConfig(..)
    ) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser
import           Control.Lens
import           Control.Monad.IO.Class                  (liftIO)
import           Control.Monad.Reader                    (asks)
import           Data.Either
import           Scientific.Workflow

import           Taiji.Pipeline.SC.DropSeq.Config
import           Taiji.Pipeline.SC.DropSeq.Functions
import           Taiji.Pipeline.RNASeq.Functions (rnaGetFastq)

builder :: Builder ()
builder = do
    nodeS "Read_Input" [| \_ -> do
        input <- asks _dropSeq_input
        liftIO $ if ".tsv" == reverse (take 4 $ reverse input)
            then readRNASeqTSV input "Drop-seq"
            else readRNASeq input "Drop-seq"
        |] $ submitToRemote .= Just False

    node' "Get_Fastq" [|
        map (\x -> x & replicates.mapped.files %~ fromRight undefined) .
        filter (\x -> isRight $ x^.replicates._2.files ) .
        rnaGetFastq
        |] $ submitToRemote .= Just False
    nodePS 1 "Extract_Barcode" 'extractBarcode $ remoteParam .= "--mem=20000 -p gpu"
    path ["Read_Input", "Get_Fastq", "Extract_Barcode"]

    nodeS "Make_Index" 'mkIndex $ remoteParam .= "--mem=40000 -p gpu"
    nodePS 1 "Align" 'tagAlign $ remoteParam .= "--ntasks-per-node=4 --mem=40000 -p gpu"
    nodePS 1 "Filter_Bam" 'filterBamFile $ return ()
    nodePS 1 "Quantification" 'quantification $ return ()
    path ["Extract_Barcode", "Make_Index", "Align", "Filter_Bam", "Quantification"]

