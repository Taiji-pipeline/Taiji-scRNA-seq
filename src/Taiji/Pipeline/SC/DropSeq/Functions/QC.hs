{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds         #-}

module Taiji.Pipeline.SC.DropSeq.Functions.QC
    (reportQC) where

import Data.Aeson
import           Bio.Pipeline
import           Bio.Data.Experiment
import Control.Lens
import           Control.Monad.IO.Class               (liftIO)
import           Control.Monad.Reader                 (asks)
import           Scientific.Workflow

import qualified Data.Text as T

import Taiji.Types
import           Taiji.Pipeline.SC.DropSeq.Types

reportQC :: DropSeqConfig config
         => [(RNASeq S (File '[] 'Tsv, [Double]))]
         -> WorkflowConfig config ()
reportQC x = do
    dir <- asks _dropseq_output_dir >>= getPath
    let output = dir ++ "/scRNA-seq.qc"
    liftIO $ encodeFile output $ [getDupRate x]

getDupRate :: [(RNASeq S (File '[] 'Tsv, [Double]))] -> QC
getDupRate es = QC "duplication_rate" res Violin
  where
    res = QCResult names rates
    (names, rates) = unzip $ flip map es $ \e ->
        let name = toJSON $ T.unpack (e^.eid) <> "_rep" <> show (e^.replicates._1)
        in (name, toJSON $ e^.replicates._2.files._2)
