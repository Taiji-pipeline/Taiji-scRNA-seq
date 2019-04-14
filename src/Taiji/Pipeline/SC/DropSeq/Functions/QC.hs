{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds         #-}

module Taiji.Pipeline.SC.DropSeq.Functions.QC
    (reportQC) where

import           Bio.Pipeline
import           Bio.Data.Experiment
import Control.Lens
import Data.List
import           Control.Monad.IO.Class               (liftIO)
import           Control.Monad.Reader                 (asks)
import           Scientific.Workflow

import qualified Data.Text as T
import Data.Aeson

import qualified Data.Map.Strict                  as M

import Taiji.Types
import           Taiji.Pipeline.SC.DropSeq.Types

reportQC :: DropSeqConfig config
         => [(RNASeq S (File '[] 'Tsv, [Double], M.Map Annotation Int))]
         -> WorkflowConfig config ()
reportQC x = do
    dir <- asks _dropseq_output_dir >>= getPath
    let output = dir ++ "/scRNA-seq.qc"
    liftIO $ encodeFile output $ [getDupRate x, getAnnotation x]

getDupRate :: [(RNASeq S (a, [Double], b))] -> QC
getDupRate es = QC "duplication_rate" $ Violin names dat
  where
    (names, dat) = unzip $ flip map es $ \e ->
        let name = T.unpack (e^.eid) <> "_rep" <> show (e^.replicates._1)
        in (name, e^.replicates._2.files._2)

getAnnotation :: [(RNASeq S (a, b, M.Map Annotation Int))] -> QC
getAnnotation es = 
    let (labels, values) = unzip dat
        dat' = zip labels $ transpose $ map normalize $ transpose values
    in QC "annotation" $ Bar names dat'
  where
    dat = [ ("CDS", map (M.findWithDefault 0 CDS) stats)
          , ("UTR", map (M.findWithDefault 0 UTR) stats)
          , ("Intron", map (M.findWithDefault 0 Intron) stats)
          , ("Intergenic", map (M.findWithDefault 0 Intergenic) stats)
          , ("rRNA", map (M.findWithDefault 0 Ribosomal) stats) ]
    (names, stats) = unzip $ flip map es $ \e ->
        let name = T.unpack (e^.eid) <> "_rep" <> show (e^.replicates._1)
        in (name, fmap fromIntegral $ e^.replicates._2.files._3)
    normalize xs = let n = foldl' (+) 0 xs in map (/n) xs