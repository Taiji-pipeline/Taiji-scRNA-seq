{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds         #-}

module Taiji.Pipeline.SC.DropSeq.Functions.QC
    (annoQC) where

import qualified Data.Text as T
import qualified Data.Map.Strict                  as M

import Taiji.Prelude
import Taiji.Utils.Plot
import Taiji.Utils.Plot.ECharts
import qualified Taiji.Utils.DataFrame as DF
import           Taiji.Pipeline.SC.DropSeq.Types

{-
dupRate :: DropSeqConfig config
        => [(RNASeq S (a, [Double], b))]
        -> ReaderT config IO ()
dupRate es = do
    dir <- qcDir
    let output = dir ++ "qc_duplication.html"
    
    vegaViolin $ flip map es $ \e ->
    let name = (e^.eid) <> "_rep" <> T.pack (show (e^.replicates._1))
    in (name, e^.replicates._2.files._2)
    -}

annoQC :: DropSeqConfig config
       => [(RNASeq S (a, b, M.Map Annotation Int))]
       -> ReaderT config IO ()
annoQC es = do
    dir <- qcDir
    let output = dir ++ "annotation.html"
        (labels, values) = unzip dat
        df = DF.mkDataFrame names labels $ transpose $ map normalize $
            transpose values
    liftIO $ savePlots output [] [stackBar df]
  where
    dat = [ ("CDS", map (M.findWithDefault 0 CDS) stats)
          , ("UTR", map (M.findWithDefault 0 UTR) stats)
          , ("Intron", map (M.findWithDefault 0 Intron) stats)
          , ("Intergenic", map (M.findWithDefault 0 Intergenic) stats)
          , ("rRNA", map (M.findWithDefault 0 Ribosomal) stats) ]
    (names, stats) = unzip $ flip map es $ \e ->
        let name = (e^.eid) <> "_rep" <> T.pack (show (e^.replicates._1))
        in (name, fmap fromIntegral $ e^.replicates._2.files._3)
    normalize xs = let n = foldl' (+) 0 xs in map ((100*) . (/n)) xs