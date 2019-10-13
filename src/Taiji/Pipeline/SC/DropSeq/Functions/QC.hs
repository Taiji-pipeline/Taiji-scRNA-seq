{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds         #-}

module Taiji.Pipeline.SC.DropSeq.Functions.QC
    (plotQC) where

import qualified Data.Text as T
import qualified Data.Matrix as M
import qualified Data.Vector as V

import Taiji.Prelude
import Taiji.Utils.Plot
import Taiji.Utils.Plot.Vega
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

plotQC :: DropSeqConfig config
       => [RNASeq S (File '[] 'Tsv, File '[] 'Tsv)]
       -> ReaderT config IO ()
plotQC [] = return ()
plotQC inputs = do
    dir <- qcDir
    forM_ inputs $ \input -> liftIO $ do
        let output = printf "%s/%s_rep%d_qc.html" dir (T.unpack $ input^.eid)
                (input^.replicates._1)
        df <- DF.readTable $ input^.replicates._2.files._2.location
        let plts = flip map (M.toColumns $ DF._dataframe_data df) $ \col -> hist (V.toList col) 100
        savePlots output plts []