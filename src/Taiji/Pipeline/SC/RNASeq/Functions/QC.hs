{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds         #-}

module Taiji.Pipeline.SC.RNASeq.Functions.QC
    (plotQC) where

import qualified Data.Text as T
import Control.Arrow (second)

import Taiji.Prelude
import Taiji.Utils.Plot
import Taiji.Utils.Plot.Vega
import           Taiji.Pipeline.SC.RNASeq.Types

plotQC :: SCRNASeqConfig config
       => [RNASeq S (a, File '[] 'Tsv)]
       -> ReaderT config IO ()
plotQC [] = return ()
plotQC inputs = do
    dir <- qcDir
    qc <- forM inputs $ \input -> liftIO $ do
        let name = (input^.eid) <> "_rep" <> T.pack (show (input^.replicates._1))
        qc <- fmap (filter passQC) $ readQC $ input^.replicates._2.files._2.location
        return (name, qc)
    liftIO $ savePlots (dir ++ "/QC.html")
        [ violin (map (second $ map $ fromIntegral . _num_umi) qc)
            defaultAxis{_axis_type="log", _axis_title="number of UMI"}
        , violin (map (second $ map $ fromIntegral . _uniq_gene) qc)
            defaultAxis{_axis_title="number of gene"}
        , violin (map (second $ map _dupRate) qc)
            defaultAxis{_axis_title="duplication rate"}
        ] []