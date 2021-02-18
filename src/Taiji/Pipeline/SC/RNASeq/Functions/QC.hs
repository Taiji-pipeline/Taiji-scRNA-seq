{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds         #-}

module Taiji.Pipeline.SC.RNASeq.Functions.QC
    (plotQC) where

import qualified Data.Text as T
import Control.Arrow (second)
import qualified Data.HashMap.Strict                  as M
import qualified Data.Vector as V

import Taiji.Prelude
import Taiji.Utils.Plot
import Taiji.Utils.Plot.Vega
import qualified Taiji.Utils.Plot.ECharts as ECharts
import           Taiji.Pipeline.SC.RNASeq.Types
import qualified Taiji.Utils.DataFrame as DF

plotQC :: SCRNASeqConfig config
       => [SCRNASeq S (a, File '[] 'Tsv)]
       -> ReaderT config IO ()
plotQC [] = return ()
plotQC inputs = do
    dir <- qcDir
    qc <- forM inputs $ \input -> liftIO $ do
        let name = (input^.eid) <> "_rep" <> T.pack (show (input^.replicates._1))
        qc <- readQC $ input^.replicates._2.files._2.location
        return (name, qc)
    liftIO $ savePlots (dir ++ "/QC.html")
        [ violin (map (extract (fromIntegral . _num_umi)) qc)
            defaultAxis{_axis_type="log", _axis_title="number of UMI"}
        , violin (map (extract (fromIntegral . _uniq_gene)) qc)
            defaultAxis{_axis_title="number of gene"}
        , violin (map (extract _dupRate) qc)
            defaultAxis{_axis_title="duplication rate"}
        ] [ECharts.stackBar $ mkStat qc]
  where
    extract f = second $ map f . filter passQC
    mkStat :: [(T.Text, [QC])] -> DF.DataFrame Double
    mkStat qc = DF.transpose $ DF.mapRows normalize $
        DF.mkDataFrame (map fst qc) (map (T.pack . show) feats) $
        map (f . foldl' (M.unionWith (+)) M.empty . map _count_table . snd) qc
      where
        feats = [Exon, Intron, Intergenic, Ribosomal, Mitochondrial]
        f m = map (\x -> fromIntegral $ M.lookupDefault 0 x m) feats
        normalize xs = let s = V.foldl1 (+) xs in V.map (/s) xs