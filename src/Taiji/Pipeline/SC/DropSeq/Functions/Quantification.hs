{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
module Taiji.Pipeline.SC.DropSeq.Functions.Quantification
    (quantification) where

import           Bio.Data.Bam
import           Bio.Data.Bed
import           Bio.Pipeline
import           Bio.Data.Experiment
import           Bio.RealWorld.GENCODE
import           Conduit
import           Data.Conduit.Internal (zipSinks)
import           Control.Lens
import           Control.Monad
import           Control.Monad.IO.Class               (liftIO)
import           Control.Monad.Reader                 (asks)
import qualified Data.ByteString.Char8                as B
import           Data.CaseInsensitive                 (original)
import Data.Ord
import qualified Data.HashMap.Strict                  as M
import qualified Data.IntervalMap.Strict              as IM
import qualified Data.HashSet as S
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import           Data.List
import           Data.Maybe
import qualified Data.Text                            as T
import System.IO
import           Scientific.Workflow
import           Text.Printf                          (printf)

import           Taiji.Pipeline.SC.DropSeq.Types

type DuplicationRates = [Double]

type ExonTree = BEDTree [B.ByteString]

quantification :: DropSeqConfig config
               => RNASeq S (File '[NameSorted] 'Bam)
               -> WorkflowConfig config (RNASeq S (File '[] 'Tsv, DuplicationRates))
quantification input = do
    dir <- asks ((<> "/Quantification/") . _dropseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_quant.tsv" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    anno_f <- asks _dropseq_annotation
    annotation <- liftIO $ mkExonTree <$> readGenes anno_f
    input & replicates.traverse.files %%~ liftIO . ( \bam -> do
        hdr <- getBamHeader $ bam^.location
        idxMap <- runResourceT $ runConduit $ streamBam (bam^.location) .|
            bamToBedC hdr .| getGeneIdxMap annotation
        withFile output WriteMode $ \h -> do
            B.hPutStrLn h $ B.intercalate "\t" $ map fst $
                sortBy (comparing snd) $ M.toList idxMap 
            dupRate <- runResourceT $ runConduit $ streamBam (bam^.location) .|
                bamToBedC hdr .| groupOnC (fst . getIndex . fromJust . (^.name)) .|
                quantify h annotation idxMap
            return (location .~ output $ emptyFile, dupRate)
        )

getGeneIdxMap :: Monad m
              => ExonTree -> ConduitT BED o m (M.HashMap B.ByteString Int)
getGeneIdxMap exons = do
    s <- foldlC f S.empty
    return $ M.fromList $ zip (S.toList s) [0..]
  where
    f geneSet bed = foldl' (flip S.insert) geneSet $ concat $
        IM.elems $ intersecting exons bed

quantify :: MonadIO m 
         => Handle  -- ^ output
         -> ExonTree
         -> M.HashMap B.ByteString Int
         -> ConduitT [BED] Void m DuplicationRates
quantify output exon idxMap = fmap snd $ mkGeneCount exon idxMap .|
    zipSinks (mapC fst .| outputResult output) (mapC snd .| sinkList)

-- | Produce the gene counts and duplication rate.
mkGeneCount :: Monad m
            => ExonTree
            -> M.HashMap B.ByteString Int
            -> ConduitT [BED] (U.Vector Int, Double) m ()
mkGeneCount exons idxMap = mapC $ \beds ->
    let totalCount = fromIntegral $ foldl' (+) 0 geneCount
        geneCount = fmap S.size $ foldl' fun M.empty beds
    in (sortCount geneCount, 1 - totalCount / fromIntegral (length beds))
  where
    fun table bed = foldl' f table genes
      where
        f m g = M.insertWith (<>) g (S.singleton umi) m
        genes = concat $ IM.elems $ intersecting exons bed
        (_, umi) = getIndex $ fromJust $ bed^.name
    sortCount count = U.create $ do
        vec <- UM.replicate (M.size idxMap) 0
        forM_ (M.toList count) $ \(g, c) -> do
            let idx = M.lookupDefault undefined g idxMap
            UM.unsafeWrite vec idx c
        return vec


outputResult :: MonadIO m => Handle -> ConduitT (U.Vector Int) o m ()
outputResult hdl = mapM_C $ liftIO . B.hPutStrLn hdl . B.intercalate "\t" .
    map (B.pack . show) . U.toList

mkExonTree :: [Gene] -> ExonTree
mkExonTree genes = bedToTree (++) $ concatMap f genes
  where
    f Gene{..} = map ( \(a,b) ->
        (asBed geneChrom a b :: BED3, [original geneName]) ) geneExon 


getIndex :: B.ByteString -> (B.ByteString, B.ByteString)
getIndex x = (a,b)
  where
    [a,b] = B.split '+' $ head $ B.split '_' x
{-# INLINE getIndex #-}

groupOnC :: (Monad m, Eq b)
         => (a -> b)
         -> ConduitT a [a] m ()
groupOnC fun = await >>= (maybe (return ()) $ \b -> go (fun b) [b])
  where
    go idx acc = await >>= maybe (yield acc) f
      where
        f b | idx' == idx = go idx $ b : acc 
            | otherwise = yield acc >> go idx' [b]
          where
            idx' = fun b