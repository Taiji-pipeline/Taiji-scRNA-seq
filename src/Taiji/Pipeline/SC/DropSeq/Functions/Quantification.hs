{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
module Taiji.Pipeline.SC.DropSeq.Functions.Quantification
    (quantification) where

import           Bio.Data.Bam
import           Bio.Data.Bed
import           Bio.Utils.Misc (readInt)
import           Bio.Data.Bed.Types
import Data.Char (toLower)
import           Bio.Pipeline
import           Bio.Data.Experiment
import           Bio.RealWorld.GENCODE
import           Conduit
import           Data.Conduit.Internal (zipSinks)
import Data.List.Ordered
import           Control.Lens
import           Control.Monad
import           Control.Monad.IO.Class               (liftIO)
import           Control.Monad.Reader                 (asks)
import qualified Data.ByteString.Char8                as B
import           Data.CaseInsensitive                 (original)
import Data.Ord
import qualified Data.Map.Strict                  as M
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

type ExonTree = BEDTree [B.ByteString]

quantification :: DropSeqConfig config
               => RNASeq S (File '[NameSorted] 'Bam)
               -> WorkflowConfig config ( RNASeq S
                  ( File '[] 'Tsv
                  , [Double]
                  , M.Map Annotation Int) )
quantification input = do
    dir <- asks ((<> "/Quantification/") . _dropseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_quant.tsv" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    anno_f <- asks _dropseq_annotation
    exons <- liftIO $ mkExonTree <$> readGenes anno_f
    anno <- liftIO $ readAnnotations anno_f
    input & replicates.traverse.files %%~ liftIO . ( \bam -> do
        hdr <- getBamHeader $ bam^.location
        idxMap <- runResourceT $ runConduit $ streamBam (bam^.location) .|
            bamToBedC hdr .| getGeneIdxMap exons
        withFile output WriteMode $ \h -> do
            B.hPutStrLn h $ B.intercalate "\t" $ map fst $
                sortBy (comparing snd) $ M.toList idxMap 
            let quant = fmap snd $ mkGeneCount exons idxMap .| zipSinks
                    (mapC fst .| outputResult h) (mapC snd .| sinkList)
            (dupRate, annoCount) <- runResourceT $ runConduit $
                streamBam (bam^.location) .| bamToBedC hdr .|
                    groupOnC (fst . getIndex . fromJust . (^.name)) .|
                    zipSinks quant (concatC .| annotate anno)
            return (location .~ output $ emptyFile, dupRate, annoCount)
        )

getGeneIdxMap :: Monad m
              => ExonTree -> ConduitT BED o m (M.Map B.ByteString Int)
getGeneIdxMap exons = do
    s <- foldlC f S.empty
    return $ M.fromList $ zip (S.toList s) [0..]
  where
    f geneSet bed = foldl' (flip S.insert) geneSet $ concat $
        IM.elems $ intersecting exons bed

-- | Produce the gene counts
mkGeneCount :: Monad m
            => ExonTree
            -> M.Map B.ByteString Int
            -> ConduitT [BED] (U.Vector Int, Double) m ()
mkGeneCount exons idxMap = mapC $ \tags -> 
    let geneCount = M.fromListWith (++) $ concatMap fun tags
        geneCountUniq = fmap (length . nubSort) geneCount
        totalUMI = fromIntegral $ foldl' (+) 0 $ fmap length geneCount
        uniqUMI = fromIntegral $ foldl' (+) 0 geneCountUniq
        dupRate = 1 - uniqUMI / totalUMI
    in (sortCount geneCountUniq, dupRate)
  where
    fun bed = zip genes $ repeat [umi]
      where
        genes = concat $ IM.elems $ intersecting exons bed
        (_, umi) = getIndex $ fromJust $ bed^.name
    sortCount count = U.create $ do
        vec <- UM.replicate (M.size idxMap) 0
        forM_ (M.toList count) $ \(g, c) -> do
            let idx = M.findWithDefault undefined g idxMap
            UM.write vec idx c
        return vec

readAnnotations :: FilePath -> IO (BEDTree (S.HashSet Annotation))
readAnnotations input = do
    res <- fmap (bedToTree (++)) $ runResourceT $ runConduit $
        sourceFile input .| linesUnboundedAsciiC .|
        filterC (not . (=='#') . B.head) .| concatMapC f .| sinkList
    return $ (fmap . fmap) S.fromList res
  where
    f l = do
          a <- basicAnno
          return ( BED3 chr (readInt start - 1) (readInt end - 1)
              , a : (mito ++ rRNA) )
      where
        basicAnno = case B.map toLower name of
            "gene" -> Just Genic
            "cds" -> Just CDS
            "utr" -> Just UTR
            "exon" -> Just Exon
            _ -> Nothing
        mito | chr == "chrM" || chr == "M" = [Mitochondrial]
             | otherwise = []
        rRNA | ty == "rRNA" = [Ribosomal]
             | otherwise = []
        ty = B.init $ B.drop 2 $ fromJust $ lookup "gene_type" $
            map (B.break isSpace . strip) $ B.split ';' info
        [chr,_,name,start,end,_,_,_,info] = B.split '\t' l
        strip = fst . B.spanEnd isSpace . B.dropWhile isSpace
        isSpace = (== ' ')

annotate :: Monad m
         => BEDTree (S.HashSet Annotation)
         -> ConduitT BED o m (M.Map Annotation Int)
annotate anno = foldlC fun M.empty
  where
    fun m bed
        | Ribosomal `S.member` res = M.insertWith (+) Ribosomal 1 m
        | not (Genic `S.member` res) = M.insertWith (+) Intergenic 1 m 
        | not (Exon `S.member` res) = M.insertWith (+) Intron 1 m 
        | otherwise = foldl' (\acc x -> M.insertWith (+) x 1 acc) m res
      where
        res = S.unions $ IM.elems $ intersecting anno bed

outputResult :: MonadIO m => Handle -> ConduitT (U.Vector Int) o m ()
outputResult hdl = mapM_C $ liftIO . B.hPutStrLn hdl . B.intercalate "\t" .
    map (B.pack . show) . U.toList

mkExonTree :: [Gene] -> ExonTree
mkExonTree genes = bedToTree (++) $ concatMap f genes
  where
    f Gene{..} = map ( \(a,b) ->
        (asBed geneChrom a b :: BED3, [original geneName]) ) geneExon 


-- | Get barcode and UMI
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

            {-
groupOnC' :: (Monad m, Eq b)
          => (a -> b)
          -> ConduitT a (ConduitT i a m r) m ()
groupOnC' fun = await >>= (maybe (return ()) $ \b -> go (fun b) )
  where
    go idx acc = await >>= maybe (yield acc) f
      where
        f b | idx' == idx = go idx $ b : acc 
            | otherwise = yield acc >> go idx' [b]
          where
            idx' = fun b
            -}