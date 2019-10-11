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
import           Bio.RealWorld.GENCODE
import           Data.Conduit.Internal (zipSinks)
import Data.ByteString.Lex.Integral (packDecimal)
import qualified Data.Vector.Unboxed as U
import Data.List.Ordered
import qualified Data.ByteString.Char8                as B
import           Data.CaseInsensitive                 (original)
import qualified Data.Map.Strict                  as M
import qualified Data.IntMap.Strict as I
import qualified Data.IntervalMap.Strict              as IM
import qualified Data.HashSet as S
import qualified Data.Text                            as T

import Taiji.Prelude
import Taiji.Utils
import           Taiji.Pipeline.SC.DropSeq.Types
import Taiji.Pipeline.SC.DropSeq.Functions.Utils

type ExonTree = BEDTree [Int]

mkExonTree :: [Gene] -> ExonTree
mkExonTree genes = bedToTree (++) $ concat $ zipWith f [0..] genes
  where
    f i Gene{..} = map ( \(a,b) -> (asBed geneChrom a b :: BED3, [i]) ) $
        nubSort $ concatMap transExon geneTranscripts

data QC = QC
    { _cell_barcode :: B.ByteString
    , _dupRate :: Double
    , _num_umi :: Int
    , _genomics_context :: M.Map Annotation Int }

showQC :: QC -> B.ByteString
showQC QC{..} = B.intercalate "\t" $ (_cell_barcode:) $ map (B.pack . show) $
    _dupRate : map fromIntegral (_num_umi : dat)
  where
    dat = map (\x -> M.findWithDefault 0 x _genomics_context)
        [CDS, UTR, Intron, Intergenic, Ribosomal, Mitochondrial]

quantification :: DropSeqConfig config
               => RNASeq S (File '[NameSorted] 'Bam)
               -> ReaderT config IO ( RNASeq S
                  ( File '[] 'Tsv
                  , File '[] 'Tsv ) )
quantification input = do
    dir <- asks ((<> "/Quantification/") . _dropseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_quant.tsv" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        qc = printf "%s/%s_rep%d_qc.tsv" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    anno_f <- asks _dropseq_annotation

    input & replicates.traverse.files %%~ liftIO . ( \bam -> do
        hdr <- getBamHeader $ bam^.location

        genes <- readGenes anno_f
        anno <- readAnnotations anno_f
        let exons = mkExonTree genes
            geneNames = B.intercalate "\t" $ "" : map (original . geneName) genes
            printRow (nm, val) = B.intercalate "\t" $ nm :
                map (fromJust . packDecimal)
                (U.toList $ U.replicate (length genes) 0 U.// val)
            outputMat = mapC fst .| (yield geneNames >> mapC printRow) .|
                unlinesAsciiC .| sinkFile output
            outputQC = (yield header >> mapC (showQC . snd)) .|
                unlinesAsciiC .| sinkFile qc
            header = "\tduplication\tUMI\tCDS\tUTR\tIntron\tIntergenic\tRibosomal\tMitochondrial"

        _ <- runResourceT $ runConduit $
            streamBam (bam^.location) .| bamToBedC hdr .|
            groupOnC (fst . getIndex . fromJust . (^.name)) .|
            mapMC ( \x -> do
                ((row, umi, dupRate), annoCount) <- runConduit $ x .|
                    zipSinks (getGeneCount exons) (annotate anno)
                return (row, QC (fst row) dupRate umi annoCount)
            ) .| zipSinks outputMat outputQC
        return ( location .~ output $ emptyFile
               , location .~ qc $ emptyFile )
        )

{-
mkTable :: DropSeqConfig config
        => RNASeq S (File '[] 'Other , File '[] 'Tsv)
        -> RNASeq S (File '[] 'Tsv)
mkTable input = do
    dir <- asks ((<> "/Quantification/") . _dropseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_quant.tsv" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    anno_f <- asks _dropseq_annotation

    input & replicates.traverse.files %%~ liftIO . ( \(fl, _) -> do
        anno <- readAnnotations anno_f
        -}

getGeneCount :: Monad m
             => ExonTree
             -> ConduitT BED o m (Row Int, Int, Double)
getGeneCount exons = do
    cellBc <- fst . getIndex . fromJust . (^.name) . fromJust <$> peekC

    geneCount <- concatMapC count .| foldlC reduce I.empty

    let c = I.toList $ fmap M.size geneCount
        uniqUMI = foldl' (+) 0 $ map snd c
        totalUMI = fromIntegral $ foldl' (+) 0 $ fmap (foldl' (+) 0) geneCount
        dupRate = 1 - fromIntegral uniqUMI / totalUMI

    return ((cellBc, c), uniqUMI, dupRate)
  where
    reduce m (idx, umi) = I.alter f idx m
      where
        f Nothing = Just $ M.singleton umi (1 :: Int)
        f (Just x) = Just $ M.insertWith (+) umi 1 x
    count :: BED -> [(Int, B.ByteString)]
    count bed = zip (concat $ IM.elems $ intersecting exons bed) $ repeat $
        snd $ getIndex $ fromJust $ bed^.name


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
        basicAnno = case B.map toLower nm of
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
        [chr,_,nm,start,end,_,_,_,info] = B.split '\t' l
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