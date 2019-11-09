{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE QuasiQuotes #-}
module Taiji.Pipeline.SC.RNASeq.Functions.Quantification
    ( quantification
    , removeDoublet
    ) where

import           Bio.Data.Bam
import           Bio.Data.Bed
import           Bio.Utils.Misc (readInt)
import           Bio.Data.Bed.Types
import Data.Char (toLower)
import Language.Javascript.JMacro
import           Bio.Pipeline
import           Bio.RealWorld.GENCODE
import           Data.Conduit.Internal (zipSinks)
import Data.ByteString.Lex.Integral (packDecimal)
import Data.List.Ordered
import qualified Data.ByteString.Char8                as B
import           Data.CaseInsensitive                 (original)
import qualified Data.Map.Strict                  as M
import qualified Data.IntMap.Strict as I
import qualified Data.IntervalMap.Strict              as IM
import qualified Data.HashSet as S
import qualified Data.Text                            as T
import Shelly (shelly, run_)

import Taiji.Prelude hiding (_cell_barcode)
import Taiji.Utils
import           Taiji.Pipeline.SC.RNASeq.Types
import Taiji.Pipeline.SC.RNASeq.Functions.Utils
import Taiji.Utils.Plot
import Taiji.Utils.Plot.Vega

type ExonTree = BEDTree [Int]

passQC :: QC -> Bool
passQC QC{..} = _mitoRate <= 0.05 && _uniq_gene >= 200

mkExonTree :: [Gene] -> ExonTree
mkExonTree genes = bedToTree (++) $ concat $ zipWith f [0..] genes
  where
    f i Gene{..} = map ( \(a,b) -> (asBed geneChrom a b :: BED3, [i]) ) $
        nubSort $ concatMap transExon geneTranscripts

quantification :: SCRNASeqConfig config
               => RNASeq S (File '[NameSorted] 'Bam)
               -> ReaderT config IO ( RNASeq S
                  ( ( File '[] 'Tsv       
                    , File '[Gzip] 'Other )-- ^ quantification
                  , File '[] 'Tsv  ) )   -- ^ QC
quantification input = do
    dir <- asks ((<> "/Quantification/") . _scrnaseq_output_dir) >>= getPath
    dir2 <- qcDir
    let output = printf "%s/%s_rep%d.mat.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        idx = printf "%s/%s_rep%d_col_names.tsv" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        qc = printf "%s/%s_rep%d_qc.tsv" dir2 (T.unpack $ input^.eid)
            (input^.replicates._1)
    anno_f <- asks _scrnaseq_annotation

    input & replicates.traverse.files %%~ liftIO . ( \bam -> do
        hdr <- getBamHeader $ bam^.location

        genes <- readGenes anno_f
        anno <- readAnnotations anno_f

        B.writeFile idx $ B.unlines $ map (original . geneName) genes

        let exons = mkExonTree genes
            outputMat = filterC (passQC . snd) .| mapC fst .|
                sinkRows' (length genes) (fromJust . packDecimal) output
            outputQC = (yield qcFileHeader >> mapC (showQC . snd)) .|
                unlinesAsciiC .| sinkFile qc
        _ <- runResourceT $ runConduit $
            streamBam (bam^.location) .| bamToBedC hdr .|
            groupOnC (fst . getIndex . fromJust . (^.name)) .|
            mapMC ( \x -> do
                ((row, umi, dupRate), annoCount') <- runConduit $ x .|
                    zipSinks (getGeneCount exons) (annotate anno)
                let mitoRate = M.findWithDefault undefined Mitochondrial annoCount / fromIntegral umi
                    annoCount = fmap (((1 - dupRate)*) . fromIntegral) annoCount'
                return (row, QC (fst row) umi (length $ snd row) dupRate
                    mitoRate 0 annoCount)
            ) .| zipSinks outputMat outputQC
        return ( (location .~ idx $ emptyFile, location .~ output $ emptyFile)
               , location .~ qc $ emptyFile )
        )

removeDoublet :: SCRNASeqConfig config
              => RNASeq S ( ( File '[] 'Tsv, File '[Gzip] 'Other)
                            , File '[] 'Tsv  )
              -> ReaderT config IO ( RNASeq S
                 ( ( File '[] 'Tsv       
                   , File '[Gzip] 'Other ) 
                 , File '[] 'Tsv  ) )
removeDoublet input = do
    outdir <- asks ((<> "/Quantification/") . _scrnaseq_output_dir) >>= getPath
    dir <- qcDir
    let qcFile = dir <> "qc_with_dsc_" <> T.unpack (input^.eid) <> "_rep" <>
            show (input^.replicates._1) <> ".tsv"
        qcPlot = dir <> "doublet_" <> T.unpack (input^.eid) <> "_rep" <>
            show (input^.replicates._1) <> ".html"
        output = printf "%s/%s_rep%d_no_doublet.mat.gz" outdir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \((idx, matFl), qcFl) -> do
        doubletScore <- detectDoublet qcFile qcPlot (qcFl^.location) (matFl^.location)
        mat <- mkSpMatrix id (matFl^.location)
        let f x = M.findWithDefault 0 x doubletScore <= 0.5
            n = M.size $ M.filter (<=0.5) doubletScore
        runResourceT $ runConduit $
            streamRows mat .| filterC (f . fst) .| sinkRows n (_num_col mat) id output
        return ( (idx, location .~ output $ emptyFile)
               , location .~ qcFile $ emptyFile )
        )

detectDoublet :: FilePath
              -> FilePath
              -> FilePath
              -> FilePath
              -> IO (M.Map B.ByteString Double)
detectDoublet qcFile qcPlot oldQC matFl = withTemp Nothing $ \tmp -> do
    shelly $ run_ "taiji-utils" ["doublet", T.pack matFl, T.pack tmp]
    [probs, threshold, sc, sim_sc] <- B.lines <$> B.readFile tmp

    let th = readDouble threshold
        dProbs = map readDouble $ B.words probs
        ds = map readDouble $ B.words sc
        ds_sim = map readDouble $ B.words sim_sc
        rate = fromIntegral (length $ filter (>0.5) dProbs) /
            fromIntegral (length dProbs) * 100 :: Double
    savePlots qcPlot [ mkHist ds th <> title (printf "doublet percentage: %.1f%%" rate)
        , mkHist ds_sim th ] []

    mat <- mkSpMatrix id matFl
    bcs <- runResourceT $ runConduit $ streamRows mat .| mapC fst .| sinkList
    let doubletScore = M.fromList $ zip bcs dProbs

    qc <- readQC oldQC

    B.writeFile qcFile $ B.unlines $ (qcFileHeader:) $ flip map qc $ \s ->
        let v = M.findWithDefault 0 (_cell_barcode s) doubletScore
        in showQC s{_doubletScore = v}
    return doubletScore
  where
    mkHist xs ref = plt <> rule
      where
        plt = hist xs 100
        rule = option [jmacroE| {
            layer: {
                mark: "rule",
                data: {values: [{ref: `ref`}]},
                encoding: {
                    x: { field: "ref"},
                    color: {value: "red"},
                    size: {value: 1}
                }
            }
       } |]

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