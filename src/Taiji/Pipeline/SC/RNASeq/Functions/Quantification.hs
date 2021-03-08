{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE QuasiQuotes #-}
module Taiji.Pipeline.SC.RNASeq.Functions.Quantification
    ( quantification
    , filterCells
    , removeDoublet
    ) where

import           Bio.Data.Bam
import           Bio.Data.Bed
import Data.Conduit.Zlib (gzip)
import Data.Char (toLower)
import Language.Javascript.JMacro
import           Bio.Pipeline
import           Bio.RealWorld.GENCODE
import           Data.Conduit.Internal (zipSinks, zipSources)
import Data.List.Ordered
import qualified Data.ByteString.Char8                as B
import           Data.CaseInsensitive                 (original)
import qualified Data.HashMap.Strict                  as M
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

mkExonTree :: [Gene] -> ExonTree
mkExonTree genes = bedToTree (++) $ concat $ zipWith f [0..] genes
  where
    f i Gene{..} = zip (runIdentity $ runConduit $ mergeBed exon .| sinkList) $
        repeat [i]
      where
        exon = map (\(a,b) -> asBed geneChrom a b :: BED3) $
            concatMap transExon geneTranscripts
{-# INLINE mkExonTree #-}

quantification :: SCRNASeqConfig config
               => SCRNASeq S (File '[NameSorted] 'Bam)
               -> ReaderT config IO ( SCRNASeq S
                  ( ( File '[ColumnName, Gzip] 'Tsv       
                    , File '[Gzip] 'Other ) -- ^ quantification
                    , File '[] 'Tsv  ) )   -- ^ QC
quantification input = do
    dir <- asks ((<> "/Quantification/") . _scrnaseq_output_dir) >>= getPath
    dir2 <- tempDir
    let output = printf "%s/%s_rep%d.mat.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        idx = printf "%s/%s_rep%d_col_names.tsv.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        qc = printf "%s/%s_rep%d_qc.tsv" dir2 (T.unpack $ input^.eid)
            (input^.replicates._1)
    anno_f <- asks _scrnaseq_annotation

    input & replicates.traverse.files %%~ liftIO . ( \bam -> do
        hdr <- getBamHeader $ bam^.location

        genes <- readGenesValidated anno_f
        anno <- readAnnotations anno_f

        runResourceT $ runConduit $ yieldMany genes .|
            mapC (original . geneName) .| unlinesAsciiC .| gzip .| sinkFile idx

        let exons = mkExonTree genes
            outputMat = mapC fst .|
                sinkRows' (length genes) (fromJust . packDecimal) output
            outputQC = (yield qcFileHeader >> mapC (showQC . snd)) .|
                unlinesAsciiC .| sinkFile qc
        _ <- runResourceT $ runConduit $
            streamBam (bam^.location) .| bamToBedC hdr .|
            groupOnC (fst . getIndex . fromJust . (^.name)) .|
            mapMC (quantify exons anno) .| zipSinks outputMat outputQC
        return ( (location .~ idx $ emptyFile, location .~ output $ emptyFile)
               , location .~ qc $ emptyFile )
        )

filterCells :: SCRNASeqConfig config
            => SCRNASeq S
                  ( ( File '[ColumnName, Gzip] 'Tsv       
                    , File '[Gzip] 'Other )-- ^ quantification
                    , File '[] 'Tsv  )  -- ^ QC
            -> ReaderT config IO ( SCRNASeq S
                  ( ( File '[ColumnName, Gzip] 'Tsv       
                    , File '[Gzip] 'Other )
                    , File '[] 'Tsv  ) )
filterCells input = do
    dir <- asks ((<> "/Quantification/") . _scrnaseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_filt.mat.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \((idx, matFl), qcFl) -> do
        qc <- readQC $ qcFl^.location
        mat <- mkSpMatrix id (matFl^.location)
        let n = length $ filter passQC qc
        runResourceT $ runConduit $
            zipSources (streamRows mat) (yieldMany qc) .|
            filterC (passQC . snd) .| mapC fst .|
            sinkRows n (_num_col mat) id output
        return ( (idx, location .~ output $ emptyFile)
               , qcFl )
        )

quantify :: Monad m
         => ExonTree
         -> BEDTree (S.HashSet Annotation)
         -> ConduitT () BED m ()
         -> m (Row Int, QC)
quantify exons anno input = do
    (n, (bc, exonCount), counts) <- runConduit $ input .| sink
    let nUMI = foldl' (+) 0 counts
        mito =
            let s = nUMI - M.lookupDefault 0 Intergenic counts -
                    M.lookupDefault 0 Ribosomal counts
            in if s == 0
                then 0
                else fromIntegral (M.lookupDefault 0 Mitochondrial counts) /
                    fromIntegral s
        qc = QC
            { _cell_barcode = bc
            , _num_umi = nUMI
            , _uniq_gene = length exonCount
            , _dupRate = 1 - fromIntegral nUMI / fromIntegral (n :: Int)
            , _mitoRate = mito
            , _doubletScore = 0
            , _count_table = counts }
    return ((bc, exonCount), qc)
  where
    sink = getZipSink $ (,,) <$> ZipSink lengthC <*>
        ZipSink (countExon exons) <*> ZipSink (countFeat anno)
{-# INLINE quantify #-}

countExon :: Monad m
          => ExonTree
          -> ConduitT BED o m (Row Int)
countExon exons = do
    cellBc <- fst . getIndex . fromJust . (^.name) . fromJust <$> peekC
    geneCount <- concatMapC count .| foldlC reduce I.empty
    return (cellBc, I.toList $ fmap M.size geneCount)
  where
    reduce m (idx, umi) = I.alter f idx m
      where
        f Nothing = Just $ M.singleton umi (1 :: Int)
        f (Just x) = Just $ M.insertWith (+) umi 1 x
    count :: BED -> [(Int, B.ByteString)]
    count bed = zip (concat $ IM.elems $ intersecting exons bed) $ repeat $
        snd $ getIndex $ fromJust $ bed^.name
{-# INLINE countExon #-}

readAnnotations :: FilePath -> IO (BEDTree (S.HashSet Annotation))
readAnnotations input = do
    res <- fmap (bedToTree (++) . concatMap mergeRegion . groupBy ((==) `on` fst) . sortBy (comparing fst)) $
        runResourceT $ runConduit $
        sourceFile input .| linesUnboundedAsciiC .|
        filterC (not . (=='#') . B.head) .| concatMapC f .| concatC .| sinkList
    return $ (fmap . fmap) S.fromList res
  where
    mergeRegion xs =
        let (a,b) = unzip xs
        in zip (runIdentity $ runConduit $ mergeBed b .| sinkList) $ map return a
    f l = do
          a <- basicAnno
          return $ zip (a : mito ++ rRNA) $ repeat region
      where
        region = BED3 chr (readInt start - 1) $ readInt end - 1
        basicAnno = case B.map toLower nm of
            "gene" -> Just Genic
            "exon" -> Just Exon
            _ -> Nothing
        mito | chr == "chrM" || chr == "M" = [Mitochondrial]
             | otherwise = []
        rRNA | ty == "rRNA" = [Ribosomal]
             | otherwise = []
        ty = B.init $ B.drop 2 $ fromJust $ lookup "gene_type" $
            map (B.break isSpace . strip) $ B.split ';' anno
        [chr,_,nm,start,end,_,_,_,anno] = B.split '\t' l
        strip = fst . B.spanEnd isSpace . B.dropWhile isSpace
        isSpace = (== ' ')
{-# INLINE readAnnotations #-}

countFeat :: Monad m
          => BEDTree (S.HashSet Annotation)
          -> ConduitT BED o m (M.HashMap Annotation Int)
countFeat anno = fmap S.size <$> foldlC go M.empty
  where
    go accum bed
        | null feats = addUMI Intergenic accum
        | Ribosomal `S.member` feats = addUMI Ribosomal accum
        | Mitochondrial `S.member` feats = addUMI Mitochondrial accum
        | Exon `S.member` feats = addUMI Exon accum
        | otherwise = addUMI Intron accum
      where
        feats = S.unions $ IM.elems $ intersecting anno bed
        addUMI :: Annotation 
               -> M.HashMap Annotation (S.HashSet B.ByteString)
               -> M.HashMap Annotation (S.HashSet B.ByteString)
        addUMI = M.alter f
          where
            f Nothing = Just $ S.singleton umi 
            f (Just s) = Just $ S.insert umi s
            umi = snd $ getIndex $ fromJust $ bed^.name
{-# INLINE countFeat #-}


removeDoublet :: SCRNASeqConfig config
              => SCRNASeq S ( ( File '[ColumnName, Gzip] 'Tsv, File '[Gzip] 'Other)
                            , File '[] 'Tsv  )
              -> ReaderT config IO ( SCRNASeq S
                 ( ( File '[ColumnName, Gzip] 'Tsv       
                   , File '[Gzip] 'Other ) 
                 , File '[] 'Tsv  ) )
removeDoublet input = do
    outdir <- asks ((<> "/Quantification/") . _scrnaseq_output_dir) >>= getPath
    dir <- qcDir
    thres <- asks _scrnaseq_doublet_score_cutoff 
    let qcFile = dir <> "qc_" <> T.unpack (input^.eid) <> "_rep" <>
            show (input^.replicates._1) <> ".tsv"
        qcPlot = dir <> "doublet_" <> T.unpack (input^.eid) <> "_rep" <>
            show (input^.replicates._1) <> ".html"
        output = printf "%s/%s_rep%d_no_doublet.mat.gz" outdir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \((idx, matFl), qcFl) -> do
        doubletScore <- detectDoublet qcFile qcPlot (qcFl^.location) (matFl^.location)
        mat <- mkSpMatrix id (matFl^.location)
        let f x = M.lookupDefault 0 x doubletScore <= thres
            n = M.size $ M.filter (<=thres) doubletScore
        runResourceT $ runConduit $
            streamRows mat .| filterC (f . fst) .|
            sinkRows n (_num_col mat) id output
        return ( (idx, location .~ output $ emptyFile)
               , location .~ qcFile $ emptyFile )
        )

detectDoublet :: FilePath
              -> FilePath
              -> FilePath
              -> FilePath
              -> IO (M.HashMap B.ByteString Double)
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
        let v = M.lookupDefault 0 (_cell_barcode s) doubletScore
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