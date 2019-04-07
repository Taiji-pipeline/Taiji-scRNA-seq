{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE LambdaCase #-}

module Taiji.Pipeline.SC.DropSeq.Functions
    ( extractBarcode
    , mkIndex
    , tagAlign
    , filterBamFile
    , quantification
    ) where

import           Bio.Data.Bam
import           Bio.Data.Bed
import           Bio.Data.Experiment
import           Bio.Pipeline.NGS.STAR
import           Bio.Pipeline.NGS.Utils
import           Bio.Pipeline.Utils
import           Bio.RealWorld.GENCODE
import           Conduit
import           Control.Lens
import           Control.Monad
import           Control.Monad.IO.Class               (liftIO)
import           Control.Monad.Reader                 (asks)
import qualified Data.ByteString.Char8                as B
import           Data.CaseInsensitive                 (original)
import           Data.Either                          (fromLeft)
import qualified Data.HashMap.Strict                  as M
import qualified Data.IntervalMap.Strict              as IM
import qualified Data.IntSet                          as S
import           Data.List
import           Data.Maybe
import qualified Data.Text                            as T
import System.IO
import           Scientific.Workflow
import Shelly hiding (FilePath)
import           Text.Printf                          (printf)

import           Taiji.Pipeline.SC.DropSeq.Config

extractBarcode :: DropSeqConfig config
               => RNASeq S (SomeTags 'Fastq, SomeTags 'Fastq)
               -> WorkflowConfig config (RNASeq S (File '[Gzip] 'Fastq))
extractBarcode input = input & replicates.traverse.files %%~ fun
  where
    fun (flRead1, flRead2) = do
        outdir <- asks _dropSeq_output_dir >>= getPath
        liftIO $ do
            shelly $ test_px "umi_tools" >>= \case
                True -> return ()
                False -> error "Please install umi_tools: https://github.com/CGATOxford/UMI-tools"
            let whitelist = printf "%s/%s_whitelist.txt" outdir (T.unpack $ input^.eid)
                output = printf "%s/%s_extract.fastq.gz" outdir (T.unpack $ input^.eid)
            getWhiteList (read1^.location) whitelist
            extract (read1^.location) (read2^.location) whitelist output
            return $ emptyFile & location .~ output
      where
        read1 = fromSomeTags flRead1 :: File '[] 'Fastq
        read2 = fromSomeTags flRead2 :: File '[] 'Fastq

getWhiteList :: FilePath -> FilePath -> IO ()
getWhiteList input output = shelly $ escaping False $ run_ "umi_tools"
    [ "whitelist", "--stdin", T.pack input
    , "--plot-prefix=" <> T.init (fst $ T.breakOnEnd "." $ T.pack output)
    , "--bc-pattern=CCCCCCCCCCCCNNNNNNNN"
    , "--log2stderr", ">", T.pack output ]

extract :: FilePath -> FilePath -> FilePath -> FilePath -> IO ()
extract input1 input2 whitelist output = shelly $ run_ "umi_tools"
    [ "extract"
    , "--quality-filter-threshold=10"
    , "--quality-encoding=phred33"
    , "--bc-pattern=CCCCCCCCCCCCNNNNNNNN"
    , "--stdin", T.pack input1
    , "--stdout=/dev/null"
    , "--read2-in", T.pack input2
    , "--read2-out=" <> T.pack output
    , "--filter-cell-barcode"
    , "--error-correct-cell"
    , "--whitelist=" <> T.pack whitelist ]

mkIndex :: DropSeqConfig config => [a] -> WorkflowConfig config [a]
mkIndex input
    | null input = return input
    | otherwise = do
        genome <- asks _dropSeq_genome_fasta
        starIndex <- asks _dropSeq_star_index
        anno <- asks _dropSeq_annotation
        liftIO $ do
            _ <- starMkIndex "STAR" starIndex [genome] anno 100
            return input

tagAlign :: DropSeqConfig config
         => RNASeq S (File '[Gzip] 'Fastq)
         -> WorkflowConfig config (RNASeq S (File '[] 'Bam))
tagAlign input = do
    dir <- asks _dropSeq_output_dir >>= getPath
    idx <- asks _dropSeq_star_index
    let outputGenome = printf "%s/%s_rep%d_genome.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        f fl = starAlign outputGenome idx (Left fl) opt >>=
            return . fst . fromLeft undefined
    input & replicates.traverse.files %%~ liftIO . f
  where
    opt = defaultSTAROpts & starCores .~ 4 & starTranscriptome .~ Nothing

-- | Filter bad quality reads and name sort Bam file.
filterBamFile :: DropSeqConfig config
              => RNASeq S (File '[] 'Bam)
              -> WorkflowConfig config (RNASeq S (File '[] 'Bam))
filterBamFile input = do
    dir <- asks _dropSeq_output_dir >>= getPath
    let output = printf "%s/%s_rep%d_filt.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \fl ->
            filterBam "./" output fl )

quantification :: DropSeqConfig config
               => RNASeq S (File tags 'Bam)
               -> WorkflowConfig config (RNASeq S (File '[] 'Tsv))
quantification input = do
    dir <- asks ((<> "/Quantification/") . _dropSeq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_quant.tsv" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    anno_f <- asks _dropSeq_annotation
    genes <- liftIO $ readGenes anno_f
    let annotation = bedToTree (++) $ map (\Gene{..} ->
            (asBed geneChrom geneLeft geneRight :: BED3, [original geneName])) genes
    input & replicates.traverse.files %%~ liftIO . ( \bam -> do
        hdr <- getBamHeader $ bam^.location
        runResourceT $ runConduit $ streamBam (bam^.location) .|
            bamToBedC hdr .| quantify output annotation
        )

quantify :: MonadIO m
         => FilePath
         -> BEDTree [B.ByteString]  -- ^ genes associated with regions
         -> ConduitT BED o m (File '[] 'Tsv)
quantify output annotation = do
    res <- mkCountTable annotation
    liftIO $ outputTable output res
    return $ emptyFile & location .~ output
{-# INLINE quantify #-}

type GeneCount = M.HashMap B.ByteString S.IntSet
type CountTable = M.HashMap Int GeneCount

outputTable :: FilePath -> CountTable -> IO ()
outputTable output t = withFile output WriteMode $ \h -> do
    B.hPutStrLn h $ B.intercalate "\t" genes
    forM_ (M.toList t) $ \(_, geneCount) -> do
        let geneCount' = fmap S.size geneCount
            res = flip map genes $ \g -> B.pack $ show $ M.lookupDefault 0 g geneCount'
        B.hPutStrLn h $ B.intercalate "\t" res
  where
    genes = nub $ concatMap M.keys $ M.elems t
  

mkCountTable :: Monad m
             => BEDTree [B.ByteString]  -- ^ genes associated with regions
             -> ConduitT BED o m CountTable
mkCountTable geneMap = foldlC f M.empty
  where
    f m x = addToMap geneMap x m

addToMap :: BEDTree [B.ByteString]  -- ^ genes associated with regions
         -> BED
         -> CountTable
         -> CountTable
addToMap geneMap bed t = M.alter fun barcode t
  where
    fun Nothing = Just new
    fun (Just old) = Just $ M.unionWith S.union old new
    genes = concat $ IM.elems $ intersecting geneMap bed
    new = M.fromList $ zip genes $ repeat $ S.singleton umi
    (barcode, umi) = getIndex bed
{-# INLINE addToMap #-}

getIndex :: BED -> (Int, Int)
getIndex bed = (dnaToInt b, dnaToInt a)
  where
    (a:b:_) = reverse $ B.split '_' $ head $ B.words $ fromJust $ bed^.name
{-# INLINE getIndex #-}

-- | Convert quinary DNA sequence (A: 0, C: 1, G: 2, T: 3, N: 4) to decimal.
dnaToInt :: B.ByteString -> Int
dnaToInt = fst . B.foldl' f (0, 1)
  where
    f (acc, i) 'A' = (acc, i * 5)
    f (acc, i) 'C' = (acc + 1 * i, i * 5)
    f (acc, i) 'G' = (acc + 2 * i, i * 5)
    f (acc, i) 'T' = (acc + 3 * i, i * 5)
    f (acc, i) 'N' = (acc + 4 * i, i * 5)
    f _ _          = error "Unexpected character!"
{-# INLINE dnaToInt #-}

-- | Convert decimal to quinary DNA sequence.
intToDna :: Int   -- ^ length of the resulting bytestring
         -> Int
         -> B.ByteString
intToDna n = B.pack . reverse . go 1 []
  where
    go !i !acc !x
        | m == 0 && i >= n = c : acc
        | otherwise = go (i+1) (c:acc) m
      where
        c = case x `mod` 5 of
            0 -> 'A'
            1 -> 'C'
            2 -> 'G'
            3 -> 'T'
            4 -> 'N'
            _ -> undefined
        m = x `div` 5
{-# INLINE intToDna #-}
