{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}

module Taiji.Pipeline.SC.DropSeq.Functions.Preprocess
    ( readInput
    , getFastq
    , extractBarcode
    ) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Types
import           Bio.Data.Experiment.Parser
import           Bio.Pipeline
import           Control.Monad.IO.Class               (liftIO)
import qualified Data.Text as T
import Shelly hiding (FilePath)
import Control.Lens
import           Control.Monad.Reader                 (asks)
import           Text.Printf                          (printf)
import Data.Either (rights)
import qualified Data.ByteString.Char8                as B
import Conduit
import Data.Conduit.Internal (zipSources)
import Bio.Pipeline.Barcode (deBarcode, unBarCode)
import qualified Bio.Data.Fastq as F
import           Scientific.Workflow
import qualified Data.HashMap.Strict                  as M

import           Taiji.Pipeline.SC.DropSeq.Types

type RAWInput = RNASeq N [Either SomeFile (SomeFile, SomeFile)]

readInput :: DropSeqConfig config
          => () -> WorkflowConfig config [RAWInput]
readInput _ = do
    input <- asks _dropseq_input 
    liftIO $ simpleInputReader input "drop-seq" RNASeq

getFastq :: [RAWInput]
         -> [RNASeq S (SomeTags 'Fastq, SomeTags 'Fastq)]
getFastq inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map (\(x,y) -> (castFile x, castFile y)) $
        filter (\(x,y) -> getFileType x == Fastq && getFileType y == Fastq) $
        rights fls

extractBarcode :: DropSeqConfig config
               => RNASeq S (SomeTags 'Fastq, SomeTags 'Fastq)
               -> WorkflowConfig config (RNASeq S (File '[Gzip] 'Fastq))
extractBarcode input = input & replicates.traverse.files %%~ fun
  where
    fun (flRead1, flRead2) = do
        outdir <- asks ((<> "/Preprocess") . _dropseq_output_dir) >>= getPath
        lenCellBc <- asks _dropseq_cell_barcode_length
        lenUmi <- asks _dropseq_molecular_barcode_length
        liftIO $ do
            shelly $ test_px "umi_tools" >>= \case
                True -> return ()
                False -> error "Please install umi_tools: https://github.com/CGATOxford/UMI-tools"
            let output = printf "%s/%s_extract.fastq.gz" outdir (T.unpack $ input^.eid)
                plt = printf "%s/%s_whitelist" outdir (T.unpack $ input^.eid)
            whitelist <- getWhiteList lenCellBc lenUmi (read1^.location) plt
            runResourceT $ runConduit $ zipSources
                (F.streamFastqGzip $ read1^.location)
                (F.streamFastqGzip $ read2^.location) .| 
                    extract whitelist lenCellBc lenUmi .| F.sinkFastqGzip output
            return $ emptyFile & location .~ output
      where
        read1 = fromSomeTags flRead1 :: File '[] 'Fastq
        read2 = fromSomeTags flRead2 :: File '[] 'Fastq

type Whitelist = M.HashMap B.ByteString B.ByteString

getWhiteList :: Int   -- ^ cell barcode length
             -> Int   -- ^ umi length
             -> FilePath
             -> FilePath
             -> IO Whitelist
getWhiteList nBc nUmi input plt = fmap
    (M.fromList . concatMap (f . B.split '\t') . B.lines . B.pack . T.unpack) $
    shelly $ run "umi_tools"
        [ "whitelist", "--stdin", T.pack input
        , "--plot-prefix=" <> T.pack plt
        , "--bc-pattern=" <> T.pack (replicate nBc 'C' <> replicate nUmi 'N')
        , "--log2stderr" ]
  where
    f (a:b:_) = zip (a : B.split ',' b) $ repeat a

extract :: Monad m
        => Whitelist
        -> Int    -- ^ cell barcode length
        -> Int    -- ^ umi length
        -> ConduitT (F.Fastq, F.Fastq) F.Fastq m ()
extract whitelist nBc nUmi = concatMapC f
  where
    f (read1, read2) = case deBarcode [nBc, nUmi] read1 of
        Just ([bc, umi], _) -> case M.lookup (unBarCode bc) whitelist of
            Nothing -> Nothing
            Just bc' ->
                let nm = B.intercalate "+" [bc', unBarCode umi] <> "_" <> F.fastqSeqId read2
                in Just $ read2{F.fastqSeqId = nm}
        _ -> Nothing

{-
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
    -}

