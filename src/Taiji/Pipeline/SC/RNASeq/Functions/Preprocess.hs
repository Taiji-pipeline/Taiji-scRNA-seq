{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}

module Taiji.Pipeline.SC.RNASeq.Functions.Preprocess
    ( readInput
    , getFastq
    , extractBarcode
    ) where

import Data.Conduit.Zlib (ungzip, multiple, gzip)
import           Bio.Data.Experiment.Parser
import           Bio.Data.Experiment.Types
import           Bio.Pipeline
import qualified Data.Text as T
import Shelly hiding (FilePath)
import Data.Either (rights)
import qualified Data.ByteString.Char8                as B
import Data.Conduit.Internal (zipSources)
import qualified Data.HashMap.Strict                  as M

import           Taiji.Pipeline.SC.RNASeq.Types
import Taiji.Prelude

type RAWInput = RNASeq N [Either SomeFile (SomeFile, SomeFile)]

readInput :: SCRNASeqConfig config
          => () -> ReaderT config IO [RAWInput]
readInput _ = do
    input <- asks _scrnaseq_input 
    liftIO $ simpleInputReader input "scRNA-seq" RNASeq

getFastq :: [RAWInput]
         -> [RNASeq S (SomeTags 'Fastq, SomeTags 'Fastq)]
getFastq inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map (\(x,y) -> (castFile x, castFile y)) $
        filter (\(x,y) -> getFileType x == Fastq && getFileType y == Fastq) $
        rights fls

extractBarcode :: SCRNASeqConfig config
               => RNASeq S (SomeTags 'Fastq, SomeTags 'Fastq)
               -> ReaderT config IO (RNASeq S (File '[Gzip] 'Fastq))
extractBarcode input = input & replicates.traverse.files %%~ fun
  where
    fun (flRead1, flRead2) = do
        outdir <- asks ((<> "/Preprocess") . _scrnaseq_output_dir) >>= getPath
        qdir <- qcDir
        lenCellBc <- asks _scrnaseq_cell_barcode_length
        lenUmi <- asks _scrnaseq_molecular_barcode_length
        liftIO $ do
            shelly $ test_px "umi_tools" >>= \case
                True -> return ()
                False -> error "Please install umi_tools: https://github.com/CGATOxford/UMI-tools"
            let output = printf "%s/%s_demux.fastq.gz" outdir (T.unpack $ input^.eid)
                whitelistFl = printf "%s/%s_whitelist.tsv" outdir (T.unpack $ input^.eid)
                plt = printf "%s/%s_knee_plot" qdir (T.unpack $ input^.eid)
            whitelist <- getWhiteList lenCellBc lenUmi (read1^.location) plt
            saveWhiteList whitelistFl whitelist
            runResourceT $ runConduit $ zipSources
                (streamFile $ read1^.location)
                (streamFile $ read2^.location) .| 
                    extract whitelist lenCellBc lenUmi .| unlinesAsciiC .|
                    gzip .| sinkFile output
            return $ emptyFile & location .~ output
      where
        read1 = fromSomeTags flRead1 :: File '[] 'Fastq
        read2 = fromSomeTags flRead2 :: File '[] 'Fastq
        streamFile fl = sourceFileBS fl .| multiple ungzip .|
            linesUnboundedAsciiC .| conduit
          where
            conduit = await >>= \case
                Nothing -> return ()
                Just l1 -> do
                    l2 <- fromMaybe (error "truncated fastq") <$> await
                    l3 <- fromMaybe (error "truncated fastq") <$> await
                    l4 <- fromMaybe (error "truncated fastq") <$> await
                    yield (l1, l2, l3, l4) >> conduit

type Whitelist = M.HashMap B.ByteString B.ByteString

saveWhiteList :: FilePath -> Whitelist -> IO ()
saveWhiteList output = B.writeFile output . B.unlines .
    map (\(a,b) -> a <> "\t" <> B.intercalate "," b) . M.toList .
    M.fromListWith (++) . map (\(a,b) -> (b, [a])) . M.toList
{-# INLINE saveWhiteList #-}

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

type FQ = (B.ByteString, B.ByteString, B.ByteString, B.ByteString)
    
extract :: Monad m
        => Whitelist
        -> Int    -- ^ cell barcode length
        -> Int    -- ^ umi length
        -> ConduitT (FQ, FQ) B.ByteString m ()
extract whitelist nBc nUmi = concatMapC f
  where
    f ((l1',sq,_,_), (l1,l2,_,l4)) = case M.lookup bc whitelist of
        Nothing -> []
        Just bc' ->
            let nm = "@" <> bc' <> "+" <> umi <> "_" <> B.tail l1
            in [nm,l2,"+",l4]
      where
        bc = B.take nBc sq
        umi = B.take nUmi $ B.drop nBc sq
        
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
