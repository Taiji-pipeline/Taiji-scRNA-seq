{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE LambdaCase #-}

module Taiji.Pipeline.SC.RNASeq.Functions.Preprocess
    ( readInput
    , getFastq
    , getDemultiplexedFastq
    , getMatrix
    , demultiplex
    ) where

import           Bio.Data.Experiment.Parser
import           Bio.Data.Experiment.Types
import           Bio.Pipeline
import qualified Data.Text as T
import qualified Data.Vector.Unboxed as U
import Bio.Data.Fastq (streamFastqGzip)
import Control.Arrow ((***))
import Shelly hiding (FilePath)
import Data.Either (lefts, rights)
import qualified Data.ByteString.Char8                as B
import Language.Javascript.JMacro

import           Taiji.Pipeline.SC.RNASeq.Types
import Taiji.Pipeline.SC.RNASeq.Functions.Internal
import Taiji.Utils.Plot
import Taiji.Utils.Plot.ECharts
import Taiji.Prelude

type RAWInput = SCRNASeq N [Either SomeFile (SomeFile, SomeFile)]

readInput :: SCRNASeqConfig config
          => () -> ReaderT config IO [RAWInput]
readInput _ = do
    input <- asks _scrnaseq_input 
    liftIO $ mkInputReader input "scRNA-seq" SCRNASeq

getFastq :: [RAWInput]
         -> [SCRNASeq S (File '[Gzip] 'Fastq, File '[Gzip] 'Fastq)]
getFastq inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~
        map (fromSomeFile *** fromSomeFile) . filter f . rights
  where
    f (x,y) = getFileType x == Fastq && getFileType y == Fastq &&
        x `hasTag` Gzip && y `hasTag` Gzip

getDemultiplexedFastq :: [RAWInput] -> [SCRNASeq S (File '[Demultiplexed, Gzip] 'Fastq)]
getDemultiplexedFastq inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ map fromSomeFile . filter f . lefts
  where
    f x = getFileType x == Fastq && x `hasTag` Demultiplexed && x `hasTag` Gzip

getMatrix :: [RAWInput]
          -> [SCRNASeq S ( File '[RowName, Gzip] 'Tsv
                         , File '[ColumnName, Gzip] 'Tsv
                         , File '[Gzip] 'MatrixMarket ) ]
getMatrix inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f . lefts
  where
    f fls = case (getRow fls, getCol fls, getMat fls) of
        (Just row, Just col, Just mat) -> [(row, col, mat)]
        (Nothing, Nothing, Nothing) -> []
        _ -> error "Incomplete matrix input"
    getMat fls = case filter (\x -> x `hasTag` Gzip && getFileType x == MatrixMarket) fls of
        [] -> Nothing
        [x] -> Just $ fromSomeFile x
        _ -> error "Found multiple matrix files in the input"
    getRow fls = case filter (\x -> x `hasTag` Gzip && x `hasTag` RowName) fls of
        [] -> Nothing
        [x] -> Just $ fromSomeFile x
        _ -> error "Found multiple row name files in the input"
    getCol fls = case filter (\x -> x `hasTag` Gzip && x `hasTag` ColumnName) fls of
        [] -> Nothing
        [x] -> Just $ fromSomeFile x
        _ -> error "Found multiple column name files in the input"

demultiplex :: SCRNASeqConfig config
            => SCRNASeq S (File '[Gzip] 'Fastq, File '[Gzip] 'Fastq)
            -> ReaderT config IO (SCRNASeq S (File '[Demultiplexed, Gzip] 'Fastq))
demultiplex input = do
    bcLen <- fromMaybe (error "cell barcode length was not provided") <$> asks _scrnaseq_cell_barcode_length
    umiLen <- fromMaybe (error "UMI length was not provided") <$> asks _scrnaseq_molecular_barcode_length 
    dir <- asks ((<> "/Fastq") . _scrnaseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_demulti.fastq.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        outputKnee = printf "%s/%s_rep%d_knee.html" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ ( \(fq1, fq2) -> liftIO $ withTemp Nothing $ \tmp -> do
        stat <- runResourceT $ runConduit $ streamFastqGzip (fq1^.location) .|
            takeC 100000000 .| barcodeStat bcLen
        B.writeFile tmp $ B.unlines $ map (B.pack . show . snd) $ U.toList stat
        thres <- fmap (read . T.unpack . head . T.lines) $ shelly $
            run "taiji-utils" ["barcode", T.pack tmp]
        kneePlot outputKnee (truncate thres) $ U.map snd stat
        let bcMap = mkBarcodeMap $ map fst $ take (truncate (thres :: Double)) $ U.toList stat
        demulti output bcMap bcLen umiLen (fq1^.location) (fq2^.location)
        return $ location .~ output $ emptyFile
        )

kneePlot :: FilePath
         -> Int        -- ^ threshold
         -> U.Vector Int -- ^ Sorted list
         -> IO ()
kneePlot output thres input = savePlots output [] [plt']
  where
    plt' = addAttr [jmacroE| {
        yAxis: {
            type: "log"
        },
        xAxis: {
            type: "log"
        },
        visualMap: {
            show: true,
            dimension: 0,
            pieces:
                [ {min: 1, max:`thres`, color:"red", label:"Cell"}
                , {min: `thres`, color: "gray", label:"Background"} ]
        }
        }|] plt
    plt = line' $ map (fromIntegral *** fromIntegral) $ reverse $
        fst $ U.ifoldl f ([(1, U.head input)], (0, U.head input)) input
    f (acc, (i', x')) i x
        | x' == x = (acc, (i', x'))
        | x' /= x && i == i' + 1 = ((i+1, x) : acc, (i, x))
        | otherwise = ((i+1, x) : (i, x') : acc, (i, x))

{-
-- | Extract barcode: @barcode+umi:seq_name
extractBarcode :: SCRNASeqConfig config
               => SCRNASeq S (SomeTags 'Fastq, SomeTags 'Fastq)
               -> ReaderT config IO (SCRNASeq S (File '[Gzip] 'Fastq))
extractBarcode input = input & replicates.traverse.files %%~ fun
  where
    fun (flRead1, flRead2) = do
        outdir <- asks ((<> "/Fastq") . _scrnaseq_output_dir) >>= getPath
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
            demultiplex output (read1^.location) (read2^.location)
                whitelist lenCellBc lenUmi
            return $ emptyFile & location .~ output
      where
        read1 = fromSomeTags flRead1 :: File '[] 'Fastq
        read2 = fromSomeTags flRead2 :: File '[] 'Fastq

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
    -}