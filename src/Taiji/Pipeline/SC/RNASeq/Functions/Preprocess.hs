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
    , getValidBarcode
    , demultiplex
    ) where

import qualified Data.IntMap.Strict as I
import           Bio.Data.Experiment.Parser
import           Bio.Data.Experiment.Types
import           Bio.Pipeline
import qualified Data.Text as T
import qualified Data.Vector.Unboxed as U
import qualified Bio.Data.Fastq as Fq
import Data.Conduit.Internal (zipSources)
import Control.Arrow ((***))
import Shelly hiding (FilePath)
import Data.Either (lefts, rights)
import qualified Data.ByteString.Char8                as B
import Language.Javascript.JMacro
import Control.DeepSeq (force)
import Data.Conduit.Async (runCConduit, (=$=&))

import           Taiji.Pipeline.SC.RNASeq.Types
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

getValidBarcode :: SCRNASeqConfig config
                => SCRNASeq S (File '[Gzip] 'Fastq, File '[Gzip] 'Fastq)
                -> ReaderT config IO ( SCRNASeq S
                    (File '[Gzip] 'Fastq, File '[Gzip] 'Fastq, File '[] 'Tsv) )
getValidBarcode input = do
    bcLen <- fromIntegral . fromMaybe (error "cell barcode length was not provided") <$>
        asks _scrnaseq_cell_barcode_length
    dir <- asks ((<> "/Fastq/Barcode/") . _scrnaseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_barcode_count.tsv" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        outputKnee = printf "%s/%s_rep%d_knee.html" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ ( \(fq1, fq2) -> liftIO $ do
        stat <- runResourceT $ runConduit $ Fq.streamFastqGzip (fq1^.location) .|
            takeC 100000000 .| barcodeStat bcLen
        B.writeFile output $ B.unlines $
            map (\(a, b) -> intToDna a <> "\t" <> (B.pack . show) b) $ U.toList stat
        thres <- fmap (read . T.unpack . head . T.lines) $ shelly $
            run "taiji-utils" ["barcode", T.pack output]
        B.writeFile output $ B.unlines $
            (("Number of valid cells" <> "\t" <> B.pack (show (thres :: Double))) :) $
            map (\(a, b) -> intToDna a <> "\t" <> (B.pack . show) b) $ U.toList stat
        kneePlot outputKnee (truncate thres) $ U.map snd stat
        return (fq1, fq2, location .~ output $ emptyFile)
        )

demultiplex :: SCRNASeqConfig config
            => SCRNASeq S (File '[Gzip] 'Fastq, File '[Gzip] 'Fastq, File '[] 'Tsv)
            -> ReaderT config IO (SCRNASeq S (File '[Demultiplexed, Gzip] 'Fastq))
demultiplex input = do
    bcLen <- fromIntegral . fromMaybe (error "cell barcode length was not provided") <$>
        asks _scrnaseq_cell_barcode_length
    umiLen <- fromIntegral . fromMaybe (error "UMI length was not provided") <$>
        asks _scrnaseq_molecular_barcode_length 
    dir <- asks ((<> "/Fastq") . _scrnaseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_demulti.fastq.gz" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ ( \(fq1, fq2, statFl) -> liftIO $ do
        content <- B.lines <$> B.readFile (statFl^.location)
        let thres = truncate $ readDouble $ last $ B.split '\t' $ head content
            stat = map ((\[a,b] -> (dnaToInt a, readInt b)) . B.split '\t') $ tail content
        let bcMap = mkBarcodeMap $ map fst $ take thres stat
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

demulti :: FilePath -> I.IntMap Int -> Int -> Int -> FilePath -> FilePath -> IO ()
demulti out bcMap bcLen umiLen fqidxFl fqFl = runResourceT $ runCConduit $
    zipSources (Fq.streamFastqGzip fqidxFl .| mapC getBc) (Fq.streamFastqGzip fqFl) =$=&
    concatMapC f =$=& Fq.sinkFastqGzip out
  where
    getBc x = let bc = dnaToInt $ B.take bcLen $ Fq.fastqSeq x
                  umi = B.take umiLen $ B.drop bcLen $ Fq.fastqSeq x
              in (bc, umi)
    f :: ((Int, B.ByteString), Fq.Fastq) -> Maybe Fq.Fastq
    f ((bc, umi), fq) = force $ case I.lookup bc bcMap of
        Nothing -> Nothing
        Just bc' -> Just fq{Fq.fastqSeqId = B.concat [intToDna bc', "_", umi, ":", Fq.fastqSeqId fq]}