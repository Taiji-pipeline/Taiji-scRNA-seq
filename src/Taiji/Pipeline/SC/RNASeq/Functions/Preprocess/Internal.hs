{-# OPTIONS_GHC -fno-full-laziness #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.SC.RNASeq.Functions.Preprocess.Internal
    (demultiplex) where

import Data.Conduit.Zlib (ungzip, multiple, gzip)
import Data.Conduit.Internal (zipSources)
import qualified Data.ByteString.Char8                as B
import qualified Data.HashMap.Strict                  as M
import Data.Maybe
import Conduit

type Whitelist = M.HashMap B.ByteString B.ByteString

demultiplex :: FilePath
            -> FilePath
            -> FilePath
            -> Whitelist
            -> Int
            -> Int
            -> IO ()
demultiplex output read1 read2 whitelist lenCellBc lenUmi = runResourceT $
    runConduit $ source .| extract whitelist lenCellBc lenUmi .|
    unlinesAsciiC .| gzip .| sinkFile output
  where
    source = zipSources (streamFile read1) (streamFile read2)
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

type FQ = (B.ByteString, B.ByteString, B.ByteString, B.ByteString)
    
extract :: Monad m
        => Whitelist
        -> Int    -- ^ cell barcode length
        -> Int    -- ^ umi length
        -> ConduitT (FQ, FQ) B.ByteString m ()
extract whitelist nBc nUmi = concatMapC f
  where
    f ((_,sq,_,_), (l1,l2,_,l4)) = case M.lookup bc whitelist of
        Nothing -> []
        Just bc' ->
            let nm = "@" <> bc' <> "+" <> umi <> "_" <> B.tail l1
            in [nm,l2,"+",l4]
      where
        bc = B.take nBc sq
        umi = B.take nUmi $ B.drop nBc sq