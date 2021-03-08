{-# OPTIONS_GHC -fno-full-laziness #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.SC.RNASeq.Functions.Internal
    ( demulti ) where

import Bio.Data.Fastq
import Data.Conduit.Internal (zipSources)
import qualified Data.ByteString.Char8                as B
import qualified Data.IntMap.Strict as I
import Conduit
import Bio.Pipeline

demulti :: FilePath -> I.IntMap Int -> Int -> Int -> FilePath -> FilePath -> IO ()
demulti out bcMap bcLen umiLen fqidxFl fqFl = runResourceT $ runConduit $
    zipSources (streamFastqGzip fqidxFl .| mapC getBc) (streamFastqGzip fqFl) .|
    concatMapC f .| sinkFastqGzip out
  where
    getBc x = let bc = dnaToInt $ B.take bcLen $ fastqSeq x
                  umi = B.take umiLen $ B.drop bcLen $ fastqSeq x
              in (bc, umi)
    f :: ((Int, B.ByteString), Fastq) -> Maybe Fastq
    f ((bc, umi), fq) = case I.lookup bc bcMap of
        Nothing -> Nothing
        Just bc' -> Just fq{fastqSeqId = B.concat [intToDna bc', "_", umi, ":", fastqSeqId fq]}

