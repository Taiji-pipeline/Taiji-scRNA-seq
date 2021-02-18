{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}

module Taiji.Pipeline.SC.RNASeq.Functions.Align
    ( mkIndex
    , tagAlign
    , filterNameSortBam
    ) where

import           Bio.Pipeline
import           Data.Either                          (fromLeft)
import System.IO.Temp (withTempFile)
import qualified Data.Text                            as T
import Shelly hiding (FilePath)

import           Taiji.Pipeline.SC.RNASeq.Types
import Taiji.Prelude

mkIndex :: SCRNASeqConfig config => [a] -> ReaderT config IO [a]
mkIndex input
    | null input = return input
    | otherwise = do
        genome <- fromMaybe (error "genome fasta not found") <$>
            asks _scrnaseq_genome_fasta
        starIndex <- asks _scrnaseq_star_index
        anno <- asks _scrnaseq_annotation
        liftIO $ do
            _ <- starMkIndex "STAR" starIndex [genome] anno 100
            return input

tagAlign :: SCRNASeqConfig config
         => SCRNASeq S (File '[Gzip] 'Fastq)
         -> ReaderT config IO (SCRNASeq S (File '[] 'Bam))
tagAlign input = do
    dir <- asks ((<> "/Bam") . _scrnaseq_output_dir) >>= getPath
    idx <- asks _scrnaseq_star_index
    let outputGenome = printf "%s/%s_rep%d_genome.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        f fl = starAlign outputGenome idx (Left fl) opt >>=
            return . fst . fromLeft undefined
    input & replicates.traverse.files %%~ liftIO . f
  where
    opt = defaultSTAROpts & starCores .~ 8 & starTranscriptome .~ Nothing

-- | Filter bad quality reads and name sort Bam file.
filterNameSortBam :: SCRNASeqConfig config
                  => SCRNASeq S (File '[] 'Bam)
                  -> ReaderT config IO (SCRNASeq S (File '[NameSorted] 'Bam))
filterNameSortBam input = do
    dir <- asks ((<> "/Bam") . _scrnaseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_filt.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \fl -> do
        fun fl output
        return $ location .~ output $ emptyFile )
  where
    fun fl output = withTempFile "./" "tmp_sort." $ \tmp_sort _ ->
        shelly $ escaping False $ silently $ bashPipeFail bash_ "samtools"
            [ "view", "-F", "0x70c", "-u", T.pack $ fl^.location, "|"
            ,  "samtools", "sort", "-", "-n", "-T", T.pack tmp_sort
            , "-l", "9", "-o", T.pack output ]