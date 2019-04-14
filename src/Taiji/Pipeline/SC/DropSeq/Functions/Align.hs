{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}

module Taiji.Pipeline.SC.DropSeq.Functions.Align
    ( mkIndex
    , tagAlign
    , filterBamSort
    ) where

import           Bio.Data.Experiment
import           Bio.Pipeline
import           Control.Lens
import           Control.Monad.IO.Class               (liftIO)
import           Control.Monad.Reader                 (asks)
import           Data.Either                          (fromLeft)
import System.IO.Temp (withTempFile)
import           Data.Maybe
import qualified Data.Text                            as T
import           Scientific.Workflow
import Shelly hiding (FilePath)
import           Text.Printf                          (printf)

import           Taiji.Pipeline.SC.DropSeq.Types

mkIndex :: DropSeqConfig config => [a] -> WorkflowConfig config [a]
mkIndex input
    | null input = return input
    | otherwise = do
        genome <- fromMaybe (error "genome fasta not found") <$>
            asks _dropseq_genome_fasta
        starIndex <- asks _dropseq_star_index
        anno <- asks _dropseq_annotation
        liftIO $ do
            _ <- starMkIndex "STAR" starIndex [genome] anno 100
            return input

tagAlign :: DropSeqConfig config
         => RNASeq S (File '[Gzip] 'Fastq)
         -> WorkflowConfig config (RNASeq S (File '[] 'Bam))
tagAlign input = do
    dir <- asks ((<> "/Bam") . _dropseq_output_dir) >>= getPath
    idx <- asks _dropseq_star_index
    let outputGenome = printf "%s/%s_rep%d_genome.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        f fl = starAlign outputGenome idx (Left fl) opt >>=
            return . fst . fromLeft undefined
    input & replicates.traverse.files %%~ liftIO . f
  where
    opt = defaultSTAROpts & starCores .~ 4 & starTranscriptome .~ Nothing

-- | Filter bad quality reads and name sort Bam file.
filterBamSort :: DropSeqConfig config
              => RNASeq S (File '[] 'Bam)
              -> WorkflowConfig config (RNASeq S (File '[NameSorted] 'Bam))
filterBamSort input = do
    dir <- asks ((<> "/Bam") . _dropseq_output_dir) >>= getPath
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