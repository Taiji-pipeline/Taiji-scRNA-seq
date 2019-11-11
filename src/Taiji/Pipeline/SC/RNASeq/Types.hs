{-# LANGUAGE OverloadedStrings     #-}
module Taiji.Pipeline.SC.RNASeq.Types
    ( SCRNASeqConfig(..)
    , qcDir
    ) where

import           Bio.Pipeline.Utils

import Taiji.Prelude

class SCRNASeqConfig config where
    _scrnaseq_input :: config -> FilePath
    _scrnaseq_output_dir :: config -> Directory
    _scrnaseq_tmp_dir :: config -> Maybe FilePath
    _scrnaseq_cell_barcode_length :: config -> Int
    _scrnaseq_molecular_barcode_length :: config -> Int
    _scrnaseq_star_index :: config -> FilePath
    _scrnaseq_genome_fasta :: config -> Maybe FilePath
    _scrnaseq_annotation :: config -> FilePath

qcDir :: SCRNASeqConfig config => ReaderT config IO FilePath
qcDir = asks _scrnaseq_output_dir >>= getPath . (<> "/QC/")