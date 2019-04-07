module Taiji.Pipeline.SC.DropSeq.Types
    ( DropSeqConfig(..)
    ) where

import           Bio.Pipeline.Utils

class DropSeqConfig config where
    _dropseq_input :: config -> FilePath
    _dropseq_output_dir :: config -> Directory
    _dropseq_cell_barcode_length :: config -> Int
    _dropseq_molecular_barcode_length :: config -> Int
    _dropseq_star_index :: config -> FilePath
    _dropseq_genome_fasta :: config -> Maybe FilePath
    _dropseq_annotation :: config -> FilePath
