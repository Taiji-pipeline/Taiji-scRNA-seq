module Taiji.Pipeline.SC.DropSeq.Config
    ( DropSeqConfig(..)
    ) where

import           Bio.Pipeline.Utils

class DropSeqConfig config where
    _dropSeq_input :: config -> FilePath
    _dropSeq_output_dir :: config -> Directory
    _dropSeq_cell_barcode_length :: config -> Int
    _dropSeq_molecular_barcode_length :: config -> Int
    _dropSeq_star_index :: config -> FilePath
    _dropSeq_genome_fasta :: config -> FilePath
    _dropSeq_annotation :: config -> FilePath
