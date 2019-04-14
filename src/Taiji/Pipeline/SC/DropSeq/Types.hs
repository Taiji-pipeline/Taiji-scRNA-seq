{-# LANGUAGE DeriveGeneric #-}
module Taiji.Pipeline.SC.DropSeq.Types
    ( DropSeqConfig(..)
    , Annotation(..)
    ) where

import Data.Hashable
import           Bio.Pipeline.Utils
import Data.Aeson
import Data.Serialize
import           GHC.Generics (Generic)

class DropSeqConfig config where
    _dropseq_input :: config -> FilePath
    _dropseq_output_dir :: config -> Directory
    _dropseq_cell_barcode_length :: config -> Int
    _dropseq_molecular_barcode_length :: config -> Int
    _dropseq_star_index :: config -> FilePath
    _dropseq_genome_fasta :: config -> Maybe FilePath
    _dropseq_annotation :: config -> FilePath

-- | A region may have multiple annotations.
data Annotation = CDS
                | UTR
                | Exon
                | Intron
                | Ribosomal
                | Mitochondrial
                | Genic
                | Intergenic
                deriving (Show, Eq, Ord, Generic)

instance Hashable Annotation
instance FromJSON Annotation
instance ToJSON Annotation
instance FromJSONKey Annotation
instance ToJSONKey Annotation
instance Serialize Annotation

