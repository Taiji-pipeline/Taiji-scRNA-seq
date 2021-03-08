{-# LANGUAGE OverloadedStrings     #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE StandaloneDeriving #-}
module Taiji.Pipeline.SC.RNASeq.Types
    ( SCRNASeq(..)
    , SCRNASeqConfig(..)
    , qcDir
    , tempDir

    , QC(..)
    , passQC
    , Annotation(..)
    , showQC
    , readQC
    , qcFileHeader
    ) where

import Bio.Data.Experiment.Types
import Bio.Data.Experiment.Replicate
import           Bio.Pipeline.Utils
import qualified Data.ByteString.Char8                as B
import Data.Hashable
import qualified Data.HashMap.Strict                  as M
import           GHC.Generics (Generic)
import Data.Binary
import Data.Aeson

import Taiji.Prelude

newtype SCRNASeq container file = SCRNASeq (CommonFields container file)
     deriving (Generic, Experiment)

deriving instance Show (container (Replicate file)) => Show (SCRNASeq container file)

instance Binary (container (Replicate file)) =>
    Binary (SCRNASeq container file)

class SCRNASeqConfig config where
    _scrnaseq_input :: config -> FilePath
    _scrnaseq_output_dir :: config -> Directory
    _scrnaseq_batch_info :: config -> Maybe FilePath
    _scrnaseq_tmp_dir :: config -> Maybe FilePath
    _scrnaseq_cell_barcode_length :: config -> Maybe Int
    _scrnaseq_molecular_barcode_length :: config -> Maybe Int
    _scrnaseq_star_index :: config -> FilePath
    _scrnaseq_genome_fasta :: config -> Maybe FilePath
    _scrnaseq_annotation :: config -> FilePath
    _scrnaseq_doublet_score_cutoff :: config -> Double
    _scrnaseq_cluster_resolution_list :: config -> [Double]
    _scrnaseq_cluster_resolution :: config -> Maybe Double
    _scrnaseq_cluster_optimizer :: config -> Optimizer

qcDir :: SCRNASeqConfig config => ReaderT config IO FilePath
qcDir = asks _scrnaseq_output_dir >>= getPath . (<> "/QC/")

tempDir :: SCRNASeqConfig config => ReaderT config IO FilePath
tempDir = asks _scrnaseq_output_dir >>= getPath . (<> "/Temporary/")

data QC = QC
    { _cell_barcode :: B.ByteString
    , _num_umi :: Int
    , _uniq_gene :: Int
    , _dupRate :: Double
    , _mitoRate :: Double
    , _doubletScore :: Double
    , _count_table :: M.HashMap Annotation Int }

passQC :: QC -> Bool
passQC QC{..} = _mitoRate <= 0.2 && _uniq_gene >= 200

-- | A region may have multiple annotations.
data Annotation = Exon
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
instance Binary Annotation

qcFileHeader :: B.ByteString
qcFileHeader = B.intercalate "\t" $
  [ "Barcode", "Num_UMI", "Num_Gene", "duplication_rate", "chrM_rate", "doublet_score"
  , "Exon", "Intron", "Intergenic", "Ribosomal", "Mitochondrial" ]

showQC :: QC -> B.ByteString
showQC QC{..} = B.intercalate "\t" $
    [ _cell_barcode
    , B.pack $ show _num_umi
    , B.pack $ show _uniq_gene
    , toShortest _dupRate
    , toShortest _mitoRate
    , toShortest _doubletScore
    ] ++ map (fromJust . packDecimal) dat
  where
    dat = map (\x -> M.lookupDefault 0 x _count_table)
        [Exon, Intron, Intergenic, Ribosomal, Mitochondrial]

readQC :: FilePath -> IO [QC]
readQC fl = do
    (_:xs) <- map (B.split '\t') . B.lines <$> B.readFile fl
    return $ map toQC xs
  where
    toQC (f1:f2:f3:f4:f5:f6:rest) = QC f1 (readInt f2) (readInt f3)
        (readDouble f4) (readDouble f5) (readDouble f6) $ M.fromList $
        zip [Exon, Intron, Intergenic, Ribosomal, Mitochondrial] $
        map readInt rest
    toQC _ = error "formatting error"