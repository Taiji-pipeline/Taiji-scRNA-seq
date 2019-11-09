{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE OverloadedStrings #-}
module Taiji.Pipeline.SC.RNASeq.Functions.Utils
    ( groupOnC
    , getIndex
    , QC(..)
    , Annotation(..)
    , showQC
    , readQC
    , qcFileHeader
    ) where

import Data.Hashable
import qualified Data.Map.Strict                  as M
import qualified Data.ByteString.Char8                as B
import Data.Aeson
import Conduit
import Data.Binary
import           GHC.Generics (Generic)

import Taiji.Prelude

groupOnC :: (Monad m, Eq b)
         => (a -> b)
         -> ConduitT a (ConduitT i a m ()) m ()
groupOnC fun = await >>= (maybe (return ()) $ \b -> go (fun b) $ yield b)
  where
    go idx acc = await >>= maybe (yield acc) f
      where
        f b | idx' == idx = go idx (acc >> yield b)
            | otherwise = yield acc >> go idx' (yield b)
          where
            idx' = fun b
{-# INLINE groupOnC #-}

-- | Get barcode and UMI
getIndex :: B.ByteString -> (B.ByteString, B.ByteString)
getIndex x = (a,b)
  where
    [a,b] = B.split '+' $ head $ B.split '_' x
{-# INLINE getIndex #-}

data QC = QC
    { _cell_barcode :: B.ByteString
    , _num_umi :: Int
    , _uniq_gene :: Int
    , _dupRate :: Double
    , _mitoRate :: Double
    , _doubletScore :: Double
    , _genomics_context :: M.Map Annotation Int }

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
instance Binary Annotation

qcFileHeader :: B.ByteString
qcFileHeader = B.intercalate "\t" $
  [ "Num_UMI", "Num_Gene", "duplication_rate", "chrM_rate", "doublet_score"
  , "tCDS", "UTR", "Intron", "Intergenic", "Ribosomal", "Mitochondrial" ]

showQC :: QC -> B.ByteString
showQC QC{..} = B.intercalate "\t" $
    [ _cell_barcode
    , B.pack $ show _num_umi
    , B.pack $ show _uniq_gene
    , toShortest _dupRate
    , toShortest _mitoRate
    , toShortest _doubletScore
    ] ++ map (B.pack . show) dat
  where
    dat = map (\x -> M.findWithDefault 0 x _genomics_context)
        [CDS, UTR, Intron, Intergenic, Ribosomal, Mitochondrial]

readQC :: FilePath -> IO [QC]
readQC fl = do
    (_:xs) <- map (B.split '\t') . B.lines <$> B.readFile fl
    return $ map toQC xs
  where
    toQC x = QC (x!!0) (readInt $ x!!1) (readInt $ x!!2)
        (readDouble $ x!!3) (readDouble $ x!!4) (readDouble $ x!!5) $ M.fromList $
        zip [CDS, UTR, Intron, Intergenic, Ribosomal, Mitochondrial] $
        map readInt $ drop 6 x
    
