{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE TemplateHaskell   #-}

module Main where

import           Bio.Pipeline.Utils
import           Control.Lens                  ((.=))
import           Data.Aeson                    (FromJSON, ToJSON)
import Data.Binary
import           Data.Default
import           Data.Maybe                    (fromJust)
import           GHC.Generics                  (Generic)
import Data.Proxy (Proxy(..))

import           Control.Workflow
import qualified Control.Workflow.Coordinator.Drmaa as D
import Control.Workflow.Main (defaultMainWith)
import Control.Workflow.Main.Command (run, view, remote)

import qualified Taiji.Pipeline.SC.DropSeq as DropSeq
import Taiji.Pipeline.SC.DropSeq.Types (DropSeqConfig(..))

data RNASeqOpts = RNASeqOpts
    { output_dir     :: Directory
    , star_index     :: Maybe FilePath
    , rsem_index     :: Maybe FilePath
    , genome         :: Maybe FilePath
    , input          :: FilePath
    , annotation     :: Maybe FilePath
    , barcode_length :: Int
    , umi_length :: Int
    } deriving (Generic)

instance FromJSON RNASeqOpts
instance ToJSON RNASeqOpts
instance Binary RNASeqOpts

instance Default RNASeqOpts where
    def = RNASeqOpts
        { output_dir = asDir "output"
        , star_index = Nothing
        , rsem_index = Nothing
        , genome = Nothing
        , input = "input.yml"
        , annotation = Nothing
        , barcode_length = 12
        , umi_length = 8
        }

instance DropSeqConfig RNASeqOpts where
    _dropseq_input = input
    _dropseq_output_dir = output_dir
    _dropseq_cell_barcode_length = barcode_length
    _dropseq_molecular_barcode_length = umi_length
    _dropseq_star_index = fromJust . star_index
    _dropseq_annotation = fromJust . annotation
    _dropseq_genome_fasta = genome

build "wf" [t| SciFlow RNASeqOpts |] DropSeq.builder

decodeDrmaa :: String -> Int -> FilePath -> IO D.DrmaaConfig
decodeDrmaa ip port _ = do
    config <- D.getDefaultDrmaaConfig
    let exePath = fst $ D._cmd config
    return $ config { D._cmd = (exePath,
        ["remote", "--ip", ip, "--port", show port]) }

main :: IO ()
main = defaultMainWith cmd wf
  where
    cmd = [ ("run", run decodeDrmaa)
          , ("view", view)
          , ("remote", remote (Proxy :: Proxy D.Drmaa))
          ]

