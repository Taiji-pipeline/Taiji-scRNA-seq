{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.RNASeq.Functions.Clustering
    ( filterMatrix
    , spectral
    ) where 

import Data.Singletons.Prelude (Elem)
import qualified Data.ByteString.Char8 as B
import Data.ByteString.Lex.Integral (packDecimal)
import Bio.Utils.Functions (scale)
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import qualified Data.Text as T
import Shelly (shelly, run_)

import           Taiji.Pipeline.SC.RNASeq.Types
import Taiji.Prelude
import Taiji.Utils

filterMatrix :: (Elem 'Gzip tags ~ 'True, SCRNASeqConfig config)
             => FilePath
             -> RNASeq S (File tags 'Other)
             -> ReaderT config IO (RNASeq S (File '[] 'Tsv, File tags 'Other))
filterMatrix prefix input = do
    dir <- asks ((<> asDir prefix) . _scrnaseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_filt.mat.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
        rownames = printf "%s/%s_rep%d_rownames.txt" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ ( \fl -> liftIO $ do
        sp <- mkSpMatrix readInt $ fl^.location
        runResourceT $ runConduit $
            streamRows sp .| mapC f .| unlinesAsciiC .| sinkFile rownames
        counts <- do
            v <- UM.replicate (_num_col sp) 0
            runResourceT $ runConduit $ streamRows sp .| concatMapC snd .|
                mapC fst .| mapM_C (UM.unsafeModify v (+1))
            U.unsafeFreeze v
        let (zeros, nonzeros) = U.partition ((==0) . snd) $
                U.zip (U.enumFromN 0 (U.length counts)) counts
            (i, v) = U.unzip nonzeros
            idx = U.toList $ fst $ U.unzip $ U.filter ((>1.65) . snd) $ U.zip i $ scale v
        filterCols output (idx ++ U.toList (fst $ U.unzip zeros)) $ fl^.location
        return ( location .~ rownames $ emptyFile
               , location .~ output $ emptyFile ) )
  where
    f (nm, xs) = nm <> "\t" <> fromJust (packDecimal $ foldl1' (+) $ map snd xs)

-- | Reduce dimensionality using spectral clustering
spectral :: (Elem 'Gzip tags ~ 'True, SCRNASeqConfig config)
         => FilePath  -- ^ directory
         -> Maybe Int  -- ^ seed
         -> RNASeq S (a, File tags 'Other)
         -> ReaderT config IO (RNASeq S (a, File '[Gzip] 'Tsv))
spectral prefix seed input = do
    dir <- asks ((<> asDir prefix) . _scrnaseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_spectral.tsv.gz" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \(rownames, fl) -> do
        shelly $ run_ "taiji-utils" $ ["reduce", "--method", "spectral"
            , "--distance", "cosine"
            , T.pack $ fl^.location, T.pack output] ++ maybe []
            (\x -> ["--seed", T.pack $ show x]) seed
        return (rownames, location .~ output $ emptyFile)
        )

{-
clustering :: SCATACSeqConfig config
           => FilePath
           -> ClustOpt
           -> SCATACSeq S [(File '[] 'Tsv, File '[Gzip] 'Tsv)]
           -> ReaderT config IO (SCATACSeq S (File '[] 'Other))
clustering prefix opt input = do
    tmp <- asks _scatacseq_tmp_dir
    dir <- asks ((<> asDir ("/" ++ prefix)) . _scatacseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_clusters.bin" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \fl -> do
        clustering' opt tmp fl >>= encodeFile output
        return $ location .~ output $ emptyFile )
        -}

{-
clustering' :: Maybe FilePath      -- ^ temp dir
            -> [(File '[] 'Tsv, File '[Gzip] 'Tsv)]   -- ^ lsa input matrix
            -> IO [CellCluster]
clustering' dir fls = withTempDir dir $ \tmpD -> do
      let sourceCells = getZipSource $ (,,) <$>
              ZipSource (iterateC succ 0) <*>
              ZipSource seqDepthC <*>
              ZipSource ( sourceFile (tmpD <> "/embed") .| linesUnboundedAsciiC .|
                mapC (map readDouble . B.split '\t') )
          input = T.pack $ intercalate "," $ map (^.location) mats
      shelly $ run_ "taiji-utils" $
          [ "clust", input, T.pack tmpD <> "/clust"
          , "--embed", T.pack tmpD <> "/embed"
          , "--embed-method", "umap"
          , "--dim" 50 ]
      cells <- runResourceT $ runConduit $ sourceCells .| mapC f .| sinkVector
      clusters <- readClusters $ tmpD <> "/clust"
      return $ zipWith (\i -> CellCluster $ B.pack $ "C" ++ show i) [1::Int ..] $
          map (map (cells V.!)) clusters
  where
    coverage = head $ fst $ unzip fls
    mats = snd $ unzip fls
    readClusters fl = map (map readInt . B.split ',') . B.lines <$>
        B.readFile fl
    seqDepthC = sourceFile (coverage^.location) .| linesUnboundedAsciiC .|
        mapC ((\[a,b] -> (a,b)) . B.split '\t')
    f (i, (bc, dep), [d1,d2,d3,d4,d5]) = Cell i (d1,d2) (d3,d4,d5) bc $ readInt dep
    f _ = error "formatting error"
{-# INLINE clustering' #-}
-}