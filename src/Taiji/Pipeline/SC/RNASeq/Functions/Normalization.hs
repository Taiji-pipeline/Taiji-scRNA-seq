{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.RNASeq.Functions.Normalization
    ( downSample
    ) where 

import Data.Singletons.Prelude (Elem)
import Data.Conduit.Zlib (multiple, ungzip, gzip)
import qualified Data.ByteString.Char8 as B
import Data.List.Ordered (nubSort)
import Data.Conduit.Internal (zipSinks)
import           Data.CaseInsensitive                 (original)
import Data.ByteString.Lex.Integral (packDecimal)
import Bio.Utils.Functions (scale, kdeWeight)
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import qualified Data.HashSet as S
import Data.Binary (decodeFile)
import Control.Arrow (first, second)
import qualified Data.HashMap.Strict as M
import qualified Data.Text as T
import Data.Binary (encodeFile)
import Shelly (shelly, run_)
import Data.Char (toUpper)
import Numeric.Sampling (psampleIO)
import System.Random.MWC

import           Taiji.Pipeline.SC.RNASeq.Types (SCRNASeqConfig(..))
import qualified Taiji.Utils.DataFrame as DF
import Taiji.Prelude
import Taiji.Utils

downSample :: Int -> SpMatrix Int -> IO (SpMatrix Int)
downSample n mat 
    | n >= c = return mat
    | otherwise = do
        weights <- U.map (\x -> 1 / x) . kdeWeight 1024 . U.map (\x -> exp (x / fromIntegral r) - 1) <$>
            colSum (fmap (\x -> log $ fromIntegral x + 1) mat)
        let probs = let s = U.sum weights in U.map (/s) weights
        gen <- create
        sample <- fmap (S.fromList . fromJust) $ psampleIO n $ zip (U.toList probs) [0..]
        let idx = filter (not . (`S.member` sample)) [0..c-1]
        return $ deleteCols idx mat
  where
    r = _num_row mat
    c = _num_col mat

{-
-- | Perform log normalization
logNorm :: 
    f (nm, xs) = ( nm <> "\t" <> fromJust (packDecimal totalReads)
                 , (nm, map normalize xs) )
      where
        totalReads = foldl1' (+) $ map snd xs
        normalize (i, x) = (i, log1p $ fromIntegral x / (fromIntegral totalReads / 10000))
        log1p x = log $ 1 + x :: Double
-}

