{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.Normalise #-}
module Taiji.Pipeline.SC.RNASeq.Functions.Spectral
    ( Spectral(..)
    , spectral
    , spectralEmbedC
    , getSpectral
    , rbf
    ) where

import Data.Matrix.Static.LinearAlgebra
import qualified Data.ByteString.Char8 as B
import Control.Monad.ST (runST)
import qualified Data.Vector.Storable as V
import qualified Data.Matrix.Static.Generic as G
import qualified Data.Matrix.Static.Generic.Mutable as GM
import qualified Data.Matrix.Static.Dense as D
import qualified Data.Matrix.Static.Sparse as S
import qualified Data.Matrix.Storable as MS
import Data.Matrix.Dynamic (Dynamic(..), withDyn, fromVector)
import qualified Data.Vector.Unboxed as U
import Data.Singletons.Prelude hiding ((@@), type (==))
import Data.Complex
import Data.Conduit.List (chunksOf)
import Conduit
import Data.Type.Equality
import Data.Singletons.TypeLits
import Data.Singletons.Decide (decideEquality)
import Numeric.Sampling (sample)
import System.Random.MWC
import Flat (flat, unflat)

import Taiji.Prelude
import           Taiji.Pipeline.SC.RNASeq.Types (SCRNASeqConfig(..), SCRNASeq)

getSpectral :: SCRNASeqConfig config
            => Int     -- ^ Number of samples
            -> ([SCRNASeq S [(Int, FilePath, FilePath)]], [Int])
            -> ReaderT config IO (Maybe (FilePath, FilePath, FilePath))
getSpectral _ ([], _) = return Nothing
getSpectral sampleSize (input, cidx') = do
    dir <- asks ((<> "/Spectral/") . _scrnaseq_output_dir) >>= getPath
    let output1 = dir <> "seed_matrix.bin"
        output2 = dir <> "spectral_1.bin"
        output3 = dir <> "spectral_2.bin"
    liftIO $ do
      gen <- create
      mat <- fmap (MS.fromBlocks 0 . map return) $ runConduit $
          yieldMany (zip sampleCells matFls) .| mapMC (f gen) .| sinkList
      let mat' = fromVector (MS.dim mat) $ MS.flatten mat :: Dynamic D.Matrix V.Vector Double
      B.writeFile output1 $ flat mat'
      withDyn mat' $ \(m@(D.Matrix _) :: Matrix n k Double)  -> do
          let dm = rbf m m :: Matrix n n Double
              left = sing :: Sing 31
              right = (sing :: Sing n) %- (sing :: Sing 2)
          case decideEquality (left %<=? right) STrue of
              Just Refl -> do
                  let Spectral s1 s2 = spectral (sing :: Sing 30) dm
                  B.writeFile output2 $ flat s1
                  B.writeFile output3 $ flat s2
                  return $ Just (output1, output2, output3)
              Nothing -> error ""
  where
    f gen (n, fl) = do
        Dynamic (mat :: Matrix m n Double) <- either (error . show) id . unflat <$> B.readFile fl
        let r = D.rows mat
        if n >= r
          then return $! MS.generate (r, U.length cidx) $ \(i, j) -> mat `D.unsafeIndex` (i, cidx U.! j)
          else do
              ridx <- U.fromList . fromJust <$> sample n [0..r-1] gen
              return $! MS.generate (n, U.length cidx) $ \(i, j) -> mat `D.unsafeIndex` (ridx U.! i, cidx U.! j)
    matFls = input^..folded.replicates._2.files.folded._2
    sampleCells = map (\r -> ceiling $ fromIntegral sampleSize * (fromIntegral r / total)) rs
    rs = input^..folded.replicates._2.files.folded._1
    total = fromIntegral $ foldl1' (+) rs :: Double
    cidx = U.fromList cidx'

{-
nystromExtend :: SCRNASeqConfig config
            -> ([SCRNASeq S [(Int, FilePath, FilePath)]], [Int], Maybe (FilePath, FilePath, FilePath))
            -> ReaderT config IO ()
nystromExtend (_, _, Just (fl1, fl2, fl3)) = do
-}

data Spectral n k = Spectral (SparseMatrix k k Double) (Matrix n k Double)

-- | Perform spectral embedding using random-walk Lapacian.
spectral :: forall k n. (SingI k, SingI n, (k + 1 <=? n - 2) ~ 'True)
         => Sing k
         -> Matrix n n Double    -- ^ Similarity matrix
         -> Spectral n k
spectral k _S@(D.Matrix _) = Spectral
    (S.diag $ G.unsafeFromVector $ V.tail $ G.flatten $ G.map (\x -> 1 / realPart x) eval)
    (G.fromColumns $ tail $ G.toColumns $ G.map realPart evec)
  where
    (eval, evec) = case k' of SNat -> eigS k' _P
    k' = k %+ (SNat :: Sing 1)
    _P = _D' @@ _S
    _D' = S.diag $ G.map recip $ rowSum _S
{-# INLINE spectral #-}

spectralEmbedC :: Monad m
               => Dynamic D.Matrix V.Vector Double  -- ^ Seed matrix used to generate the spectral
               -> Dynamic S.SparseMatrix V.Vector Double
               -> Dynamic D.Matrix V.Vector Double
               -> ConduitT (Dynamic D.Matrix V.Vector Double) 
                           (Dynamic D.Matrix V.Vector Double) m ()
spectralEmbedC sm s1 s2 = withSpectral sm s1 s2 nystromC
{-# INLINE spectralEmbedC #-}

withSpectral :: Dynamic D.Matrix V.Vector Double  -- ^ Seed matrix used to generate the spectral
             -> Dynamic S.SparseMatrix V.Vector Double
             -> Dynamic D.Matrix V.Vector Double
             -> (forall n p k. (SingI n, SingI p, SingI k) => Matrix n p Double -> Spectral n k -> a)
             -> a
withSpectral sm' s1' s2' f = withDyn sm' $ \sm@(D.Matrix _ :: Matrix r1 c1 Double) ->
    withDyn s1' $ \s1@(S.SparseMatrix _ _ _ :: SparseMatrix r2 c2 Double) ->
    withDyn s2' $ \s2@(D.Matrix _ :: Matrix r3 c3 Double) ->
        case decideEquality (sing :: Sing r2) (sing :: Sing c2) of
            Just Refl -> case decideEquality (sing :: Sing c2) (sing :: Sing c3) of
                Just Refl -> case decideEquality (sing :: Sing r1) (sing :: Sing r3) of
                    Just Refl -> f sm $ Spectral s1 s2
                    Nothing -> error ""
                Nothing -> error ""
            Nothing -> error ""
{-# INLINE withSpectral #-}
  
nystromC :: forall n p k monad. (SingI n, SingI p, SingI k, Monad monad)
         => Matrix n p Double   -- ^ Seed matrix
         -> Spectral n k
         -> ConduitT (Dynamic D.Matrix V.Vector Double) 
                     (Dynamic D.Matrix V.Vector Double) monad ()
nystromC mat sp = mapC $ \(Dynamic m@(D.Matrix _ :: Matrix r c Double)) -> 
    case decideEquality (sing :: Sing c) (sing :: Sing p) of
        Just Refl -> Dynamic $ nystrom sp $ rbf m mat
        Nothing -> error ""
{-# INLINE nystromC #-}
      
-- | Nystrom extension.
nystrom :: (SingI k, SingI n, SingI m)
        => Spectral n k
        -> Matrix m n Double   -- ^ Out-of-bag similarity matrix
        -> Matrix m k Double
nystrom (Spectral eval evec) _S = _P @@ evec @@ eval
  where
    _P = _D' @@ _S
    _D' = S.diag $ G.map recip $ rowSum _S
{-# INLINE nystrom #-}

-- | Compute the rbf (gaussian) kernel between X and Y::
-- K(x, y) = exp(-gamma ||x-y||^2)
rbf :: (SingI k, SingI n, SingI m)
    => Matrix n k Double
    -> Matrix m k Double
    -> Matrix n m Double
rbf a b = D.map (\x -> exp $ (negate gamma) * x) $ euclideanSquared a b
  where
    gamma = 1 / fromIntegral (D.cols a)
{-# INLINE rbf #-}

-- | Squared euclidean distance
euclideanSquared :: (SingI k, SingI n, SingI m) => Matrix n k Double -> Matrix m k Double -> Matrix n m Double
euclideanSquared a b = a2 + b2 - D.map (2*) ab
  where
    a2 = rowSum (D.map (\x -> x * x) a) @@ ones
    b2 = ones @@ D.transpose (rowSum (D.map (\x -> x * x) b))
    ab = a @@ D.transpose b
{-# INLINE euclideanSquared #-}