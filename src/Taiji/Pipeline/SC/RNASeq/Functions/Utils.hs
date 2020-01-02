{-# LANGUAGE OverloadedStrings #-}
module Taiji.Pipeline.SC.RNASeq.Functions.Utils
    ( groupOnC
    , getIndex
    ) where

import qualified Data.ByteString.Char8                as B
import Conduit

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