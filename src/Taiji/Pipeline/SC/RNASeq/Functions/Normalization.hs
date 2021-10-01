{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.SC.RNASeq.Functions.Normalization
    ( fitNB
    , normalization
    , selectFeatures
    , logNormalize
    ) where 

import qualified Data.ByteString.Char8 as B
import Data.Conduit.Internal (zipSources)
import Bio.Utils.Functions (gaussianKDE)
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Storable as SV
import qualified Data.Text as T
import Control.Arrow (second)
import Shelly (shelly, run_, mkdir_p)
import Numeric.Sampling (psample)
import System.Random.MWC
import Language.Javascript.JMacro
import Data.Conduit.List (chunksOf)

import qualified Data.Matrix.Static.Dense as D
import qualified Data.Matrix.Static.Sparse as S
import Data.Matrix.Dynamic (fromTriplet, withDyn, Dynamic(..))
import Data.Matrix.Static.LinearAlgebra hiding (colSum)
import qualified Data.Matrix.Static.LinearAlgebra as L
import Data.Singletons.Prelude hiding ((@@))
import Flat (flat, unflat)

import           Taiji.Pipeline.SC.RNASeq.Types (figDir, SCRNASeqConfig(..), SCRNASeq)
import Taiji.Prelude
import Taiji.Utils
import Taiji.Utils.Plot
import Taiji.Utils.Plot.ECharts

logNormalize :: (Elem 'Gzip tags ~ 'True, SCRNASeqConfig config)
             => SCRNASeq S (File '[ColumnName, Gzip] 'Tsv, File tags 'Matrix)
             -> ReaderT config IO (SCRNASeq S (File '[ColumnName, Gzip] 'Tsv, File tags 'Matrix))
logNormalize input = do
    dir <- asks ((<> "/Quantification/Normalization/") . _scrnaseq_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_log_norm.mat.gz" dir (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \(r, fl) -> do
        fmap (mapRows (second f)) (mkSpMatrix readDouble $ fl^.location) >>=
            saveMatrix output toShortest
        return (r, location .~ output $ fl)
        )
  where
    f xs = let s = foldl1' (+) (map snd xs) / 10000
           in map (second (\x -> x / s)) xs

fitNB :: (Elem 'Gzip tags ~ 'True, SCRNASeqConfig config)
      => SCRNASeq S (File '[ColumnName, Gzip] 'Tsv, File tags 'Matrix)
      -> ReaderT config IO (SCRNASeq S (FilePath, FilePath))
fitNB input = do
    dir <- asks ((<> "/Quantification/Normalization/") . _scrnaseq_output_dir) >>= getPath
    let output = outdir <> "/model_parameters.txt"
        outputGeneMean = outdir <> "/log_gene_mean.txt"
        outdir = printf "%s/%s_rep%d" dir (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traversed.files %%~ liftIO . ( \(r, fl) -> withTempDir Nothing $ \tmpdir -> do
        let tmpMat = tmpdir <> "/tmp.mat.gz"
            tmpRow = tmpdir <> "/row.txt"
            tmpCol = tmpdir <> "/col.txt"
        shelly $ mkdir_p outdir
        gen <- create

        cellByGene <- mkSpMatrix readInt $ fl^.location
        let totalCell = _num_row cellByGene
            totalGene = _num_col cellByGene
            nGene = 2000
            nCell = 30000

        -- Compute geometric mean for each gene
        logGeneMean <- U.map (\x -> logBase 10 $ exp (x / fromIntegral totalCell) - 1) <$>
            colSum (fmap (\x -> log $ fromIntegral x + 1) cellByGene)
        (col, mat) <- if totalGene > nGene
            then do
                geneIdx <- sort <$> kdeSampling gen nGene logGeneMean
                return (map (logGeneMean U.!) geneIdx, selectCols geneIdx cellByGene)
            else
                return (U.toList logGeneMean, cellByGene)

        -- Compute log10(totalReads) for each cell
        cellReads <- runResourceT $ runConduit $ streamRows cellByGene .|
            mapC (logBase 10 . fromIntegral . foldl1' (+) . map snd . snd) .|
            sinkVector
        (row, cellByGene') <- if totalCell > nCell
            then do
                cellIdx <- sort <$> kdeSampling gen nCell cellReads
                return (map (cellReads U.!) cellIdx, selectRows cellIdx mat)
            else 
                return (U.toList cellReads, mat)

        B.writeFile tmpRow $ B.unlines $ map toShortest row
        B.writeFile tmpCol $ B.unlines $ map toShortest col
        B.writeFile outputGeneMean $ B.unlines $ map toShortest $ U.toList logGeneMean
        saveMatrix tmpMat (fromJust . packDecimal) cellByGene'
        shelly $ run_ "taiji-utils" $ map T.pack
            [ "normalize", tmpMat, tmpCol, tmpRow, outputGeneMean, output
            , "--plot-dir", outdir ]
        return (output, outputGeneMean)
        )

normalization :: SCRNASeqConfig config
              => SCRNASeq S ((FilePath, FilePath), (Int, Int), File tags 'Matrix)
              -> ReaderT config IO (SCRNASeq S [(Int, FilePath, FilePath)])
normalization input = do
    dir <- asks ((<> "/Quantification/Normalization/") . _scrnaseq_output_dir) >>= getPath
    input & replicates.traversed.files %%~ liftIO . ( \((params, _), (i,j), fl) -> do
        [beta_0, beta_1, alpha] <- map (SV.fromList . map readDouble . B.words) . B.lines <$>
            B.readFile params
        mat <- mkSpMatrix readDouble $ fl^.location
        let f (b, Dynamic m@(D.Matrix _)) = do
                let output = printf "%s/normalized_gene_count_%d-%d_%d.bin" dir i j b
                    outputSummary = printf "%s/normalized_gene_count_summary_%d-%d_%d.bin" dir i j b
                let [s] = D.toRows $ L.colSum m
                    [s_2] = D.toRows $ L.colSum $ D.map (\x -> x * x) m
                liftIO $ do
                    B.writeFile output $ flat m
                    B.writeFile outputSummary $ flat $ U.zip (U.convert s) (U.convert s_2)
                return (D.rows m, output, outputSummary)
            batchSize = 2000
        runResourceT $ runConduit $
            zipSources (iterateC succ (0 :: Int)) (predictBatch i j batchSize mat beta_0 beta_1 alpha) .|
            mapMC f .| sinkList
        )

predictBatch :: Int -> Int -> Int
             -> SpMatrix Double
             -> SV.Vector Double    -- ^ beta_0
             -> SV.Vector Double    -- ^ beta_1
             -> SV.Vector Double    -- ^ alpha
             -> ConduitT () (Dynamic D.Matrix SV.Vector Double) (ResourceT IO) ()
predictBatch i j bs mat beta0' beta1' alpha' = readSparseMatrixChunk mat i j bs .| mapC f
  where
    f (Dynamic (sp@(S.SparseMatrix _ _ _) :: S.SparseMatrix c g SV.Vector Double)) = do
        let beta0 = D.fromVector beta0' :: Matrix 1 g Double
            beta1 = D.fromVector beta1' :: Matrix 1 g Double
            alpha = D.fromVector alpha' :: Matrix 1 g Double
            mu = expectedCount beta0 beta1 $ D.map (logBase 10) $ rowSum sp
            sigma = expectedStd alpha mu
         in Dynamic $ D.map (max (negate clipValue) . min clipValue) $
                pearsonResidual mu sigma sp
    clipValue = sqrt 30000
{-# INLINE predictBatch #-}

readSparseMatrixChunk :: SpMatrix Double
                      -> Int   -- ^ Starting row index, inclusive
                      -> Int   -- ^ End row index, exclusive
                      -> Int   -- ^ Chunk size
                      -> ConduitT () (Dynamic S.SparseMatrix SV.Vector Double) (ResourceT IO) ()
readSparseMatrixChunk mat s e size = streamRows mat .| (dropC s >> takeC n) .| chunksOf size .| mapC f
  where
    f chunk = fromTriplet (length chunk, _num_col mat) $ U.fromList $ concat $ zipWith g [0..] chunk
      where
        g i (_, xs) = map (\(j,x) -> (i,j,x)) xs
    n = e - s
{-# INLINE readSparseMatrixChunk #-}

-- | mu
expectedCount ::(SingI c, SingI g) 
              => Matrix 1 g Double   -- ^ beta_0
              -> Matrix 1 g Double   -- ^ beta_1
              -> Matrix c 1 Double   -- ^ log10 of total reads of the cell
              -> Matrix c g Double
expectedCount beta0 beta1 lg_m = D.map exp $ (ones @@ beta0) + (lg_m @@ beta1)
{-# INLINE expectedCount #-}

-- | sigma
expectedStd :: (SingI c, SingI g)
            => Matrix 1 g Double   -- ^ alpha
            -> Matrix c g Double
            -> Matrix c g Double
expectedStd alpha mu = D.map sqrt $ mu + ((ones @@ alpha) * D.map (\x -> x * x) mu)
{-# INLINE expectedStd #-}

pearsonResidual :: (SingI c, SingI g)
                => Matrix c g Double
                -> Matrix c g Double
                -> SparseMatrix c g Double
                -> Matrix c g Double
pearsonResidual mu sigma x = (x %-% mu) / sigma 
{-# INLINE pearsonResidual #-}

kdeSampling :: PrimMonad m
            => Gen (PrimState m)
            -> Int   -- ^ Number of samples
            -> U.Vector Double  -- ^ Input
            -> m [Int]   -- ^ index
kdeSampling gen n xs = fromJust <$>
    psample n (zip (U.toList $ normalize $ U.map (recip . kde) xs) [0..]) gen
  where
    kde = gaussianKDE 10000 xs
    normalize vs = let s = U.foldl1 (+) vs in U.map (/s) vs
{-# INLINE kdeSampling #-}

selectFeatures :: SCRNASeqConfig config
               => Int   -- ^ Number of features
               -> ( Maybe (SCRNASeq S (File '[ColumnName, Gzip] 'Tsv, File tags 'Matrix))
                  , [SCRNASeq S [(Int, FilePath, FilePath)]]
                  , Maybe (SCRNASeq S (FilePath, FilePath)) )
               -> ReaderT config IO [Int]
selectFeatures n (Just input1, input2, Just input3) = do
    dir <- figDir
    let output = dir <> "/variable_gene_selection.html"
    liftIO $ do
        genes <- fmap (map (head . words) . lines) $ readFile $
            input1^.replicates._2.files._1.location
        mat <- mkSpMatrix id $ input1^.replicates._2.files._2.location
        geneVar <- map snd . U.toList <$> meanVarC
            (input2^..folded.replicates._2.files.folded._3) (_num_row mat) (_num_col mat)
        geneMean <- fmap (map ((\x -> 10**x) . readDouble) . B.lines) $ B.readFile $ input3^.replicates._2.files._2

        let (selected, remaining) = splitAt n $
                sortBy (flip (comparing (snd . snd))) $ zip [0..] $ zip geneMean geneVar
            plt = toolbox >+> attr >+> scatter'
                [ ("Selected", map snd selected)
                , ("Unused", map snd remaining) ]
            attr = [jmacroE| {
                grid: {
                    height: 400,
                    width: 400
                },
                xAxis: {
                    type: "log",
                    axisTick: {show:true},
                    axisLabel: {show:true},
                    splitLine: {show:true},
                    name: "Gene mean"
                },
                yAxis: {
                    type: "log",
                    axisTick: {show:true},
                    axisLabel: {show:true},
                    splitLine: {show:true},
                    name: "Residue variance"
                } }|]
 
        savePlots output [] [plt]

        return $ map fst selected
selectFeatures _ _ = return []

meanVarC :: [FilePath] -> Int -> Int -> IO (U.Vector (Double, Double))
meanVarC fls r c = runConduit $ fmap (U.map g) $ yieldMany fls .| foldMC f (U.replicate c (0, 0))
  where
    g (m, v) = let m' = m / fromIntegral r
                in (m', v / fromIntegral r - m' * m')
    f vec fl = do
        vec' <- either (error . show) id . unflat <$> B.readFile fl
        return $ U.zipWith (\(a,b) (a',b') -> (a+a', b+b')) vec' vec