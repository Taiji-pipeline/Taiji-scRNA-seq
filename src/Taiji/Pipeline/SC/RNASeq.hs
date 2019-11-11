{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.SC.RNASeq ( builder) where

import           Control.Workflow
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B

import Taiji.Prelude
import Taiji.Utils
import           Taiji.Pipeline.SC.RNASeq.Functions

builder :: Builder ()
builder = do
    node "Read_Input" 'readInput $ return ()

    node "Get_Fastq" [| return . getFastq |] $ return ()
    nodePar "Demultiplex" 'extractBarcode $ return ()
    path ["Read_Input", "Get_Fastq", "Demultiplex"]

    node "Make_Index" 'mkIndex $ return ()
    nodePar "Align" 'tagAlign $ do
        nCore .= 8
        memory .= 50
    nodePar "Filter_Bam" 'filterNameSortBam $ return ()
    nodePar "Quantification" 'quantification $ memory .= 8
    nodePar "Remove_Doublet" 'removeDoublet $ return ()
    path ["Demultiplex", "Make_Index", "Align", "Filter_Bam"
        , "Quantification", "Remove_Doublet"]

    {-
    nodePar "Reduce_Dimension" [| \input -> do
        let input' = input & replicates.traverse.files %~ snd . fst
            prefix = "/Cluster/"
        filterMatrix prefix input' >>= spectral prefix Nothing 
        |] $ return ()
    nodePar "Cluster" [| \input -> clustering "/Cluster/" $
        input & replicates.traverse.files %~ return
        |] $ return ()
    path ["Remove_Doublet", "Reduce_Dimension", "Cluster"]
    -}

    node "Merge_Matrix" [| \mats -> if length mats < 1
        then return Nothing
        else do
            dir <- asks _scrnaseq_output_dir >>=
                getPath . (<> (asDir "/Cluster/"))
            let output = dir <> "Merged.mat.gz"
            liftIO $ concatMatrix output $ flip map mats $ \mat ->
                ( Just $ B.pack $ T.unpack $ mat^.eid
                , mat^.replicates._2.files._1._2.location )
            return $ Just $ (head mats & eid .~ "Merged") &
                replicates._2.files %~ (\((_,fl),_) -> location .~ output $ fl)
            |] $ return ()
    node "Merged_Reduce_Dimension" [| \case
        Nothing -> return Nothing
        Just input -> do
            let prefix = "/Cluster/"
            fmap Just $ filterMatrix prefix input >>= spectral prefix Nothing
        |] $ return ()
    node "Merged_Cluster" [| \case
        Nothing -> return Nothing
        Just input -> fmap Just $ clustering "/Cluster/" $
            input & replicates.traverse.files %~ return
        |] $ return ()
    path ["Remove_Doublet", "Merge_Matrix", "Merged_Reduce_Dimension", "Merged_Cluster"]

    node "Extract_Sub_Matrix" [| \(mats, cl) -> case cl of
        Nothing -> return []
        Just clFl -> do
            let mats' = flip map mats $ \x ->
                    x & replicates.traverse.files %~ snd . fst
            subMatrix "/temp/" mats' $ clFl^.replicates._2.files
        |] $ return ()
    node "Make_Expr_Table" [| mkExprTable "/Quantification/" |] $ return ()
    node "Marker_Gene" [| \case
        Nothing -> return Nothing
        Just x -> fmap Just $ specificGene "/Marker/" x
        |] $ return ()
    ["Remove_Doublet", "Merged_Cluster"] ~> "Extract_Sub_Matrix"
    path ["Extract_Sub_Matrix", "Make_Expr_Table", "Marker_Gene"]

    node "QC" 'plotQC $ return ()
    ["Quantification"] ~> "QC"
