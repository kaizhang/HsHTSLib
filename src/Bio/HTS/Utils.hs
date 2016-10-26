{-# LANGUAGE OverloadedStrings #-}
module Bio.HTS.Utils where

import           Conduit
import           Control.Arrow               (first, second)
import           Control.Monad
import           Control.Monad.Base          ()
import           Control.Monad.State
import qualified Data.ByteString.Char8       as B
import qualified Data.HashMap.Strict         as M
import qualified Data.Vector.Unboxed         as U
import qualified Data.Vector.Unboxed.Mutable as UM
import           Data.Word                   (Word8)
import           System.IO
import           Text.Printf

import           Bio.HTS
import           Bio.HTS.Types

-- | Compute read counts
bamToCounts :: [(B.ByteString, Int)]   -- ^ Chr sizes
            -> Sink Bam HeaderState (M.HashMap B.ByteString (U.Vector (Word8, Word8)))
bamToCounts chrs = do
    counts <- liftBase $ fmap M.fromList $ forM chrs $ \(chr, s) -> do
        v <- UM.replicate s (0,0)
        return (chr, v)
    mapM_C $ \b -> do
        BamHeader hdr <- lift get
        case getChr hdr b >>= flip M.lookup counts of
            Nothing -> return ()
            Just v -> liftBase $ if isRev b
                then UM.modify v (second add1) $ fromIntegral $ endPos b
                else UM.modify v (first add1) $ fromIntegral $ position b
    mapM (liftBase . U.unsafeFreeze) counts
  where
    add1 x | x < 255 = x+1
           | otherwise = x


hCountsToWig :: (Show a, Eq a, Num a, U.Unbox a)
             => Handle
             -> Int   -- ^ step size
             -> M.HashMap B.ByteString (U.Vector a)
             -> IO ()
hCountsToWig h step rcs = do
    _ <- flip M.traverseWithKey rcs $ \chr rc -> do
        hPutStrLn h $ printf "variableStep  chrom=%s span=%d" (B.unpack chr) step
        flip U.imapM_ rc $ \i x -> unless (x == 0) $
            hPutStrLn h $ (show $ i+1) ++ " " ++ show x
    return ()
