{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE PartialTypeSignatures #-}
module Bio.HTS.Utils
    ( -- * Mark duplicates
      markDupBy
    , Orientation(..)
    , BAMKey(..)
    , makeKey
    , makeKeySingle
    , makeKeyPair
    
      -- * Other utilities
    , fragmentSizeDistr
    ) where

import           Conduit
import qualified Data.Map.Strict as M
import Data.List
import Data.Maybe
import qualified Data.Sequence as S
import qualified Data.ByteString.Char8 as B
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM

import           Bio.HTS.BAM
import           Bio.HTS.Types

data Orientation = F | R | FF | FR | RR | RF deriving (Eq, Ord, Show)

data BAMKey = Single { _ref_id1 :: Int
                     , _loc1 :: Int
                     , _orientation :: Orientation
                     , _barcode :: Maybe B.ByteString }
            | Pair { _ref_id1 :: Int  -- ^ ref id of this tag
                   , _ref_id2 :: Int  -- ^ ref id of the paired tag
                   , _loc1 :: Int     -- ^ location of this tag
                   , _loc2 :: Int     -- ^ location of the paired tag
                   , _orientation :: Orientation
                   , _leftmost :: Bool   -- ^ Is this tag leftmost
                   , _barcode :: Maybe B.ByteString }
            deriving (Eq, Ord, Show)

-- | Generate fingerprint for single-end reads.
makeKeySingle :: (BAM -> Maybe B.ByteString)   -- ^ Barcode extraction function
              -> BAM
              -> BAMKey
makeKeySingle fn bam = Single ref1 (if isFwd1 then lloc1 else rloc1)
    (if isFwd1 then F else R) $ fn bam
  where
    ref1 = refId bam
    isFwd1 = not $ isRC flg
    lloc1 = startLoc bam - fst clipped1 + 1
    rloc1 = endLoc bam + snd clipped1
    flg = flag bam
    clipped1 = getClipped $ fromJust $ cigar bam
{-# INLINE makeKeySingle #-}

-- | Create a pair of keys (single-end and paired-end).
makeKeyPair :: (BAM -> Maybe B.ByteString)   -- ^ Barcode extraction function
            -> (BAM, BAM)
            -> BAMKey
makeKeyPair fn (bam1, bam2) = Pair ref1 ref2 loc1 loc2 orientation isLeftMost $ fn bam1
  where
    ref1 = refId bam1
    lloc1 = startLoc bam1 - fst clipped1 + 1
    rloc1 = endLoc bam1 + snd clipped1
    clipped1 = getClipped $ fromJust $ cigar bam1
    isFwd1 = not $ isRC $ flag bam1

    ref2 = mateRefId bam2
    lloc2 = startLoc bam2 - fst clipped2 + 1
    rloc2 = endLoc bam2 + snd clipped2
    clipped2 = getClipped $ fromJust $ cigar bam2
    isFwd2 = not $ isRC $ flag bam2

    isLeftMost
        | ref1 /= ref2 = ref1 < ref2
        | otherwise = if isFwd1 == isFwd2
            then if isFwd1 then lloc1 <= lloc2 else rloc1 <= rloc2
            else if isFwd1 then lloc1 <= rloc2 else rloc1 <= lloc2
    orientation
        | isLeftMost = if isFwd1 == isFwd2
            then if isFwd1
                then if isFirstSegment flg then FF else RR
                else if isFirstSegment flg then RR else FF
            else if isFwd1 then FR else RF
        | otherwise = if isFwd1 == isFwd2
            then if isFwd1
                then if isFirstSegment flg then RR else FF
                else if isFirstSegment flg then FF else RR
            else if isFwd1 then RF else FR
    loc1 | isFwd1 == isFwd2 = if isLeftMost then lloc1 else rloc1
         | otherwise = if isFwd1 then lloc1 else rloc1
    loc2 | isFwd1 == isFwd2 = if isLeftMost then rloc2 else lloc2
         | otherwise = if isFwd1 then rloc2 else lloc2
    flg = flag bam1


-- | Create a pair of keys (single-end and paired-end).
makeKey :: (BAM -> Maybe B.ByteString)   -- ^ Barcode extraction function
        -> BAM
        -> (BAMKey, BAMKey)
makeKey fn bam = (single, pair)
  where
    single = Single ref1 (if isFwd1 then lloc1 else rloc1)
        (if isFwd1 then F else R) bc
    pair = Pair ref1 ref2 loc1 loc2 orientation isLeftMost bc
    bc = fn bam
    ref1 = refId bam
    lloc1 = startLoc bam - fst clipped1 + 1
    rloc1 = endLoc bam + snd clipped1
    clipped1 = getClipped $ fromJust $ cigar bam
    isFwd1 = not $ isRC flg
    ref2 = mateRefId bam
    lloc2 = mateStartLoc bam - fst clipped2 + 1
    rloc2 = mateStartLoc bam + ciglen cig + snd clipped2
    clipped2 = getClipped cig
    isFwd2 = not $ isNextRC flg
    cig = case queryAuxData ('M', 'C') bam of
        Just (AuxString x) -> string2Cigar x
        _ -> error "No MC tag. Please run samtools fixmate on file first."
    isLeftMost
        | ref1 /= ref2 = ref1 < ref2
        | otherwise = if isFwd1 == isFwd2
            then if isFwd1 then lloc1 <= lloc2 else rloc1 <= rloc2
            else if isFwd1 then lloc1 <= rloc2 else rloc1 <= lloc2
    orientation
        | isLeftMost = if isFwd1 == isFwd2
            then if isFwd1
                then if isFirstSegment flg then FF else RR
                else if isFirstSegment flg then RR else FF
            else if isFwd1 then FR else RF
        | otherwise = if isFwd1 == isFwd2
            then if isFwd1
                then if isFirstSegment flg then RR else FF
                else if isFirstSegment flg then FF else RR
            else if isFwd1 then RF else FR
    loc1 | isFwd1 == isFwd2 = if isLeftMost then lloc1 else rloc1
         | otherwise = if isFwd1 then lloc1 else rloc1
    loc2 | isFwd1 == isFwd2 = if isLeftMost then rloc2 else lloc2
         | otherwise = if isFwd1 then rloc2 else lloc2
    flg = flag bam
    ciglen (CIGAR c) = foldl' f 0 c
      where f acc (n,x) = if x `elem` "MDN=X" then n + acc else acc

-- | Get the number of clips at both ends.
getClipped :: CIGAR -> (Int, Int)
getClipped (CIGAR c) | length c <= 1 = (0,0)
                     | otherwise = (getC $ head c, getC $ last c)
  where
    getC x = case x of
        (n, 'S') -> n
        (n, 'H') -> n
        _ -> 0
{-# INLINE getClipped #-}


-- | Remove duplicated reads. Duplicates are determined by
-- checking for matching keys. The Key is comprised of:
--
-- 1. Chromosome
-- 2. Orientation (forward/reverse)
-- 3. Unclipped Start(forward)/End(reverse)
-- 4. Barcode
--
-- Keep the read that has a higher base quality sum (sum of all
-- base qualities in the record above 15).
markDupBy :: MonadIO m
          => (BAM -> Maybe B.ByteString)   -- ^ Barcode extraction function, if any.
          -> ConduitT BAM BAM m ()
markDupBy bcFn = go (-1,-1) M.empty S.empty
  where
    go prev keyMap readBuf = await >>= \case
        Nothing -> mapM_ (\(_,_,x) -> yield x) readBuf
        Just bam -> markDup prev keyMap readBuf bam >>=
            (\(a,b,c) -> go a b c)
    markDup prev keyMap readBuf bam
        | isUnmapped flg = yield bam >> return (prev, keyMap, readBuf)
        | not (isSorted prev cur) = error $
            "bad coordinate order: " ++ show (prev, cur)
        | isDup flg = return (cur, keyMap, readBuf')
        | otherwise = do
            keyMap' <- addAndMark key bam keyMap
            update keyMap' readBuf' cur
      where
        readBuf' = readBuf S.|> (key, _loc1 singleKey, bam)
        cur = (refId bam, startLoc bam)
        key | not (hasMultiSegments flg) || isNextUnmapped flg = singleKey
            | otherwise = pairKey
        (singleKey, pairKey) = makeKey bcFn bam
        flg = flag bam
    update keyMap readBuf (chr, loc) = mapM_ (\(_,_,b) -> yield b) exclude >>
        return ((chr, loc), keyMap', kept)
      where
        keyMap' = foldl' (\m (k,_,_) -> M.delete k m) keyMap exclude
        (exclude, kept) = S.breakl toKeep readBuf
        toKeep (key, pos,_) = _ref_id1 key == chr && pos + max_len > loc
    addAndMark key bam keyMap = liftIO $ case M.lookup key keyMap of
        Nothing -> return $ M.insert key (bam, sc) keyMap
        Just (old, old_sc) -> if sc > old_sc -- || (sc == old_sc && queryName bam < queryName old)
            then setDup old >> return (M.insert key (bam, sc) keyMap)
            else setDup bam >> return keyMap
      where
        sc = case queryAuxData ('m', 's') bam of
            Just (AuxInt x) -> x + fromJust (sumQual 15 bam)
            _ -> fromJust $ sumQual 15 bam
    max_len = 500
    isSorted (chr1, loc1) (chr2, loc2) = chr1 < chr2 ||
        (chr1 == chr2 && loc1 <= loc2)
{-# INLINE markDupBy #-}


-- | Compute fragment size distribution from paired end BAM records.
fragmentSizeDistr :: PrimMonad m
                  => Int    -- ^ Largest fragment size
                  -> ConduitT BAM o m (U.Vector Double)
fragmentSizeDistr n = do
    vec <- lift $ UM.replicate n 0
    mapM_C $ f vec
    vec' <- lift $ U.unsafeFreeze vec
    return $ U.map (/ (U.sum vec')) vec'
  where
    f v x | s >= n = return ()
          | otherwise = UM.modify v (+1) s
      where
        s = abs $ tLen x
{-# INLINE fragmentSizeDistr #-}