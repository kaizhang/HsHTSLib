{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE PartialTypeSignatures #-}
module Bio.HTS.Utils
    ( Orientation(..)
    , BAMKey(..)
    , markDupBy
    ) where

import           Conduit
import qualified Data.Map.Strict as M
import Data.List
import Data.Maybe
import qualified Data.Sequence as S

import           Bio.HTS.BAM
import           Bio.HTS.Types

data Orientation = F | R | FF | FR | RR | RF deriving (Eq, Ord, Show)

data BAMKey = Single { _ref_id1 :: Int
                     , _loc1 :: Int
                     , _orientation :: Orientation }
            | Pair { _ref_id1 :: Int
                   , _ref_id2 :: Int
                   , _loc1 :: Int
                   , _loc2 :: Int
                   , _orientation :: Orientation
                   , _leftmost :: Bool
                   }
            deriving (Eq, Ord, Show)

makeKey :: BAM -> (BAMKey, BAMKey)
makeKey bam = (single, pair)
  where
    single = Single ref1 (if isFwd1 then lloc1 else rloc1)
        (if isFwd1 then F else R)
    pair = Pair ref1 ref2 loc1 loc2 orientation isLeftMost
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

-- Duplicates are determined by checking for matching keys.
-- The Key is comprised of:
-- Chromosome
-- Orientation (forward/reverse)
-- Unclipped Start(forward)/End(reverse)
-- Library
-- Rules:
-- Skip Unmapped Reads, they are not marked as duplicate
-- Reads whose mate is unmapped are treated as single-end
-- Mark a Single-End Read Duplicate (or remove it if configured to do so) if:
-- A paired-end record has the same key (even if the pair is not proper/the mate is not found)
-- -OR-
-- A single-end record has the same key and a higher base quality sum (sum of all base qualities in the record above --minBaseQual)
-- Mark both Paired-End Reads Duplicate if:


-- | Remove duplicated reads according to comparison function.
-- Keep the read that has a higher base quality sum (sum of all
-- base qualities in the record above 15).
markDupBy :: MonadIO m
        => (BAM -> Int)
        -> ConduitT BAM BAM m ()
markDupBy scFn = go (-1,-1) M.empty S.empty
  where
    go prev keyMap readBuf = await >>= \case
        Nothing -> mapM_ (yield . fst) $ M.elems keyMap
        Just bam -> markDup prev keyMap readBuf bam >>=
            (\(a,b,c) -> go a b c)
    markDup prev keyMap readBuf bam
        | isUnmapped flg = yield bam >> return (prev, keyMap, readBuf)
        | not (isSorted prev cur) = error $
            "bad coordinate order: " ++ show (prev, cur)
        | isDup flg = return (cur, keyMap, readBuf')
        | otherwise = do
            keyMap' <- addBamKey key bam keyMap
            update keyMap' readBuf' cur
      where
        readBuf' = readBuf S.|> (key, _loc1 singleKey, bam)
        cur = (refId bam, startLoc bam)
        key | not (hasMultiSegments flg) || isNextUnmapped flg = singleKey
            | otherwise = pairKey
        (singleKey, pairKey) = makeKey bam
        flg = flag bam
    update keyMap readBuf (chr, loc) = mapM_ (\(_,_,b) -> yield b) exclude >>
        return ((chr, loc), keyMap', kept)
      where
        keyMap' = foldl' (\m (k,_,_) -> M.delete k m) keyMap exclude
        (exclude, kept) = S.breakl toKeep readBuf
        toKeep (key, pos,_) = _ref_id1 key == chr && pos + max_len > loc
    addBamKey key bam keyMap = case M.insertLookupWithKey cmp key (bam, sc) keyMap of
        (Nothing, m) -> return m
        (Just old, m) -> do
            let b = if sc > snd old then fst old else bam
            liftIO $ setDup b
            return m
      where
        sc = scFn bam
        cmp _ new old = if snd new > snd old then new else old
    max_len = 500
    isSorted (chr1, loc1) (chr2, loc2) = chr1 < chr2 ||
        (chr1 == chr2 && loc1 <= loc2)
{-# INLINE markDupBy #-}
