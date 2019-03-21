{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}

module Bio.HTS.BAM
    ( getBamHeader
    , getSortOrder
    , streamBam
    , sinkBam

      -- * Field accessors
    , refId
    , refName
    , startLoc
    , endLoc
    , readLen
    , isRev
    , flag
    , mapq
    , getSeq
    , queryName
    , qualityS
    , quality
    , sumQual
    , cigar
    , mateRefId
    , mateRefName
    , mateStartLoc
    , tLen
    , auxData
    , queryAuxData
    , bamToSam

      -- * Flag interpretation
    , hasMultiSegments
    , isProperAligned
    , isUnmapped
    , isNextUnmapped
    , isRC
    , isNextRC
    , isFirstSegment
    , isLastSegment
    , isSecondary
    , isBadQual
    , isDup
    , isSupplementary

      -- * Modify BAM
    , appendAux
    , setDup
    ) where

import           Conduit
import           Data.Bits                (testBit)
import qualified Data.ByteString.Char8    as B
import qualified Data.ByteString as BS
import           Data.Int
import           Data.Word
import           Foreign.C.String
import           Foreign.ForeignPtr
import           Foreign.Marshal.Array
import Foreign.Marshal.Alloc
import Foreign.Storable (poke, peekByteOff)
import           Foreign.Ptr
import           System.IO
import           System.IO.Unsafe         (unsafePerformIO)

import           Bio.HTS.Internal
import           Bio.HTS.Types

#include "htslib/sam.h"

getBamHeader :: FilePath -> IO BAMHeader
getBamHeader input = withHTSFile input ReadMode $ \hts ->
    getBgzf hts >>= bamHdrRead >>= newForeignPtr bamHdrDestroy >>=
    return . BAMHeader
{-# INLINE getBamHeader #-}
foreign import ccall unsafe "&bam_hdr_destroy"
    bamHdrDestroy :: FunPtr (Ptr BamHdr -> IO ())

-- | Get the sort information.
getSortOrder :: BAMHeader -> SortOrder
getSortOrder header = case lookup "SO" fields of
    Just "unknown" -> Unknown
    Just "unsorted" -> Unsorted
    Just "queryname" -> Queryname
    Just "coordinate" -> Coordinate
    _ -> Unknown
  where
    fields = map (f . B.split ':') $ tail $ B.split '\t' $ head $ B.lines $
        showBamHeader header
    f [a,b] = (a,b)
    f _ = error "Auxiliary field parsing failed!"

-- | Create a Bam stream from the file.
streamBam :: MonadResource m => FilePath -> ConduitT i BAM m ()
streamBam input = bracketP (htsOpen input "r") htsClose $ \hts -> do
    bgzf <- liftIO $ getBgzf hts
    _ <- liftIO $ bamHdrRead bgzf >>= newForeignPtr bamHdrDestroy
    loop bgzf
  where
    loop bgzf = do
        ptr <- liftIO $ callocBytes {# sizeof bam1_t #}
        code <- liftIO $ bamRead1 bgzf ptr
        case () of
            _ | code > 0 -> liftIO (newForeignPtr bamDestory1 ptr) >>=
                    yield . BAM >> loop bgzf
              | code == -1 -> return ()
              | code == -2 -> error "truncated file"
              | otherwise -> error $ "read bam failed with code: " ++ show code
{-# INLINE streamBam #-}

foreign import ccall unsafe "&bam_destroy1"
    bamDestory1 :: FunPtr (Ptr Bam1 -> IO ())

-- | Write Bam records to a file.
sinkBam :: MonadResource m => FilePath -> BAMHeader -> ConduitT BAM o m ()
sinkBam output header = bracketP (htsOpen output "wb") htsClose $ \hts -> do
    bgzf <- liftIO $ getBgzf hts
    liftIO (withForeignPtr (unbamHeader header) $ bamHdrWrite bgzf) >>= \case
        0 -> mapM_C $ \bam -> liftIO $ do
            code <- withForeignPtr (unbam bam) (bamWrite1 bgzf)
            if code < 0
                then error "bam_write1 failed"
                else return ()
        _ -> error "'bam_hdr_write' failed."


-- | Reference sequence ID, −1 <= refId < n_ref; -1 for a read
-- without a mapping position.
refId :: BAM -> Int
refId = unsafePerformIO . flip withForeignPtr fun . unbam
  where
    fun = fmap fromIntegral . {#get bam1_t->core.tid #}
{-# INLINE refId #-}

-- | Return the reference sequence name (chromosome name) given the
-- bam file header.
refName :: BAMHeader -> BAM -> Maybe B.ByteString
refName header bam
    | chr < 0 = Nothing
    | otherwise = Just $ unsafePerformIO $
        withForeignPtr (unbamHeader header) $ \h ->
            bamChr h chr >>= B.packCString
  where
    chr = refId bam
{-# INLINE refName #-}

-- | Return the 0-based leftmost coordinate.
startLoc :: BAM -> Int
startLoc = unsafePerformIO . flip withForeignPtr fun . unbam
  where
    fun = fmap fromIntegral . {#get bam1_t->core.pos #}
{-# INLINE startLoc #-}

-- | Return the end location according to a 0-based coordinate system.
-- For a mapped read, this is just position + cigar2rlen.
-- For an unmapped read (either according to its flags or if it has no cigar
-- string), we return position + 1 by convention.
endLoc :: BAM -> Int
endLoc = unsafePerformIO . flip withForeignPtr bamEndpos .unbam
{-# INLINE endLoc #-}

-- | Return the query length (read length).
readLen :: BAM -> Int
readLen = unsafePerformIO . flip withForeignPtr fun . unbam
  where
    fun = fmap fromIntegral . {#get bam1_t->core.l_qseq #}
{-# INLINE readLen #-}

-- | Whether the query is on the reverse strand.
isRev :: BAM -> Bool
isRev = unsafePerformIO . flip withForeignPtr bamIsRev . unbam
{-# INLINE isRev #-}

-- | Return the flag.
flag :: BAM -> Word16
flag = unsafePerformIO . flip withForeignPtr fn . unbam
  where
    fn = fmap fromIntegral . {#get bam1_t->core.flag #}
{-# INLINE flag #-}

-- | MAPping Quality. It equals −10 log10 Pr{mapping position is wrong},
-- rounded to the nearest integer. A value 255 indicates that the
-- mapping quality is not available.
mapq :: BAM -> Word8
mapq = unsafePerformIO . flip withForeignPtr fn . unbam
  where
    fn = fmap fromIntegral . {#get bam1_t->core.qual #}
{-# INLINE mapq #-}

-- | Return the DNA sequence.
getSeq :: BAM -> Maybe B.ByteString
getSeq = unsafePerformIO . flip withForeignPtr fn . unbam
  where
    fn b = fromIntegral <$> {#get bam1_t->core.l_qseq #} b >>= \case
        0 -> return Nothing
        n -> allocaArray n $ \str -> do
            bamGetSeq b str n
            Just <$> B.packCStringLen (str, n)
{-# INLINE getSeq #-}

-- | Get the name of the query.
queryName :: BAM -> B.ByteString
queryName = unsafePerformIO . flip withForeignPtr fn . unbam
  where
    fn b = {#get bam1_t->data #} b >>= B.packCString . castPtr
{-# INLINE queryName #-}

-- | Human readable quality score which is: Phred base quality + 33.
qualityS :: BAM -> Maybe B.ByteString
qualityS = fmap (BS.map (+33)) . quality

-- | Phred base quality (a sequence of 0xFF if absent).
quality :: BAM -> Maybe B.ByteString
quality = unsafePerformIO . flip withForeignPtr fn . unbam
  where
    fn b = fromIntegral <$> {#get bam1_t->core.l_qseq #} b >>= \case
        0 -> return Nothing
        n -> allocaArray n $ \str -> bamGetQual b str n >>= \case
            0 -> Just <$> B.packCStringLen (str, n)
            _ -> return Nothing
{-# INLINE quality #-}

-- | Sum of quality scores that above certain threshold (e.g., 15).
sumQual :: Int  -- ^ Threshold
        -> BAM
        -> Maybe Int
sumQual th bam = BS.foldl' f 0 <$> quality bam 
  where
    f acc x | fromIntegral x < th = acc
            | otherwise = fromIntegral x + acc
{-# INLINE sumQual #-}

-- | Get the CIGAR string.
cigar :: BAM -> Maybe CIGAR
cigar = unsafePerformIO . flip withForeignPtr fn . unbam
  where
    fn b = fromIntegral <$> {#get bam1_t->core.n_cigar #} b >>= \case
        0 -> return Nothing
        n -> allocaArray n $ \num -> allocaArray n $ \str -> do
            bamGetCigar b num str n
            num' <- peekArray (fromIntegral n) num
            str' <- peekArray (fromIntegral n) str
            return $ Just $ CIGAR $ zip (map fromIntegral num') $
                map castCCharToChar str'
{-# INLINE cigar #-}

-- | Ref-ID of the next segment (the paired read).
mateRefId :: BAM -> Int
mateRefId = unsafePerformIO . flip withForeignPtr fn . unbam
  where
    fn = fmap fromIntegral . {#get bam1_t->core.mtid #}
{-# INLINE mateRefId #-}

mateRefName :: BAMHeader -> BAM -> Maybe B.ByteString
mateRefName header bam
    | chr < 0 = Nothing
    | otherwise = Just $ unsafePerformIO $
        withForeignPtr (unbamHeader header) $ \h ->
            bamChr h chr >>= B.packCString
  where
    chr = mateRefId bam
{-# INLINE mateRefName #-}

-- | 0-based
mateStartLoc :: BAM -> Int
mateStartLoc = unsafePerformIO . flip withForeignPtr fn . unbam
  where
    fn = fmap fromIntegral . {#get bam1_t->core.mpos #}
{-# INLINE mateStartLoc #-}

-- | Observed Template length.
-- If all segments are mapped to the same reference, the unsigned
-- observed template length equals the number of bases from the
-- leftmost mapped base to the rightmost mapped base. 
tLen :: BAM -> Int
tLen = unsafePerformIO . flip withForeignPtr fn . unbam
  where
    fn = fmap fromIntegral . {#get bam1_t->core.isize #}
{-# INLINE tLen #-}

auxData :: BAM -> [((Char, Char), AuxiliaryData)]
auxData bam = unsafePerformIO $ withForeignPtr (unbam bam) $ \b -> do
    l <- bamGetLAux b
    aux <- bamGetAux b 
    go aux l
  where
    go ptr i
        | i <= 0 = return []
        | otherwise = do
            name <- (,) <$> peekByteOff ptr 0 <*> peekByteOff ptr 1
            castCCharToChar <$> peekByteOff ptr 2 >>= \case
                'A' -> do
                    r <- AuxChar . castCCharToChar <$> peekByteOff ptr 3
                    rs <- go (plusPtr ptr 4) $ i - 4
                    return $ (name, r) : rs
                -- int8_t
                'c' -> do
                    r <- AuxInt . fromIntegral <$> (peekByteOff ptr 3 :: IO Int8)
                    rs <- go (plusPtr ptr 4) $ i - 4
                    return $ (name, r) : rs
                -- uint8_t
                'C' -> do
                    r <- AuxInt . fromIntegral <$> (peekByteOff ptr 3 :: IO Word8)
                    rs <- go (plusPtr ptr 4) $ i - 4
                    return $ (name, r) : rs
                -- int16_t
                's' -> do
                    r <- AuxInt . fromIntegral <$> (peekByteOff ptr 3 :: IO Int16)
                    rs <- go (plusPtr ptr 5) $ i - 5
                    return $ (name, r) : rs
                -- uint16_t
                'S' -> do
                    r <- AuxInt . fromIntegral <$> (peekByteOff ptr 3 :: IO Word16)
                    rs <- go (plusPtr ptr 5) $ i - 5
                    return $ (name, r) : rs
                -- int32_t
                'i' -> do
                    r <- AuxInt . fromIntegral <$> (peekByteOff ptr 3 :: IO Int32)
                    rs <- go (plusPtr ptr 7) $ i - 7
                    return $ (name, r) : rs
                -- uint32_t
                'I' -> do
                    r <- AuxInt . fromIntegral <$> (peekByteOff ptr 3 :: IO Word32)
                    rs <- go (plusPtr ptr 7) $ i - 7
                    return $ (name, r) : rs
                'f' -> do
                    r <- AuxFloat <$> peekByteOff ptr 3
                    rs <- go (plusPtr ptr 7) $ i - 7
                    return $ (name, r) : rs
                'Z' -> do
                    str <- B.packCString (plusPtr ptr 3)
                    let l = B.length str + 1 + 3
                    rs <- go (plusPtr ptr l) $ i - l
                    return $ (name, AuxString str) : rs
                'H' -> do 
                    str <- B.packCString (plusPtr ptr 3)
                    let l = B.length str + 1 + 3
                    rs <- go (plusPtr ptr l) $ i - l
                    return $ (name, AuxByteArray str) : rs
                'B' -> do
                    n <- fromIntegral <$> (peekByteOff ptr 4 :: IO Int32)
                    castCCharToChar <$> peekByteOff ptr 3 >>= \case
                        'c' -> do
                            r <- AuxIntArray . map fromIntegral <$>
                                (peekArray n $ plusPtr ptr 8 :: IO [Int8])
                            let l = 8 + n
                            rs <- go (plusPtr ptr l) $ i - l
                            return $ (name, r) : rs
                        'C' -> do
                            r <- AuxIntArray . map fromIntegral <$>
                                (peekArray n $ plusPtr ptr 8 :: IO [Word8])
                            let l = 8 + n
                            rs <- go (plusPtr ptr l) $ i - l
                            return $ (name, r) : rs
                        's' -> do
                            r <- AuxIntArray . map fromIntegral <$>
                                (peekArray n $ plusPtr ptr 8 :: IO [Int16])
                            let l = 8 + 2 * n
                            rs <- go (plusPtr ptr l) $ i - l
                            return $ (name, r) : rs
                        'S' -> do
                            r <- AuxIntArray . map fromIntegral <$>
                                (peekArray n $ plusPtr ptr 8 :: IO [Word16])
                            let l = 8 + 2 * n
                            rs <- go (plusPtr ptr l) $ i - l
                            return $ (name, r) : rs
                        'i' -> do
                            r <- AuxIntArray . map fromIntegral <$>
                                (peekArray n $ plusPtr ptr 8 :: IO [Int32])
                            let l = 8 + 4 * n
                            rs <- go (plusPtr ptr l) $ i - l
                            return $ (name, r) : rs
                        'I' -> do
                            r <- AuxIntArray . map fromIntegral <$>
                                (peekArray n $ plusPtr ptr 8 :: IO [Word32])
                            let l = 8 + 4 * n
                            rs <- go (plusPtr ptr l) $ i - l
                            return $ (name, r) : rs
                        'f' -> do
                            r <- AuxFloatArray <$>
                                (peekArray n $ plusPtr ptr 8 :: IO [Float])
                            let l = 8 + 4 * n
                            rs <- go (plusPtr ptr l) $ i - l
                            return $ (name, r) : rs
                        x -> error $ "Unknown auxiliary record type: " ++ [x]
                x -> error $ "Unknown auxiliary record type: " ++ [x]
{-# INLINE auxData #-}

queryAuxData :: (Char, Char) -> BAM -> Maybe AuxiliaryData
queryAuxData (x1,x2) bam = unsafePerformIO $
    withForeignPtr (unbam bam) $ \b -> do
        ptr <- bamAuxGet b [x1,x2]
        if ptr == nullPtr
            then return Nothing
            else Just <$> getAuxData1 ptr
{-# INLINE queryAuxData #-}

getAuxData1 :: Ptr () -> IO AuxiliaryData
getAuxData1 ptr = castCCharToChar <$> peekByteOff ptr 0 >>= \case
    'A' -> AuxChar . castCCharToChar <$> peekByteOff ptr 1
    'c' -> AuxInt . fromIntegral <$> (peekByteOff ptr 1 :: IO Int8)
    'C' -> AuxInt . fromIntegral <$> (peekByteOff ptr 1 :: IO Word8)
    's' -> AuxInt . fromIntegral <$> (peekByteOff ptr 1 :: IO Int16)
    'S' -> AuxInt . fromIntegral <$> (peekByteOff ptr 1 :: IO Word16)
    'i' -> AuxInt . fromIntegral <$> (peekByteOff ptr 1 :: IO Int32)
    'I' -> AuxInt . fromIntegral <$> (peekByteOff ptr 1 :: IO Word32)
    'f' -> AuxFloat <$> peekByteOff ptr 1
    'Z' -> AuxString <$> B.packCString (plusPtr ptr 1)
    'H' -> AuxByteArray <$> B.packCString (plusPtr ptr 1)
    'B' -> do
        n <- fromIntegral <$> (peekByteOff ptr 2 :: IO Int32)
        peekByteOff ptr 1 >>= \case
            'c' -> AuxIntArray . map fromIntegral <$>
                (peekArray n $ plusPtr ptr 6 :: IO [Int8])
            'C' -> AuxIntArray . map fromIntegral <$>
                (peekArray n $ plusPtr ptr 6 :: IO [Word8])
            's' -> AuxIntArray . map fromIntegral <$>
                (peekArray n $ plusPtr ptr 6 :: IO [Int16])
            'S' -> AuxIntArray . map fromIntegral <$>
                (peekArray n $ plusPtr ptr 6 :: IO [Word16])
            'i' -> AuxIntArray . map fromIntegral <$>
                (peekArray n $ plusPtr ptr 6 :: IO [Int32])
            'I' -> AuxIntArray . map fromIntegral <$>
                (peekArray n $ plusPtr ptr 6 :: IO [Word32])
            'f' -> AuxFloatArray <$>
                (peekArray n $ plusPtr ptr 6 :: IO [Float])
            x -> error $ "Unknown auxiliary record type: " ++ [x]
    x -> error $ "Unknown auxiliary record type: " ++ [x]
{-# INLINE getAuxData1 #-}
            
-- | Convert Bam record to Sam record.
bamToSam :: BAMHeader -> BAM -> SAM
bamToSam h b = SAM (queryName b) (flag b) (refName h b) (startLoc b + 1) (mapq b)
    (cigar b) (mateRefName h b) (mateStartLoc b + 1) (tLen b) (getSeq b) (quality b)
    (auxData b)

-- | Template having multiple segments in sequencing
hasMultiSegments :: Word16 -> Bool
hasMultiSegments f = testBit f 0

-- | Each segment properly aligned according to the aligner
isProperAligned :: Word16 -> Bool
isProperAligned f = testBit f 1

-- | Segment unmapped
isUnmapped :: Word16 -> Bool
isUnmapped f = testBit f 2

-- | Next segment in the template unmapped
isNextUnmapped :: Word16 -> Bool
isNextUnmapped f = testBit f 3

-- | SEQ being reverse complemented
isRC :: Word16 -> Bool
isRC f = testBit f 4

-- | SEQ of the next segment in the template being reverse complemented
isNextRC :: Word16 -> Bool
isNextRC f = testBit f 5

-- | The first segment in the template
isFirstSegment :: Word16 -> Bool
isFirstSegment f = testBit f 6

-- | The last segment in the template
isLastSegment :: Word16 -> Bool
isLastSegment f = testBit f 7

-- | Secondary alignment
isSecondary :: Word16 -> Bool
isSecondary f = testBit f 8

-- | Not passing filters, such as platform/vendor quality controls
isBadQual :: Word16 -> Bool
isBadQual f = testBit f 9

-- | PCR or optical duplicate
isDup :: Word16 -> Bool
isDup f = testBit f 10

-- | Supplementary alignment
isSupplementary :: Word16 -> Bool
isSupplementary f = testBit f 11


--------------------------------------------------------------------------------
-- Modify BAM 
--------------------------------------------------------------------------------

-- | Append tag data to a bam record.
appendAux :: (Char, Char)  -- ^ Tag
          -> AuxiliaryData -- ^ Data
          -> BAM
          -> IO ()
appendAux (x1,x2) aux bam = withForeignPtr (unbam bam) $ \b -> do
    code <- case aux of
        AuxChar d -> with d $ bamAuxAppend b [x1,x2] 'A' 1
        AuxInt d -> with d $ bamAuxAppend b [x1,x2] 'i' 4
        AuxFloat d -> with d $ bamAuxAppend b [x1,x2] 'f' 4
        AuxString d -> B.useAsCString d $
            bamAuxAppend b [x1,x2] 'Z' (B.length d + 1) . castPtr
        _ -> error "Not implemented"
    if code == 0 then return () else error "Append aux failed"
  where
    with x fun = alloca $ \ptr -> poke ptr x >> fun (castPtr ptr)
{-# INLINE appendAux #-}

-- | Turn on duplicate flag.
setDup :: BAM -> IO ()
setDup bam = withForeignPtr (unbam bam) bamMarkDup
{-# INLINE setDup #-}