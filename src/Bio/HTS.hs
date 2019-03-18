{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE TemplateHaskell   #-}

module Bio.HTS where

import           Conduit
import           Control.Monad
import           Control.Monad.Reader
import           Data.Bits                (testBit)
import qualified Data.ByteString.Char8    as B
import qualified Data.ByteString as BS
import           Data.Int
import           Data.Monoid              ((<>))
import           Data.Word
import           Foreign.C.String
import           Foreign.C.Types
import Control.Exception (bracket)
import           Foreign.ForeignPtr
import Foreign.Storable (peek)
import           Foreign.Marshal.Array
import           Foreign.Ptr
import           Language.C.Inline        (baseCtx, bsCtx, context, include,
                                           withPtr)
import qualified Language.C.Inline.Unsafe as CU
import           System.IO
import           System.IO.Unsafe         (unsafePerformIO)

import           Bio.HTS.Types

context (baseCtx <> bsCtx <> htsCtx)

include "htslib/sam.h"

type HeaderState = ResourceT (ReaderT FileHeader IO)

withBamFile :: FilePath -> (BamFileHandle -> HeaderState a) -> IO a
withBamFile fl action = bracket (openBamFile fl ReadMode) closeBamFile $ \h -> do
    hdr <- readBamHeader h
    runReaderT (runResourceT $ action h) hdr

readBam :: BamFileHandle -> ConduitT () Bam HeaderState ()
readBam h = do
    r <- liftIO $ bamRead1 h
    case r of
        Nothing -> return ()
        Just b -> yield b >> readBam h
{-# INLINE readBam #-}

writeBam :: FilePath -> ConduitT Bam o HeaderState ()
writeBam fn = bracketP (openBamFile fn WriteMode) closeBamFile $ \(BamFileHandle fp) -> do
    maybeBam <- await
    case maybeBam of
        Nothing -> return ()
        Just b' -> do
            leftover b'
            header <- lift ask
            case header of
                BamHeader hdr -> do
                    err <- liftIO $ [CU.exp| int {
                        bam_hdr_write($(htsFile* fp)->fp.bgzf, $(bam_hdr_t* hdr)) } |]
                    if err /= 0
                        then error "'bam_hdr_write' failed."
                        else sink fp
                _ -> error "No Bam header was found!"
  where
    sink fp = awaitForever $ \b' -> liftIO $ withForeignPtr b' $ \b ->
        [CU.exp| int { bam_write1($(htsFile* fp)->fp.bgzf, $(bam1_t* b)) } |]
{-# INLINE writeBam #-}

openBamFile :: FilePath -> IOMode -> IO BamFileHandle
openBamFile fn mode = do
    fn' <- newCString fn
    h <- case mode of
        ReadMode -> [CU.exp| htsFile* { hts_open($(char* fn'), "r") } |]
        WriteMode -> [CU.exp| htsFile* { hts_open($(char* fn'), "wb") } |]
        _ -> error ""
    return $ BamFileHandle h

readBamHeader :: BamFileHandle -> IO FileHeader
readBamHeader (BamFileHandle h) =
    BamHeader <$> [CU.exp| bam_hdr_t* { bam_hdr_read($(htsFile* h)->fp.bgzf) } |]

closeBamFile :: BamFileHandle -> IO ()
closeBamFile (BamFileHandle h) = [CU.exp| void { hts_close($(htsFile* h)) } |]


showBamHeader :: FileHeader -> B.ByteString
showBamHeader (BamHeader hdr) = unsafePerformIO $ do
    ptr <- [CU.exp| char* { $(bam_hdr_t* hdr)->text } |]
    l <- [CU.exp| uint32_t { $(bam_hdr_t* hdr)->l_text } |]
    B.packCStringLen (ptr, fromIntegral l)
showBamHeader _ = error "No Bam Header was found."

data SortOrder = Unknown
               | Unsorted
               | Queryname
               | Coordinate

-- | Get the sort information.
getSortOrder :: FileHeader -> SortOrder
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

-- | Read one record.
bamRead1 :: BamFileHandle -> IO (Maybe Bam)
bamRead1 (BamFileHandle h) = do
    (r, b) <- withPtr $ \r -> [CU.block| bam1_t* {
            bam1_t *b = bam_init1();
            *$(int* r) = bam_read1($(htsFile* h)->fp.bgzf, b);
            return b;
        } |]
    if r < 0 then return Nothing else Just <$> newForeignPtr bamDestory b
{-# INLINE bamRead1 #-}

foreign import ccall unsafe "&bam_destroy1"
   bamDestory :: FunPtr (Ptr Bam' -> IO ())

-- | Return the chromosome id.
getChrId :: Bam -> Int32
getChrId = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| int32_t { $(bam1_t* b)->core.tid } |]
{-# INLINE getChrId #-}

-- | Return the chromosome name given the bam file header.
getChr :: Ptr BamHdr -> Bam -> Maybe B.ByteString
getChr h b' | i < 0 = Nothing
            | otherwise = Just $ unsafePerformIO $ join $ B.packCString <$>
                [CU.exp| char* { $(bam_hdr_t* h)->target_name[$(int32_t i)] } |]
  where
    i = getChrId b'
{-# INLINE getChr #-}

-- | Return the 0-based starting location.
position :: Bam -> Int32
position = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| int32_t { $(bam1_t* b)->core.pos } |]
{-# INLINE position #-}

-- | For a mapped read, this is just position + cigar2rlen.
-- For an unmapped read (either according to its flags or if it has no cigar
-- string), we return position + 1 by convention.
endPos :: Bam -> Int32
endPos = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| int32_t { bam_endpos($(bam1_t* b)) } |]
{-# INLINE endPos #-}

-- | Return the query length (read length).
queryLen :: Bam -> Int32
queryLen = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| int32_t { $(bam1_t* b)->core.l_qseq } |]
{-# INLINE queryLen #-}

-- | Whether the query is on the reverse strand.
isRev :: Bam -> Bool
isRev = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = do
        r <- [CU.exp| int {bam_is_rev($(bam1_t* b)) } |]
        return $ if r == 0 then False else True
{-# INLINE isRev #-}

-- | Return the flag.
flag :: Bam -> Word16
flag = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| uint16_t { $(bam1_t* b)->core.flag } |]
{-# INLINE flag #-}

-- | MAPping Quality. It equals −10 log10 Pr{mapping position is wrong}, rounded to the nearest integer. A value 255 indicates that the mapping quality is not available.
mapq :: Bam -> Word8
mapq = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| uint8_t { $(bam1_t* b)->core.qual } |]
{-# INLINE mapq #-}


--int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
-- | Return the DNA sequence.
getSeq :: Bam -> Maybe B.ByteString
getSeq = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = do
        l <- [CU.exp| int32_t { $(bam1_t* b)->core.l_qseq } |]
        if (l == 0)
            then return Nothing
            else allocaArray (fromIntegral l) $ \str -> do
                    [CU.block| void {
                        int32_t i;
                        uint8_t *s = bam_get_seq($(bam1_t* b));
                        for (i = 0; i < $(int32_t l); ++i)
                            $(char* str)[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(s, i)];
                    } |]
                    Just <$> B.packCStringLen (str, fromIntegral l)
{-# INLINE getSeq #-}


-- | Get the name of the query.
qName :: Bam -> B.ByteString
qName = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = join $ B.packCString <$>
        [CU.exp| char* {bam_get_qname($(bam1_t* b)) } |]
{-# INLINE qName #-}

-- | Human readable quality score which is: Phred base quality + 33.
qualityS :: Bam -> Maybe B.ByteString
qualityS = fmap (BS.map (+33)) . quality

-- | Phred base quality (a sequence of 0xFF if absent).
quality :: Bam -> Maybe B.ByteString
quality = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = do
        l <- [CU.exp| int32_t { $(bam1_t* b)->core.l_qseq } |]
        if l == 0
            then return Nothing
            else do
                x <- [CU.block| int8_t {
                        uint8_t *s = bam_get_qual($(bam1_t* b));
                        return (s[0] == 0xff);
                     } |]
                if x /= 0
                    then return Nothing
                    else allocaArray (fromIntegral l) $ \str -> do
                            [CU.block| void {
                                int32_t i;
                                uint8_t *s = bam_get_qual($(bam1_t* b));
                                for (i = 0; i < $(int32_t l); ++i)
                                    $(char* str)[i] = s[i];
                            } |]
                            Just <$> B.packCStringLen (str, fromIntegral l)
{-# INLINE quality #-}

cigar :: Bam -> Maybe [(Int, Char)]
cigar = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = do
        n <- [CU.exp| uint16_t { $(bam1_t* b)->core.n_cigar } |]
        if (n == 0)
            then return Nothing
            else allocaArray (fromIntegral n) $ \num ->
                     allocaArray (fromIntegral n) $ \str -> do
                         [CU.block| void {
                            uint16_t i;
                            uint32_t *cigar = bam_get_cigar($(bam1_t* b));
                            for (i = 0; i < $(uint16_t n); ++i) {
                                $(int* num)[i] = bam_cigar_oplen(cigar[i]);
                                $(char* str)[i] = bam_cigar_opchr(cigar[i]);
                            }
                         } |]
                         num' <- peekArray (fromIntegral n) num
                         str' <- peekArray (fromIntegral n) str
                         return $ Just $
                            zip (map fromIntegral num') $ map castCCharToChar str'
{-# INLINE cigar #-}

mateChrId :: Bam -> Int32
mateChrId = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| int32_t { $(bam1_t* b)->core.mtid } |]
{-# INLINE mateChrId #-}

mateChr :: Ptr BamHdr -> Bam -> Maybe B.ByteString
mateChr h b' | i < 0 = Nothing
             | otherwise = Just $ unsafePerformIO $ join $ B.packCString <$>
                [CU.exp| char* { $(bam_hdr_t* h)->target_name[$(int32_t i)] } |]
  where
    i = mateChrId b'
{-# INLINE mateChr #-}

-- | 0-based
matePos :: Bam -> Int32
matePos = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| int32_t { $(bam1_t* b)->core.mpos } |]
{-# INLINE matePos #-}

tLen :: Bam -> Int32
tLen = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| int32_t { $(bam1_t* b)->core.isize } |]
{-# INLINE tLen #-}

auxData :: Bam -> [(B.ByteString, AuxiliaryData)]
auxData bam = unsafePerformIO $ withForeignPtr bam $ \b -> do
    l <- [CU.exp| int32_t { bam_get_l_aux($(bam1_t* b)) } |]
    aux <- [CU.exp| uint8_t* { bam_get_aux($(bam1_t* b)) } |]
    go aux $ fromIntegral l
  where
    go ptr i
        | i <= 0 = return []
        | otherwise = do
            name <- B.packCStringLen (castPtr ptr, 2)
            castCCharToChar <$> (peek $ plusPtr ptr 2) >>= \case
                'A' -> do
                    r <- AuxChar . castCCharToChar <$> peek (plusPtr ptr 3)
                    rs <- go (plusPtr ptr 4) $ i - 4
                    return $ (name, r) : rs
                -- int8_t
                'c' -> do
                    r <- AuxInt . fromIntegral <$> (peek $ plusPtr ptr 3 :: IO Int8)
                    rs <- go (plusPtr ptr 4) $ i - 4
                    return $ (name, r) : rs
                -- uint8_t
                'C' -> do
                    r <- AuxInt . fromIntegral <$> (peek $ plusPtr ptr 3 :: IO Word8)
                    rs <- go (plusPtr ptr 4) $ i - 4
                    return $ (name, r) : rs
                -- int16_t
                's' -> do
                    r <- AuxInt . fromIntegral <$> (peek $ plusPtr ptr 3 :: IO Int16)
                    rs <- go (plusPtr ptr 5) $ i - 5
                    return $ (name, r) : rs
                -- uint16_t
                'S' -> do
                    r <- AuxInt . fromIntegral <$> (peek $ plusPtr ptr 3 :: IO Word16)
                    rs <- go (plusPtr ptr 5) $ i - 5
                    return $ (name, r) : rs
                -- int32_t
                'i' -> do
                    r <- AuxInt . fromIntegral <$> (peek $ plusPtr ptr 3 :: IO Int32)
                    rs <- go (plusPtr ptr 7) $ i - 7
                    return $ (name, r) : rs
                -- uint32_t
                'I' -> do
                    r <- AuxInt . fromIntegral <$> (peek $ plusPtr ptr 3 :: IO Word32)
                    rs <- go (plusPtr ptr 7) $ i - 7
                    return $ (name, r) : rs
                'f' -> do
                    r <- AuxFloat <$> peek (plusPtr ptr 3)
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
                    n <- fromIntegral <$> (peek $ plusPtr ptr 4 :: IO Int32)
                    castCCharToChar <$> (peek $ plusPtr ptr 3) >>= \case
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
            
-- | Convert Bam record to Sam record.
bamToSam :: Ptr BamHdr -> Bam -> Sam
bamToSam h b = Sam (qName b) (flag b) (getChr h b) (position b) (mapq b)
    (cigar b) (mateChr h b) (matePos b) (tLen b) (getSeq b) (quality b)
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
