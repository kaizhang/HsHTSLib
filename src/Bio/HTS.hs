{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE TemplateHaskell   #-}

module Bio.HTS where
-- FIXME: fix memory leaks using withPtr

import           Conduit
import           Control.Monad
import           Control.Monad.State
import qualified Data.ByteString.Char8    as B
import           Data.Int
import           Data.Monoid              ((<>))
import           Data.Word
import           Foreign.C.String
import           Foreign.C.Types
import           Foreign.ForeignPtr
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

readBam :: FilePath -> Source (ResourceT (StateT FileHeader IO)) Bam
readBam fn = bracketP (openBamFile fn ReadMode) closeBamFile $ \h -> do
    hdr <- liftIO (readBamHeader h)
    lift $ put hdr
    source h
  where
    source x = do
        r <- liftIO $ bamRead1 x
        case r of
            Nothing -> return ()
            Just b -> yield b >> source x

writeBam :: FilePath -> Sink Bam (ResourceT (StateT FileHeader IO)) ()
writeBam fn = bracketP (openBamFile fn WriteMode) closeBamFile $ \(BamFileHandle fp) -> do
    maybeBam <- await
    case maybeBam of
        Nothing -> return ()
        Just b' -> do
            leftover b'
            header <- lift get
            case header of
                BamHeader hdr -> do
                    err <- liftIO $ [CU.exp| int {
                        bam_hdr_write($(htsFile* fp)->fp.bgzf, $(bam_hdr_t* hdr)) } |]
                    if err /= 0
                        then error "'bam_hdr_write' failed."
                        else sink fp
                _ -> error "No header was provided."
  where
    sink fp = awaitForever $ \b' -> liftIO $ withForeignPtr b' $ \b ->
        [CU.exp| int { bam_write1($(htsFile* fp)->fp.bgzf, $(bam1_t* b)) } |]

openBamFile :: FilePath -> IOMode -> IO BamFileHandle
openBamFile fn mode = do
    fn' <- newCString fn
    h <- case mode of
        ReadMode -> [CU.exp| htsFile* { hts_open($(char* fn'), "r") } |]
        WriteMode -> [CU.exp| htsFile* { hts_open($(char* fn'), "wb") } |]
        _ -> error ""
    return $ BamFileHandle h

closeBamFile :: BamFileHandle -> IO ()
closeBamFile (BamFileHandle h) = [CU.exp| void { hts_close($(htsFile* h)) } |]

readBamHeader :: BamFileHandle -> IO FileHeader
readBamHeader (BamFileHandle h) =
    BamHeader <$> [CU.exp| bam_hdr_t* { bam_hdr_read($(htsFile* h)->fp.bgzf) } |]

bamRead1 :: BamFileHandle -> IO (Maybe Bam)
bamRead1 (BamFileHandle h) = do
    (r, b) <- withPtr $ \r -> [CU.block| bam1_t* {
            bam1_t *b = bam_init1();
            *$(int* r) = bam_read1($(htsFile* h)->fp.bgzf, b);
            return b;
        } |]
    if r < 0 then return Nothing else Just <$> newForeignPtr bamDestory b

foreign import ccall unsafe "&bam_destroy1"
   bamDestory :: FunPtr (Ptr Bam1' -> IO ())

-- | Return the chromosome id.
getChrId :: Bam -> Int32
getChrId = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| int32_t { $(bam1_t* b)->core.tid } |]

-- | Return the chromosome name given the bam file header.
getChr :: Ptr BamHdr -> Bam -> Maybe B.ByteString
getChr h b' | i < 0 = Nothing
            | otherwise = Just $ unsafePerformIO $ join $ B.packCString <$>
                [CU.exp| char* { $(bam_hdr_t* h)->target_name[$(int32_t i)] } |]
  where
    i = getChrId b'

-- | Return the 0-based starting location.
position :: Bam -> Int32
position = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| int32_t { $(bam1_t* b)->core.pos } |]

-- | For a mapped read, this is just position + cigar2rlen.
-- For an unmapped read (either according to its flags or if it has no cigar
-- string), we return position + 1 by convention.
endPos :: Bam -> Int32
endPos = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| int32_t { bam_endpos($(bam1_t* b)) } |]

-- | Return the query length (read length).
queryLen :: Bam -> Int32
queryLen = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| int32_t { $(bam1_t* b)->core.l_qseq } |]

-- | Whether the query is on the reverse strand.
isRev :: Bam -> Bool
isRev = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = do
        r <- [CU.exp| int {bam_is_rev($(bam1_t* b)) } |]
        return $ if r == 0 then False else True

-- | Return the flag.
flag :: Bam -> Word16
flag = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| uint16_t { $(bam1_t* b)->core.flag } |]

mapq :: Bam -> Word8
mapq = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| uint8_t { $(bam1_t* b)->core.qual } |]


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


-- | Get the name of the query.
qName :: Bam -> B.ByteString
qName = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = join $ B.packCString <$>
        [CU.exp| char* {bam_get_qname($(bam1_t* b)) } |]

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

mateChrId :: Bam -> Int32
mateChrId = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| int32_t { $(bam1_t* b)->core.mtid } |]

mateChr :: Ptr BamHdr -> Bam -> Maybe B.ByteString
mateChr h b' | i < 0 = Nothing
             | otherwise = Just $ unsafePerformIO $ join $ B.packCString <$>
                [CU.exp| char* { $(bam_hdr_t* h)->target_name[$(int32_t i)] } |]
  where
    i = mateChrId b'

-- | 0-based
matePos :: Bam -> Int32
matePos = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| int32_t { $(bam1_t* b)->core.mpos } |]

tLen :: Bam -> Int32
tLen = unsafePerformIO . flip withForeignPtr fn
  where
    fn b = [CU.exp| int32_t { $(bam1_t* b)->core.isize } |]


-- | Convert Bam record to Sam record.
bamToSam :: Ptr BamHdr -> Bam -> Sam
bamToSam h b = Sam (qName b) (flag b) (getChr h b) (position b) (mapq b)
    (cigar b) (mateChr h b) (matePos b) (tLen b) (getSeq b) (quality b)
