{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE TemplateHaskell   #-}

module Bio.Data.HTS where

import           Control.Monad
import qualified Data.ByteString.Char8 as B
import           Data.Monoid           ((<>))
import           Data.Word
import           Foreign.C.String
import           Foreign.C.Types
import           Foreign.Ptr
import Language.C.Inline (context, baseCtx, bsCtx, include, withPtr, withPtr_)
import qualified Language.C.Inline.Unsafe as CU
import           System.IO.Unsafe      (unsafePerformIO)
import Foreign.Storable (peek)
import Foreign.Marshal.Array
import Data.Int

import           Bio.Data.Types

context (baseCtx <> bsCtx <> htsCtx)

include "htslib/sam.h"

openBamFile :: FilePath -> IO BamHandle
openBamFile fn = do
    fn' <- newCString fn
    [CU.exp| samFile* { sam_open($(char* fn'), "r") } |]

samHdrRead :: BamHandle -> IO BamHeader
samHdrRead fp = [CU.exp| bam_hdr_t* { sam_hdr_read($(samFile* fp)) } |]

samRead1 :: BamHandle -> BamHeader -> IO Bam
samRead1 fp hdr =
    [CU.block| bam1_t* {
        int err;
        bam1_t *b;
        b = bam_init1();
        err = sam_read1($(samFile* fp), $(bam_hdr_t* hdr), b);
        return b;
    } |]


getChrId :: Bam -> Int32
getChrId b = [CU.pure| int32_t { $(bam1_t* b)->core.tid } |]

getChr :: BamHeader -> Bam -> Maybe B.ByteString
getChr h b = if i < 0 then Nothing else Just $ unsafePerformIO $ do
    join $ B.packCString <$>
        [CU.exp| char* { $(bam_hdr_t* h)->target_name[$(int32_t i)] } |]
  where
    i = getChrId b

position :: Bam -> Int32
position b = [CU.pure| int32_t { $(bam1_t* b)->core.pos } |]

-- | For a mapped read, this is just position + cigar2rlen.
-- For an unmapped read (either according to its flags or if it has no cigar
-- string), we return position + 1 by convention.
endPos :: Bam -> Int32
endPos b = [CU.pure| int32_t { bam_endpos($(bam1_t* b)) } |]

queryLen :: Bam -> Int32
queryLen b = [CU.pure| int32_t { $(bam1_t* b)->core.l_qseq } |]

-- | Whether the query is on the reverse strand.
isRev :: Bam -> Bool
isRev b = if [CU.pure| int {bam_is_rev($(bam1_t* b)) } |] == 0 then False else True

flag :: Bam -> Word16
flag b = [CU.pure| uint16_t { $(bam1_t* b)->core.flag } |]

mapq :: Bam -> Word8
mapq b = [CU.pure| uint8_t { $(bam1_t* b)->core.qual } |]


--int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
getSeq :: Bam -> Maybe B.ByteString
getSeq b = unsafePerformIO $ do
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
qName b = unsafePerformIO $ join $ B.packCString <$>
    [CU.exp| char* {bam_get_qname($(bam1_t* b)) } |]

quality :: Bam -> Maybe B.ByteString
quality b = unsafePerformIO $ do
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
cigar b = unsafePerformIO $ do
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
mateChrId b = [CU.pure| int32_t { $(bam1_t* b)->core.mtid } |]

mateChr :: BamHeader -> Bam -> Maybe B.ByteString
mateChr h b = if i < 0 then Nothing else Just $ unsafePerformIO $ do
    join $ B.packCString <$>
        [CU.exp| char* { $(bam_hdr_t* h)->target_name[$(int32_t i)] } |]
  where
    i = mateChrId b


-- | 0-based
matePos :: Bam -> Int32
matePos b = [CU.pure| int32_t { $(bam1_t* b)->core.mpos } |]

tLen :: Bam -> Int32
tLen b = [CU.pure| int32_t { $(bam1_t* b)->core.isize } |]


bamToSam :: BamHeader -> Bam -> Sam
bamToSam h b = Sam (qName b) (flag b) (getChr h b) (position b) (mapq b)
    (cigar b) (mateChr h b) (matePos b) (tLen b) (getSeq b) (quality b)
