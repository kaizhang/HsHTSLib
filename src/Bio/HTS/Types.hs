{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}

module Bio.HTS.Types
    ( BAM(..)
    , BAMHeader(..)
    , showBamHeader
    , SortOrder(..)
    , SAM(..)
    , CIGAR(..)
    , cigar2String
    , string2Cigar
    , AuxiliaryData(..)
    , showSam
    ) where

import qualified Data.ByteString              as BS
import qualified Data.ByteString.Char8        as B
import           Data.ByteString.Lex.Integral
import           Data.Maybe                   (fromJust, fromMaybe)
import           System.IO.Unsafe         (unsafePerformIO)
import           Data.Word
import           Foreign.ForeignPtr

import Bio.HTS.Internal

-- | The BAM format.
newtype BAM = BAM { unbam :: ForeignPtr Bam1 }

-- | The BAM file header.
newtype BAMHeader = BAMHeader {unbamHeader :: ForeignPtr BamHdr}

-- | Convert bam file header to string.
showBamHeader :: BAMHeader -> B.ByteString
showBamHeader header = unsafePerformIO $
    withForeignPtr (unbamHeader header) $ \ptr -> do
        s <- getHeaderText ptr
        l <- getHeaderSize ptr
        B.packCStringLen (s, l)

-- | BAM sort order.
data SortOrder = Unknown
               | Unsorted
               | Queryname
               | Coordinate
               deriving (Show, Eq)

newtype CIGAR = CIGAR [(Int, Char)]

-- | Convert CIGAR to string.
cigar2String :: CIGAR -> B.ByteString
cigar2String (CIGAR c) = B.concat $
    concatMap (\(i, x) -> [fromJust $ packDecimal i, B.singleton x]) c

-- | Read CIGAR from string.
string2Cigar :: B.ByteString -> CIGAR
string2Cigar c = CIGAR $ go c
  where
    go x = case readDecimal x of
        Nothing -> error $ "parse cigar fail: " ++ show c
        Just (n, remain) -> if B.length remain == 1
            then [(n, B.head remain)]
            else (n, B.head remain) : go (B.tail remain)
        
-- | The SAM format. The SAM format uses 1-based coordinate system.
data SAM = SAM
    { _sam_qname :: !B.ByteString
    , _sam_flag  :: !Word16
    , _sam_rname :: !(Maybe B.ByteString)
    , _sam_pos   :: !Int
    , _sam_mapq  :: !Word8
    , _sam_cigar :: !(Maybe CIGAR)
    , _sam_rnext :: !(Maybe B.ByteString)
    , _sam_pnext :: !Int
    , _sam_tlen  :: !Int
    , _sam_seq   :: !(Maybe B.ByteString)
    , _sam_qual  :: !(Maybe B.ByteString)
    , _sam_aux   :: [((Char,Char), AuxiliaryData)]
    }

-- | Convert SAM to string.
showSam :: SAM -> B.ByteString
showSam SAM{..} = B.intercalate "\t" $
    [ _sam_qname, pack' _sam_flag, fromMaybe "*" _sam_rname, pack' _sam_pos
    , pack' _sam_mapq, fromMaybe "*" $ cigar2String <$> _sam_cigar
    , fromMaybe "*" _sam_rnext, pack' _sam_pnext, pack' _sam_tlen
    , fromMaybe "*" _sam_seq, fromMaybe "*" $ BS.map (+33) <$> _sam_qual
    ] ++ map showAuxiliaryData _sam_aux
  where
    pack' :: Integral a => a -> B.ByteString
    pack' = fromJust . packDecimal

-- | Auxiliary data
data AuxiliaryData = AuxChar Char
                   | AuxInt Int
                   | AuxFloat Float
                   | AuxString B.ByteString
                   | AuxByteArray BS.ByteString
                   | AuxIntArray [Int]
                   | AuxFloatArray [Float]
                   deriving (Show)

-- | Convert aux data to string.
showAuxiliaryData :: ((Char, Char), AuxiliaryData) -> B.ByteString
showAuxiliaryData ((x1,x2), aux) = B.pack [x1,x2] <> aux'
  where
    aux' = case aux of
        AuxChar x -> B.pack [':', 'A', ':', x]
        AuxInt x -> B.pack $ ":i:" <> show x
        AuxFloat x -> B.pack $ ":f:" <> show x
        AuxString x -> ":Z:" <> x
        _ -> error "Not implemented"