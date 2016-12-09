{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE TemplateHaskell   #-}

module Bio.HTS.Types
    ( htsCtx
    , BamFileHandle(..)
    , FileHeader(..)
    , BamHdr
    , HTSFile
    , Bam
    , Bam'
    , Sam(..)
    , Flag(..)
    , showSam
    ) where

import qualified Data.ByteString              as BS
import qualified Data.ByteString.Char8        as B
import           Data.ByteString.Lex.Integral
import           Data.Int
import qualified Data.Map                     as M
import           Data.Maybe                   (fromJust, fromMaybe)
import           Data.Word
import           Foreign.ForeignPtr
import           Foreign.Ptr
import qualified Language.C.Inline            as C
import qualified Language.C.Inline.Context    as C
import qualified Language.C.Types             as C
import qualified Language.Haskell.TH          as TH

type Bam = ForeignPtr Bam'
data Bam'

data FileHeader = Empty
                | BamHeader (Ptr BamHdr)
data BamHdr

data HTSFile
newtype BamFileHandle = BamFileHandle (Ptr HTSFile)

-- | SAM record flag
newtype Flag = Flag Word16

htsCtx :: C.Context
htsCtx = mempty
  { C.ctxTypesTable = htsTypesTable
  }

htsTypesTable :: M.Map C.TypeSpecifier TH.TypeQ
htsTypesTable = M.fromList
   [ (C.TypeName "htsFile", [t| HTSFile |])
   , (C.TypeName "bam_hdr_t", [t| BamHdr |])
   , (C.TypeName "bam1_t", [t| Bam' |])
   ]

-- | currently only contains 11 mandatory fields
data Sam = Sam
    { samQname :: !B.ByteString
    , samFlag  :: !Word16
    , samRname :: !(Maybe B.ByteString)
    , samPos   :: !Int32
    , samMapq  :: !Word8
    , samCigar :: !(Maybe [(Int, Char)])
    , samRnext :: !(Maybe B.ByteString)
    , samPnext :: !Int32
    , samTlen  :: !Int32
    , samSeq   :: !(Maybe B.ByteString)
    , samQual  :: !(Maybe B.ByteString)
    } deriving (Show)

showSam :: Sam -> B.ByteString
showSam s = B.intercalate "\t"
    [ samQname s, pack' $ samFlag s, fromMaybe "*" $ samRname s, pack' $ samPos s + 1
    , pack' $ samMapq s, fromMaybe "*" $ f <$> samCigar s, fromMaybe "*" $ samRnext s
    , pack' $ samPnext s + 1, pack' $ samTlen s, fromMaybe "*" $ samSeq s
    , fromMaybe "*" $ BS.map (+33) <$> samQual s ]
  where
    f = B.concat . concatMap (\(i, x) -> [pack' i, B.singleton x])
    pack' :: Integral a => a -> B.ByteString
    pack' = fromJust . packDecimal
