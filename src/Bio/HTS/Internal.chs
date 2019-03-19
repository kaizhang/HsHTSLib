{-# LANGUAGE ForeignFunctionInterface #-}
module Bio.HTS.Internal
    ( HTSFile
    , withHTSFile
    , htsOpen
    , htsClose

    -- * BGZF
    , BGZF
    , getBgzf

    -- * Bam Header
    , BamHdr
    , bamHdrRead
    , bamHdrWrite
    , getHeaderText
    , getHeaderSize

    -- * Bam
    , Bam1
    , bamRead1
    , bamWrite1
    , bamChr
    , bamEndpos
    , bamIsRev
    , bamGetSeq
    , bamGetQual
    , bamGetCigar
    , bamGetAux
    , bamGetLAux
    , bamAuxAppend
    ) where

import Control.Exception (bracket)
import Foreign
import Foreign.C.Types
import Foreign.C.String
import System.IO (IOMode(..))
        
#include "cbits/aux.c"
#include "htslib/sam.h"

-- | Opaque data representing the hts file.
data HTSFile

withHTSFile :: FilePath -> IOMode -> (Ptr HTSFile -> IO a) -> IO a
withHTSFile fl mode = bracket (htsOpen fl md) htsClose
  where
    md = case mode of
        ReadMode -> "r"
        WriteMode -> "wb"
        _ -> error "Not a supported mode"

-- | Open a hts file.
{#fun hts_open as ^
    { withCString* `String'
    , withCString* `String'
    } -> `Ptr HTSFile' castPtr #}

-- | Close the hts file.
{#fun hts_close as ^ { castPtr `Ptr HTSFile'} -> `CInt' void- #}


--------------------------------------------------------------------------------
-- BGZF
--------------------------------------------------------------------------------

-- | Opaque data representing the Blocked GNU Zip Format (BGZF) used
-- by Bam.
data BGZF

-- | Get the location of BGZF block.
{#fun get_bgzf as ^ { castPtr `Ptr HTSFile' } -> `Ptr BGZF' castPtr #}


--------------------------------------------------------------------------------
-- BAM
--------------------------------------------------------------------------------

-- | Opaque data representing Bam type.
data Bam1

-- | Read one Bam record.
{#fun bam_read1 as ^ { castPtr `Ptr BGZF', castPtr `Ptr Bam1'} -> `CInt' #}

-- | Save one Bam record.
{#fun bam_write1 as ^ { castPtr `Ptr BGZF', castPtr `Ptr Bam1'} -> `CInt' #}

-- | Get chromosome.
{#fun bam_chr as ^ {castPtr `Ptr BamHdr', `Int'} -> `CString' #}

-- | Get end position.
{#fun bam_endpos as ^ {castPtr `Ptr Bam1'} -> `Int' #}

-- | Is reverse.
{#fun bam_is_rev_ as bamIsRev {castPtr `Ptr Bam1'} -> `Bool' #}

{#fun bam_get_seq_ as bamGetSeq {castPtr `Ptr Bam1', `CString', `Int' } -> `()' #}

{#fun bam_get_qual_ as bamGetQual {castPtr `Ptr Bam1', `CString', `Int' } -> `CInt' #}

{#fun bam_get_cigar_ as bamGetCigar
    { castPtr `Ptr Bam1', castPtr `Ptr CInt', `CString', `Int' } -> `()' #}

{#fun bam_get_aux_ as bamGetAux { castPtr `Ptr Bam1' } -> `Ptr ()' castPtr #}

{#fun bam_get_l_aux_ as bamGetLAux { castPtr `Ptr Bam1' } -> `Int' #}

{#fun bam_aux_append as ^
    { castPtr `Ptr Bam1', `String', `Char', `Int', castPtr `Ptr ()'
    } -> `CInt' #}

--------------------------------------------------------------------------------
-- BAM Header
--------------------------------------------------------------------------------

-- | Opaque data representing Bam header type.
data BamHdr

-- | Get the bam header.
{#fun bam_hdr_read as ^ { castPtr `Ptr BGZF' } -> `Ptr BamHdr' castPtr #}

-- | Save the bam header.
{#fun bam_hdr_write as ^ { castPtr `Ptr BGZF',  castPtr `Ptr BamHdr' } -> `CInt' #}

{#fun get_header_text as ^ { castPtr `Ptr BamHdr' } -> `CString' #}
{#fun get_header_size as ^ { castPtr `Ptr BamHdr' } -> `Int' #}