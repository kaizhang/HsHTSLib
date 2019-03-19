{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleContexts #-}
import Test.Tasty
import           Test.Tasty.Golden
import Conduit
import Control.Monad.Reader

import Bio.HTS
import Bio.HTS.Types

main :: IO ()
main = defaultMain $ testGroup "Main"
    [ tests ]

tests :: TestTree
tests = goldenVsFile "BAM Read/Write Test" expect output $
    bamFileToSamFile input output
  where
    input = "tests/data/example.bam"
    output = "tests/data/output.sam"
    expect = "tests/data/example.sam"

bamFileToSamFile :: FilePath  -- ^ Input bam file
                 -> FilePath  -- ^ Output Sam file
                 -> IO ()
bamFileToSamFile input output = do
    header <- getBamHeader input
    runResourceT $ runConduit $ streamBam input .|
        mapC (showSam . bamToSam header) .| unlinesAsciiC .| sinkFile output