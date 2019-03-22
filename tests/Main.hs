{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleContexts #-}
import Test.Tasty
import Test.Tasty.Golden
import Conduit

import Bio.HTS

main :: IO ()
main = defaultMain $ testGroup "Main"
    [ bamToSamTest
    , markDupTest ]

bamToSamTest :: TestTree
bamToSamTest = goldenVsFile "BAM Read/Write Test" expect output $
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

markDupTest :: TestTree
markDupTest = testGroup "Mark duplicates"
    [ markDup "Single end" "tests/data/single_end_dedup.bam" 
        "tests/data/single_end.bam" "dedup.bam" 
    , markDup "Paired end" "tests/data/paired_end_dedup.bam" 
        "tests/data/paired_end.bam" "dedup.bam" 
    ]
  where
    markDup name expect input output = goldenVsFile name expect output $ do
        header <- getBamHeader input
        runResourceT $ runConduit $ streamBam input .|
            markDupBy (const Nothing) .|
            filterC (not . isDup . flag) .| sinkBam output header
