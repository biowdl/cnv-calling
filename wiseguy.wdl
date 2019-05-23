version 1.0

# Copyright (c) 2019 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import "tasks/common.wdl" as common
import "tasks/wiseguy.wdl" as wiseguy
import "tasks/samtools.wdl" as samtools
import "tasks/bedtools.wdl" as bedtools

workflow wiseguyCnv {
    input {
    # IndexedBamFile from common.wdl
    Array[IndexedBamFile] referenceBams
    IndexedBamFile sample
    File reference
    File referenceIndex
    File? binFile
    String outputDir = "."
    }

    # Prepare reference
    scatter (refBam in referenceBams) {
        call wiseguy.Count as wiseguyCountReference {
            input:
                inputBam = refBam.file,
                inputBamIndex = refBam.index,
                reference = reference,
                referenceIndex = referenceIndex,
                binFile = binFile,
                outputBed = outputDir + "/references/" + basename(refBam.file) + ".bed"
        }

        call wiseguy.GcCorrect as wiseguyGcCorrectReference {
            input:
                inputBed = wiseguyCountReference.bedFile,
                outputBed = outputDir + "/references/" + basename(sample.file) + ".gccorrect.bed",
                reference = reference,
                referenceIndex = referenceIndex,
                binFile = binFile,
        }

    }

    call wiseguy.Newref as wiseguyNewref {
        input:
            inputBeds = wiseguyGcCorrectReference.bedFile,
            outputBed = outputDir + "/reference.bed",
            binFile = binFile,
            reference = reference,
            referenceIndex = referenceIndex,
    }

    call bedtools.Sort as bedtoolsSortRef {
        input:
            inputBed = wiseguyNewref.bedFile
    }


    call samtools.BgzipAndIndex as wiseguyReferenceBgzip {
        input:
            type = "bed",
            outputDir = outputDir,
            inputFile = bedtoolsSortRef.bedFile
    }


    # Prepare sample
    call wiseguy.Count as wiseguyCountSample {
        input:
            inputBam = sample.file,
            inputBamIndex = sample.index,
            reference = reference,
            referenceIndex = referenceIndex,
            binFile = binFile,
            outputBed = outputDir + "/" + basename(sample.file) + ".bed"
    }

    call wiseguy.GcCorrect as wiseguyGcCorrectSample {
        input:
            inputBed = wiseguyCountSample.bedFile,
            outputBed = outputDir + "/" + basename(sample.file) + ".gccorrect.bed",
            reference = reference,
            referenceIndex = referenceIndex,
            binFile = binFile,
    }

    call bedtools.Sort as bedtoolsSortSample {
        input:
            inputBed = wiseguyGcCorrectSample.bedFile
    }

    call samtools.BgzipAndIndex as wiseguyGcCorrectSampleIndex {
        input:
            type = "bed",
            outputDir = outputDir,
            inputFile = bedtoolsSortSample.bedFile
    }

    # Calculate zscores

    call wiseguy.Zscore as wiseguyZscore {
        input:
            inputBed = wiseguyGcCorrectSampleIndex.compressed,
            inputBedIndex = wiseguyGcCorrectSampleIndex.index,
            dictionaryFile = wiseguyReferenceBgzip.compressed,
            dictionaryFileIndex = wiseguyReferenceBgzip.index,
            reference = reference,
            referenceIndex = referenceIndex,
            binFile = binFile,
    }

    output {
        File zscores = wiseguyZscore.bedFile
        File wiseguyReference = wiseguyGcCorrectSampleIndex.compressed
        File wiseguyReferenceIndex = wiseguyGcCorrectSampleIndex.index
    }
}