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
import "tasks/wisestork.wdl" as wisestork
import "tasks/samtools.wdl" as samtools
import "tasks/bedtools.wdl" as bedtools

workflow wisestorkCnv {
    input {
    # IndexedBamFile from common.wdl
    Array[IndexedBamFile] controlBams
    IndexedBamFile case
    File reference
    File referenceIndex
    File? binFile
    String outputDir = "."
    }

    # Prepare reference
    scatter (controlBam in controlBams) {
        call wisestork.Count as wisestorkCountControls {
            input:
                inputBam = controlBam.file,
                inputBamIndex = controlBam.index,
                reference = reference,
                referenceIndex = referenceIndex,
                binFile = binFile,
                outputBed = outputDir + "/references/" + basename(controlBam.file) + ".bed"
        }

        call wisestork.GcCorrect as wisestorkGcCorrectControls {
            input:
                inputBed = wisestorkCountControls.bedFile,
                outputBed = outputDir + "/references/" + basename(controlBam.file) + ".gccorrect.bed",
                reference = reference,
                referenceIndex = referenceIndex,
                binFile = binFile,
        }

    }

    call wisestork.Newref as wisestorkNewref {
        input:
            inputBeds = wisestorkGcCorrectControls.bedFile,
            outputBed = outputDir + "/reference.bed",
            binFile = binFile,
            reference = reference,
            referenceIndex = referenceIndex,
    }

    call bedtools.Sort as bedtoolsSortRef {
        input:
            inputBed = wisestorkNewref.bedFile,
            outputBed = outputDir + "/reference.sorted.bed"
    }


    call samtools.BgzipAndIndex as wisestorkReferenceBgzip {
        input:
            type = "bed",
            outputDir = outputDir,
            inputFile = bedtoolsSortRef.bedFile
    }


    # Prepare sample
    call wisestork.Count as wisestorkCountCase {
        input:
            inputBam = case.file,
            inputBamIndex = case.index,
            reference = reference,
            referenceIndex = referenceIndex,
            binFile = binFile,
            outputBed = outputDir + "/" + basename(case.file) + ".bed"
    }

    call wisestork.GcCorrect as wisestorkGcCorrectCase {
        input:
            inputBed = wisestorkCountCase.bedFile,
            outputBed = outputDir + "/" + basename(case.file) + ".gccorrect.bed",
            reference = reference,
            referenceIndex = referenceIndex,
            binFile = binFile,
    }

    call bedtools.Sort as bedtoolsSortCase {
        input:
            inputBed = wisestorkGcCorrectCase.bedFile,
            outputBed = outputDir + "/sample.gccorrect.sorted.bed"
    }

    call samtools.BgzipAndIndex as wisestorkGcCorrectCaseIndex {
        input:
            type = "bed",
            outputDir = outputDir,
            inputFile = bedtoolsSortCase.bedFile
    }

    # Calculate zscores

    call wisestork.Zscore as wisestorkZscore {
        input:
            inputBed = wisestorkGcCorrectCaseIndex.compressed,
            inputBedIndex = wisestorkGcCorrectCaseIndex.index,
            dictionaryFile = wisestorkReferenceBgzip.compressed,
            dictionaryFileIndex = wisestorkReferenceBgzip.index,
            reference = reference,
            referenceIndex = referenceIndex,
            binFile = binFile,
            outputBed = outputDir + "/" + basename(case.file)+ ".zscores.bed"
    }

    output {
        File zscores = wisestorkZscore.bedFile
        File wisestorkReference = wisestorkGcCorrectCaseIndex.compressed
        File wisestorkReferenceIndex = wisestorkGcCorrectCaseIndex.index
    }
}