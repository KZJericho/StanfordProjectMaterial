import pysam
from Bio.Seq import Seq

def GetReverseComplement(Motif):
    Reverse = Motif.reverse_complement()
    ReverseMotif = str(Reverse)
    return ReverseMotif

#Bounds --> (Lower, Upper)
#Chromosome --> "chr" + chromosome no.
def ExtractGenome(Chromosome, Bounds):
    genome = pysam.FastaFile("/home/kevin/Documents/StanfordProjectMaterial/hg19.fa")
    Sequence = genome.fetch(reference=Chromosome, start=Bounds[0], end=Bounds[1])
    Sequence = Sequence.lower()

    return Sequence

#Motif --> must be in lowercase
def GetCorrectMotif(Sequence, Motif):

    Motif = Seq(Motif)
    strMotif = str(Motif)
    MotifRepeat = Sequence.count(strMotif)

    Reverse = GetReverseComplement(Motif)
    ReverseMotif = str(Reverse)
    ReverseMotifRepeat = Sequence.count(ReverseMotif)

    if MotifRepeat > ReverseMotifRepeat:
        #"The correct motif is:" +
        return strMotif
    else:
        return ReverseMotif

def FindPerfectLength(Sequence, Motif):
    CurrLocation = 0
    BestLocation = 0

    CurrLength = 0
    BestLength = 0

    base = 0
    while base <= (len(Sequence)-len(Motif)):
        CurrChunk = Sequence[base:base+len(Motif)]

        if CurrChunk == Motif:
            CurrLength += 1
            CurrLocation = base
            base += len(Motif)

        else:
            #how to approach equal repeats?
            if CurrLength >= BestLength:
                BestLength = CurrLength
                BestLocation = CurrLocation - len(Motif)*((CurrLength)-1)
                CurrLength = 0

            base += 1
    #The longest [insert motif] stretch is:
    return BestLength

def FindPerfectCoords(Sequence, Motif):
    CurrLocation = 0
    BestLocation = 0

    CurrLength = 0
    BestLength = 0

    base = 0
    while base <= (len(Sequence)-len(Motif)):
        CurrChunk = Sequence[base:base+len(Motif)]

        if CurrChunk == Motif:
            CurrLength += 1
            CurrLocation = base
            base += len(Motif)

        else:
            #how to approach equal repeats?
            if CurrLength >= BestLength:
                BestLength = CurrLength
                BestLocation = CurrLocation - len(Motif)*((CurrLength)-1)
                CurrLength = 0

            base += 1
    #The longest [insert motif] stretch coordinates are:
    return (BestLocation, BestLocation+((len(Motif)*BestLength)-1))

#PerfectBounds --> (Lower, Upper)
def ApplyHeuristicLeft(Sequence, PerfectBounds, CorrectMotif):
    CurrChunk = Bounds[0] - len(Motif)
    MotifChunk = 0

    base = 0
    while base <= (len(Sequence)-len(Motif)):

        count = 0
        while CurrChunk >= 0 and MotifChunk < len(Motif):
            if Sequence[CurrChunk] == Motif[MotifChunk]:
                print("here", Sequence[CurrChunk])
                count += 1

            CurrChunk -= 1
            MotifChunk += 1

        if count >= len(Motif) + 1:
            Bounds[0] -= len(Motif)
            count = 0

        base += len(Motif)
    return Bounds

#print(ApplyHeuristicLeft("ABCABBABCABC", (6,11), "ABC"))

def ApplyHeuristicRight(Sequence, PerfectBounds, CorrectMotif):
    pass

import json
def GenerateVariantCatalogEntry(Chromosome, InitialBounds, FinalBounds, CorrectMotif):
    result = {
        "VariantType": "Repeat",
        "LocusId": str(Chromosome) + "_" + str(InitialBounds[0]) + "_" + str(InitialBounds[1]),
        "LocusStructure": str(CorrectMotif) + "*",
        "ReferenceRegion": str(Chromosome) + ":" + str(FinalBounds[0]) + "-" + str(FinalBounds[1]),
    }
    output = json.dumps(result, indent=4)
    return json.dumps(result)


#print(GenerateVariantCatalogEntry("chr1", (0,10), (5,7), "ABC"))

def FindFinalResult(Chromosome, Bounds, Motif):
    Sequence = ExtractGenome(Chromosome, Bounds)
    CorrectMotif = GetCorrectMotif(Sequence, Motif)
    print(CorrectMotif)
    PerfectLength = FindPerfectLength(Sequence, CorrectMotif)
    PerfectCoords = FindPerfectCoords(Sequence, CorrectMotif)
    ExpandedCoordsLeft = ApplyHeuristic(Sequence, PerfectCoords, CorrectMotif)
    ExpandedCoordsRight = ApplyHeuristicRight(Sequence, ExpandedCoordsLeft, CorrectMotif)
    return ExpandedCoordsRight

#print(FindFinalResult("chr1",(57222577,57224042), "aaag"))


'''
TO-DO:
- Fix the coordinate inaccuracy with perfect coords
- Figure out the heuristic addition --> the loops are messy
- Figure out how to add spaces to the GenerateVariantCatalogEntry
- Add print statements to FindFinalResult
- Pep8 style
'''
