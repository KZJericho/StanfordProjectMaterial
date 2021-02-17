import pysam

#Bounds --> (Lower, Upper)
#Chromosome --> "chr" + chromosome no.
def ExtractGenome(Chromosome, Bounds):
    genome = pysam.FastaFile("/home/kevin/Documents/StanfordProjectMaterial/hg19.fa")
    #start bound is minus 1 to account for closed-base system
    Sequence = genome.fetch(reference=Chromosome, start=Bounds[0]-1, end=Bounds[1])
    Sequence = Sequence.lower()

    return Sequence

def GetReverseComplement(Motif):
    ConvertDict = {"a": "t", "c": "g", "g": "c", "t": "a"}
    return "".join([ConvertDict[base] for base in Motif[::-1]])

#Motif --> must be in lowercase
def GetCorrectMotif(Sequence, Motif):

    MotifRepeat = Sequence.count(Motif)

    ReverseMotif = GetReverseComplement(Motif)
    ReverseMotifRepeat = Sequence.count(ReverseMotif)

    if MotifRepeat > ReverseMotifRepeat:
        #"The correct motif is:" +
        return strMotif
    else:
        return ReverseMotif

def FindPerfectCoords(Sequence, Motif):
    CurrLocation = 0
    BestLocation = 0

    CurrLength = 0
    BestLength = 0

    base = 0
    while base <= (len(Sequence)):
        CurrChunk = Sequence[base:base+len(Motif)]

        if CurrChunk == Motif:
            CurrLength += 1
            CurrLocation = base

            base += len(Motif)

        else:
            #how to approach equal repeats?
            if CurrLength > BestLength:
                BestLength = CurrLength
                BestLocation = CurrLocation - (len(Motif))*((CurrLength)-1)

            CurrLength = 0
            base += 1

    return (BestLocation,((BestLocation+((BestLength)*(len(Motif))))))
'''
FOR TESTING CODE ABOVE
def FindRegion(Chromosome, Bounds, Motif):
    ReverseComplement = GetReverseComplement(Motif)
    print("RevComp=", ReverseComplement)
    Sequence = ExtractGenome(Chromosome,Bounds)
    print("Sequence:", Sequence)
    CorrectMotif = GetCorrectMotif(Sequence,Motif)
    print("CorrectMotif:", CorrectMotif)

    return FindPerfectCoords(Sequence,CorrectMotif)
print(FindRegion("chr1",(57222577,57224042), "aaag"))

def FindPerfectLength(Sequence, Motif):
    res = FindPerfectCoords(Sequence,Motif)
    return (res[1]-res[0])+1
'''
def NextKmerLeft(Sequence, LowerBound, CorrectMotif):
        NextKmerLower = LowerBound - len(CorrectMotif)
        if NextKmerLower < 0:
            return None
        else:
            return Sequence[NextKmerLower:NextKmerLower+len(CorrectMotif)]

def ApplyHeuristicLeft(Sequence, PerfectBounds, CorrectMotif):
    LowerBound = PerfectBounds[0]

    while LowerBound >= 0:
        NextKmer = NextKmerLeft(Sequence, LowerBound, CorrectMotif)

        if NextKmer == CorrectMotif:
            NextKmer = NextKmerLeft(Sequence, LowerBound, CorrectMotif)
            LowerBound -= len(CorrectMotif)
        else:
            Minus2Bound = LowerBound -(len(CorrectMotif))
            if Minus2Bound < 0:
                return (LowerBound, PerfectBounds[1])
            KmerMinus2 = NextKmerLeft(Sequence,Minus2Bound,CorrectMotif)
            if KmerMinus2 == CorrectMotif:
                LowerBound -= (len(CorrectMotif)*(2))
            else:
                return (LowerBound, PerfectBounds[1])
    return (LowerBound, PerfectBounds[1])

def NextKmerRight(Sequence, UpperBound, CorrectMotif):
        NextKmerUpper = UpperBound + 1
        if NextKmerUpper >= len(Sequence):
            return None
        else:
            return Sequence[NextKmerUpper:NextKmerUpper+len(CorrectMotif)]

def ApplyHeuristicRight(Sequence, PerfectBounds, CorrectMotif):
    UpperBound = PerfectBounds[1]

    while UpperBound < len(Sequence):
        NextKmer = NextKmerRight(Sequence, UpperBound, CorrectMotif)

        if NextKmer == CorrectMotif:
            NextKmer = NextKmerRight(Sequence, UpperBound, CorrectMotif)
            UpperBound += len(CorrectMotif)
        else:
            Plus2Bound = UpperBound + (len(CorrectMotif))
            if Plus2Bound < 0:
                return (PerfectBounds[0], UpperBound)
            KmerPlus2 = NextKmerRight(Sequence,Plus2Bound,CorrectMotif)
            if KmerPlus2 == CorrectMotif:
                UpperBound += (len(CorrectMotif)*(2))
            else:
                return (PerfectBounds[0], UpperBound)
    return (PerfectBounds[0], UpperBound)

import json
import pprint

def GenerateVariantCatalogEntry(Chromosome, InitialBounds, FinalBounds, CorrectMotif):

    res = (
            f"VariantType: Repeat", \
            f"LocusId: {Chromosome}_{InitialBounds[0]}_{InitialBounds[1]}", \
            f"LocusStructure: ({CorrectMotif})*", \
            f"ReferenceRegion: {Chromosome}:{FinalBounds[0]}-{FinalBounds[1]}"

            )
    output = json.dumps(res, indent=4)
    return output

def FindFinalResult(Chromosome, InitialBounds, Motif):
    Sequence = ExtractGenome(Chromosome, InitialBounds)
    print("Sequence Identified")
    CorrectMotif = GetCorrectMotif(Sequence, Motif)
    print(f"The correct motif is {CorrectMotif}")
    PerfectCoords = FindPerfectCoords(Sequence, CorrectMotif)
    print(f"The longest {CorrectMotif} stretch coordinates are: {PerfectCoords}")
    ExpandedCoordsLeft = ApplyHeuristicLeft(Sequence, PerfectCoords, CorrectMotif)
    ExpandedCoordsRight = ApplyHeuristicRight(Sequence, ExpandedCoordsLeft, CorrectMotif)

    return GenerateVariantCatalogEntry(Chromosome, InitialBounds, ExpandedCoordsRight, CorrectMotif)

print(FindFinalResult("chr1",(57222577,57224042), "aaag"))


'''
TO-DO:
- Pep8 styles
- multiple catalog entries at once?
'''
