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
        return Motif
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
        NextKmerUpper = UpperBound
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

def FinalCoords(InitialBounds, FinalBounds):
    LowerFinal = FinalBounds[0] + InitialBounds[0]
    UpperFinal = FinalBounds[1] + InitialBounds[0] - 1
    return (LowerFinal,UpperFinal)

def GenerateVariantCatalogEntry(Chromosomes, InitBounds, FinBounds, CorrMotif):

    list_dicts = []
    for Chromosome, InitialBounds, FinalBounds, CorrectMotif in zip(Chromosomes, InitBounds, FinBounds, CorrMotif):
        res = {}
        res["VariantType"] = "Repeat"
        res["LocusId"] = "{0}_{1}_{2}".format(Chromosome, InitialBounds[0], InitialBounds[1])
        res["LocusStructure"] = "({0})*".format(CorrectMotif.upper())
        res["ReferenceRegion"] = "{0}:{1}-{2}".format(Chromosome,FinalBounds[0],FinalBounds[1])
        list_dicts.append(res)
    with open('/home/kevin/Documents/StanfordProjectMaterial/results.json', 'w') as outfile:
        json.dump(list_dicts, outfile, indent=4)
    return "done"

def FindFinalResult(Chromosomes, InitialBounds, Motifs):
    FinalCoordsList = []
    CorrectMotifList = []
    Count = len(Chromosomes)
    i = 0
    while i < Count:
        Chromosome = Chromosomes[i]
        InitialBound = InitialBounds[i]
        Motif = Motifs[i]
        print(InitialBound)
        print("here")
        Motif = Motif.lower()
        Sequence = ExtractGenome(Chromosome, InitialBound)
        #print("Sequence Identified")
        CorrectMotif = GetCorrectMotif(Sequence, Motif)
        #print(f"The correct motif is {CorrectMotif}")
        PerfectCoords = FindPerfectCoords(Sequence, CorrectMotif)
        #print(f"The longest {CorrectMotif} stretch coordinates are: {FinalCoords(InitialBounds,PerfectCoords)}")
        ExpandedCoordsLeft = ApplyHeuristicLeft(Sequence, PerfectCoords, CorrectMotif)
        ExpandedCoordsRight = ApplyHeuristicRight(Sequence, ExpandedCoordsLeft, CorrectMotif)
        FinalCoord = FinalCoords(InitialBound, ExpandedCoordsRight)
        FinalCoordsList.append(FinalCoord)
        CorrectMotifList.append(CorrectMotif)
        i += 1
    return GenerateVariantCatalogEntry(Chromosomes, InitialBounds, FinalCoordsList, CorrectMotifList)

#print(FindFinalResult(["chrX", "chrX"],[(1642529	,1644412), (1642529	,1644412)], ["ATCC", "ATCC"]))


import pandas as pd

def RunTest():
    Chromosomes = []
    InitialBounds = []
    Motifs = []

    df = pd.read_csv('/home/kevin/Documents/StanfordProjectMaterial/CancerLoci30Annotations.bed', header=None, sep='\t')
    df.columns = ["Chr","Start", "End", "STR", ""]
    df = df.drop(columns=[""])

    Chromosomes = df["Chr"].tolist()

    Start = df["Start"].tolist()

    End = df["End"].tolist()

    Motifs = df["STR"].tolist()

    Bounds = []

    i = 0
    while i < len(Start):
        Bounds.append((Start[i],End[i]))
        i += 1

    FindFinalResult(Chromosomes, Bounds, Motifs)

RunTest()


#with open("v094_all_filtered_repeatloci_merged.bed") as file:
#    lines = csv.reader(file, delimiter = “\t”)
#    for line in lines:



'''
TO-DO:
- Pep8 styles
- multiple catalog entries at once?
- do kmer minus 2 to be accepted
- look into first and second longest being acceptable
- variants allowed for longer repeat lengths

'''
