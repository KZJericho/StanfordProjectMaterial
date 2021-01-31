import pysam
from Bio.Seq import Seq

def GetReverseComplement(Motif):
    Reverse = Motif.reverse_complement()
    ReverseMotif = str(Reverse)
    return ReverseMotif

#Bounds --> (Lower, Upper)
#Chromosome --> "chr" + chromosome no.
#Motif --> must be in lowercase
def GetCorrectMotif(Chromosome, Bounds, Motif):

    genome = pysam.FastaFile("/home/kevin/Documents/StanfordProjectMaterial/hg19.fa")
    Sequence = genome.fetch(reference=Chromosome, start=Bounds[0], end=Bounds[1])
    Sequence = Sequence.lower()

    Motif = Seq(Motif)
    strMotif = str(Motif)
    MotifRepeat = Sequence.count(strMotif)

    Reverse = GetReverseComplement(Motif)
    ReverseMotif = str(Reverse)
    ReverseMotifRepeat = Sequence.count(ReverseMotif)

    if MotifRepeat > ReverseMotifRepeat:
        return "The correct motif is:" + strMotif
    else:
        return "The correct motif is: " + ReverseMotif

def FindRelativeCoords(Sequence, Motif):
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

    return "The largest perfect " + Motif + " stretch is: " + str(BestLength) + " found at index: " + str(BestLocation)

