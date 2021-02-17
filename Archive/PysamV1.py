import pysam
from Bio.Seq import Seq

#Bounds --> (Lower, Upper)
#chr = "chr" + chromosome no.
#motif --> must be in lowercase
def CountSTRs(sequence, motif):
    base = 0
    CurrLocation = 0
    BestLocation = 0
    CurrRepeat = 0
    BestRepeat = 0

    while base <= (len(sequence)-len(motif)):
        CurrChunk = sequence[base:base+len(motif)]

        if CurrChunk == motif:
            CurrRepeat += 1
            CurrLocation = base

            base += len(motif)

        else:
            if CurrRepeat > BestRepeat:
                BestRepeat = CurrRepeat
                BestLocation = CurrLocation - (len(motif))*((CurrRepeat)-1)

            CurrRepeat = 0
            base += 1

    if CurrRepeat > BestRepeat:
        BestRepeat = CurrRepeat
        BestLocation = CurrLocation - (len(motif))*((CurrRepeat)-1)


    return ["length", BestRepeat, "location", BestLocation, (BestLocation+((len(motif)*(BestRepeat))-1))]

'''
def FindFinalCoordinates(sequence, bounds, motif):
    result = CountSTRs(sequence, motif)
    result[3] += bounds[0]
    result[4] += bounds[1]

    return result
'''
def FindRegion(chromosome, Bounds, Motif):

    genome = pysam.FastaFile("/home/kevin/Documents/StanfordProjectMaterial/hg19.fa")
    Sequence = genome.fetch(reference=chromosome, start=Bounds[0], end=Bounds[1])
    Sequence = Sequence.lower()

    Motif = Seq(Motif)
    strMotif = str(Motif)
    MotifRepeat = Sequence.count(strMotif)

    Reverse = Motif.reverse_complement()
    ReverseMotif = str(Reverse)
    ReverseMotifRepeat = Sequence.count(ReverseMotif)

    if MotifRepeat > ReverseMotifRepeat:
        return Motif, CountSTRs(Sequence, Motif)
    else:
        return ReverseMotif, CountSTRs(Sequence,ReverseMotif)

print(FindRegion("chr1",(57222577,57224042), "aaag" ))



'''
IMPROVEMENTS:
    - Combine Motif Repeat and Reverse Motif reading process
    - Use the Pep8 format
    - Format the return better
    - More comments
    - Should HigherBound be without the -1?
'''