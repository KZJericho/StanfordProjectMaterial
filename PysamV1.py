import pysam
from Bio.Seq import Seq

#Bounds --> (Lower, Upper)
#chr = "chr" + chromosome no.
#motif --> must be in lowercase
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
        return ReverseMotif, CountSTRs(Sequence, ReverseMotif)

print(FindRegion("chr1",(57222577,57224042), "aaag" ))


def CountSTRs(sequence, motif):
    base = 0
    LowerBound = 0
    HigherBound = 0
    MotifRepeat = 0
    while base <= (len(sequence)-len(motif)):
        CurrChunk = sequence[base:base+len(motif)]

        if CurrChunk == motif:
            MotifRepeat += 1
            LowerBound = base - (len(motif)*MotifRepeat)
            HigherBound == base + len(motif) - 1
            base += len(motif)
        else:
            base +=1


    return MotifRepeat, base





'''
IMPROVEMENTS:
    - Combine Motif Repeat and Reverse Motif reading process
    - Use the Pep8 format
    - Format the return better
    - More comments
    - Should HigherBound be without the -1?
'''