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

print(FindRegion("chr1",(57223282,57223446), "tttc" ))


def CountSTRs(sequence, motif):
    base = 0
    MotifRepeat = 0
    while base <= (len(sequence)-len(motif)):
        curr_chunk = sequence[base:base+len(motif)]
        if curr_chunk == motif:
            MotifRepeat += 1
            base += len(motif)
        else:
            base +=1
    return MotifRepeat





'''
IMPROVEMENTS:
    - Combine Motif Repeat and Reverse Motif reading process
    - Use the Pep8 format thing
    - Format the return better
    - More comments
    - combine chromosome and bounds argument?
    - str.count DONE
    - Maybe use Karp-Rabin to speed it up? nah
'''