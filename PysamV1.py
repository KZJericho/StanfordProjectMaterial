import pysam
from Bio.Seq import Seq

#Bounds --> (Lower, Upper)
#chr = "chr" + chromosome no.
def FindRegion(chromosome, Bounds, Motif):

    genome = pysam.FastaFile("/home/kevin/Documents/StanfordProjectMaterial/hg19.fa")
    Sequence = genome.fetch(reference=chromosome, start=Bounds[0], end=Bounds[1])
    Sequence = Sequence.lower()

    Motif = Seq(Motif)
    MotifRepeat = 0

    ReverseMotif = Motif.reverse_complement()
    ReverseMotifRepeat = 0

    #Motif Repeat Length
    base = 0
    while base < (len(Sequence)-3):
       curr_chunk = Sequence[base:base+len(Motif)]
       if curr_chunk == Motif:
           MotifRepeat += 1
           base += len(Motif)
       else:
           base += 1

    #Reverse Motif Repeat Length
    base2 = 0
    while base2 < (len(Sequence)-3):
       curr_chunk = Sequence[base2:base2+len(ReverseMotif)]
       if curr_chunk == ReverseMotif:
           ReverseMotifRepeat += 1
           base2 += len(ReverseMotif)
       else:
           base2 += 1

    if MotifRepeat > ReverseMotifRepeat:
        return Motif, MotifRepeat
    else:
        return ReverseMotif, ReverseMotifRepeat

'''
IMPROVEMENTS:
    - Combine Motif Repeat and Reverse Motif reading process
    - Use the Pep8 format thing
    - Format the return better
    - More comments
    - combine chromosome and bounds argument?
    - str.count
    - Maybe use Karp-Rabin to speed it up?