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

Seq1 = "ABCABCABCCATABCABCCATCAT"
Bounds1 = (0,8)
Motif1 = "ABC"
print(ApplyHeuristicRight(Seq1,Bounds1,Motif1))