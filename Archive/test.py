def NextKmerLeft(Sequence, LowerBound, CorrectMotif):
        NextKmerLower = LowerBound - len(CorrectMotif)
        if NextKmerLower < 0:
            return None
        else:
            return Sequence[NextKmerLower:NextKmerLower+len(CorrectMotif)]

def ApplyHeuristicLeft(Sequence, PerfectBounds, CorrectMotif):
    LowerBound = PerfectBounds[0]
    print("LowerBound:", LowerBound)

    while LowerBound >= 0:
        NextKmer = NextKmerLeft(Sequence, LowerBound, CorrectMotif)
        print("NextKmer:", NextKmer)

        if NextKmer == CorrectMotif:
            print("perfectmatch")
            print(LowerBound)
            NextKmer = NextKmerLeft(Sequence, LowerBound, CorrectMotif)
            LowerBound -= len(CorrectMotif)
        else:
            print("here0")
            Minus2Bound = LowerBound -(len(CorrectMotif))
            print("Minus2Bound:", Minus2Bound)
            if Minus2Bound < 0:
                return "here1", (LowerBound, PerfectBounds[1])
            KmerMinus2 = NextKmerLeft(Sequence,Minus2Bound,CorrectMotif)
            print("KmerMinus2:", KmerMinus2)
            if KmerMinus2 == CorrectMotif:
                print("Match")
                LowerBound -= (len(CorrectMotif)*(2))
                print("LowerBound:", LowerBound)
            else:
                return "here2", (LowerBound, PerfectBounds[1])
    return "here3", (LowerBound, PerfectBounds[1])

Seq1 = "CATCATABCABCCTGABCCATABCABCABC"
Bounds1 = (21,29)
Motif1 = "ABC"
print(ApplyHeuristicLeft(Seq1,Bounds1,Motif1))

'''
def NextKmer(give next motif on the left or none)
#ask for next motif (kmer) (the repeat before)
#function either returns next kmer or none if next to bound

#check if Kmer is perfect match
#if yes, find next Kmer
    #check if next Kmer is perfect match (will only happen after first iteration)
#if no, find next Kmer
    #if next Kmer (n-2) is perfect, find Next Kmer and add
    #if not perfect terminate
#repeat and repeat recursively
'''