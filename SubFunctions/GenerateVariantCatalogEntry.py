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

print(GenerateVariantCatalogEntry("chr1", (1,100), (25,50), "CAT"))