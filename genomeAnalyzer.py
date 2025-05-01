'''
Analyzes a genome FASTA file.
Outputs:
  - Total genome size.
  - GC content percentage.
  - Relative codon usage by amino acid.
'''

from sequenceAnalysis import NucParams, FastAreader

def main(fileName=None):
    myReader = FastAreader(fileName)
    myNuc = NucParams()

    for head, seq in myReader.readFasta():
        myNuc.addSequence(seq)
    
    # sort codons in alpha order, by Amino Acid
    totalNuc = myNuc.nucCount()
    gcCount = myNuc.nucComposition()['G'] + myNuc.nucComposition()['C']
    gcContent = (gcCount / totalNuc) * 100

    print('sequence length = {:.2f} Mb'.format(totalNuc / 1e6))
    print()
    print('GC content = {:.1f}%'.format(gcContent))
    print()

    # calculate relative codon usage for each codon and print
    codonComp = myNuc.codonComposition()
    aaComp = myNuc.aaComposition()

    aa_codons = {}
    for codon, aa in NucParams.rnaCodonTable.items():
        if aa not in aa_codons:
            aa_codons[aa] = []
        aa_codons[aa].append(codon)

    for aa in sorted(aa_codons.keys()):
        for nuc in sorted(aa_codons[aa]):
            count = codonComp.get(nuc, 0)
            aa_total = aaComp.get(aa, 0) or 1  # Avoid division by zero
            val = count / aa_total
            print('{:s} : {:s} {:5.1f} ({:6d})'.format(nuc, aa, val*100, count))


if __name__ == "__main__":
    main()
    # main('C:/Users/karen/Documents/BME 160/testGenome.fa')

   