'''
compareGenomes.py

Compare GC content, amnio acid composition, and relative codon bias between
two genome FASTA files (halophile and hyperthermophile)
'''

from sequenceAnalysis import NucParams, FastAreader

def analyze (fileName=None):
  reader = FastAreader(fileName)
  nuc = NucParams()

  for head, seq in reader.readFasta():
    nuc.addSequence(seq)

  totalNuc = nuc.nucCount()
  gcCount = nuc.nucComposition()['G'] + nuc.nucComposition()['C']
  gcContent = (gcCount / totalNuc) * 100
  
  return gcContent, nuc.aaComposition(), nuc.codonComposition()


def main(genome1, genome2):
  gc1, aa1, codon1 = analyze(genome1)
  gc2, aa2, codon2 = analyze(genome2)

  print('GC content comparison:')
  print()
  print(f'{genome1}: {gc1:.2f}%')
  print(f'{genome2}: {gc2:.2f}%\n')

  print('Amino acid composition comparison:')
  print()
  print(f'{genome1}: {aa1}')
  print(f'{genome2}: {aa2}\n')

  print('Relative codon bias comparison:')
  print()
  print(f'{genome1}: {codon1}')
  print(f'{genome2}: {codon2}\n')

if __name__ == "__main__":
  main()
  # main("C:/Users/karen/Documents/BME 160/haloVolc1_1-genes.fa", "C:/Users/karen/Documents/BME 160/testGenome.fa")


