#!/usr/bin/env python3
# Name: Karen Wang (1551065)
# Group Members: Catherine Gholamipour, Autumn Seda Kong, jooYoung Gwag, Keerthana Ande

class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
      '''
      Initialize with a protein sequence, count valid amino acids.
      '''
      fixedProtein = protein.upper().strip()
      self.protein = fixedProtein
      self.proteinDict = {}
      for aa in ProteinParam.aa2mw.keys():
        self.proteinDict[aa] = 0

      for aa in fixedProtein:
        if aa in self.proteinDict:
          self.proteinDict[aa] += 1
        pass

    def aaCount (self):
      '''
      Return the total count of valid amino acids in the sequence.
      '''
      return sum(self.proteinDict.values())
      pass

    def pI (self):
      '''
      Estimate isoelectric point by scanning pH 0.00 to 14.00 in 0.01 steps.
      '''
      bestpH = 0.0
      for i in range (1, 1401):
        pH = i / 100
        if abs(self._charge_(pH)) < abs(self._charge_(bestpH)):
          bestpH = pH
      return round(bestpH, 2)
      pass

    def aaComposition (self) :
      '''
      Return a dictionary of amnio acid counts for the squence.
      '''
      return self.proteinDict
      pass

    def _charge_ (self, pH):
      '''
      Calculate the net charge at a specific pH for the protein sequence.
      '''
      posCharge = (10**ProteinParam.aaNterm) / (10**ProteinParam.aaNterm + 10**pH)
      for aa, pKa in ProteinParam.aa2chargePos.items():
        posCharge += self.proteinDict[aa] * (10**pKa) / (10**pKa + 10**pH)
      
      negCharge = (10**pH) / (10**pH + 10**ProteinParam.aaCterm)
      for aa, pKa in ProteinParam.aa2chargeNeg.items():
        negCharge += self.proteinDict[aa] * (10**pH) / (10**pH + 10**pKa)

      return posCharge - negCharge
      pass

    def molarExtinction (self, cystine=True):
      '''
      Calculate the molar extinction coefficient at 280nm for the protein sequence.
      '''
      Y = self.proteinDict['Y']
      W = self.proteinDict['W']
      C = self.proteinDict['C']
      if cystine:
        Cys = C // 2
      else:
        Cys = 0
      return Y * ProteinParam.aa2abs280['Y'] + W * ProteinParam.aa2abs280['W'] + Cys * ProteinParam.aa2abs280['C']
      pass

    def massExtinction (self, cystine=True):
      ''' 
      Calculate the mass extinction coefficient for the protein sequence.
      '''
      myMW =  self.molecularWeight()
      myME = self.molarExtinction(cystine=cystine)
      return myME / myMW if myMW else 0.0

    def molecularWeight (self):
      ''' 
      Return the total molecular weight of the protein.
      '''
      if self.aaCount() == 0:
        return 0.0
      total = sum(count*ProteinParam.aa2mw[aa] for aa, count in self.proteinDict.items())
      total -= (self.aaCount() - 1) * ProteinParam.mwH2O
      return total
      pass

# Please do not modify any of the following.  This will produce a standard output that can be parsed

import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")

        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present

        for aa,n in sorted(myParamMaker.aaComposition().items(),
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))

        inString = input('protein sequence?')

if __name__ == "__main__":
    main()