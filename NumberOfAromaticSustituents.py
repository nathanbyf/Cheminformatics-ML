#This script calculates the of aromatic substituents and unique aromatic substituents
#You have to execute this script where your casnumbersmiles.csv file is saved

from identifyfunctionalgroups import *
from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem import Mol
import pandas
import numpy as np

# Read csv files containing casnumber and smiles into a dataframe
df = pandas.read_csv('casnumbersmiles.csv', sep =',')

# Add a column numberoffunctionalgroups functionalgroup to the dataframe
df['functionalgroupssubstituents'] = np.nan
df['substituents'] = np.nan

# Get the number of substituents (functional groups)
for (x,row) in df.iterrows():

   functionalgroupssubstituents = []
   substituents = 0

   molecule = Chem.MolFromSmiles(row['smiles'])
   aromaticatoms = Chem.Mol.GetAromaticAtoms(molecule)
   aromaticatomindexes = [atom.GetIdx() for atom in aromaticatoms]

   #Using functionalgroups to calculate the number of substituents
   functionalgroups = identify_functional_groups(molecule)
   print(functionalgroups)
   numberfunctionalgroups = len(functionalgroups)
   if (numberfunctionalgroups > 0):
       for i in range(0, numberfunctionalgroups):
           atomIds = functionalgroups[i].atomIds
           for atomId in atomIds:
               for aromaticatomindex in aromaticatomindexes:
                   aromaticneighborindexes = [neighbor.GetIdx() for neighbor in molecule.GetAtomWithIdx(aromaticatomindex).GetNeighbors()]
                   if atomId in aromaticneighborindexes:
                       #functionalgroupssubstituents = functionalgroupssubstituents + 1
                       functionalgroupssubstituents += [i]

print(functionalgroupssubstituents)

   #Determining substituents by checking if aromaticatom neighbors are not aromatic atoms themself
numberofaromaticindexes = len(aromaticatomindexes)
if (numberofaromaticindexes > 0):
    for z in range(0, numberofaromaticindexes):
        atomId = aromaticatomindexes[z]
        aromaticneighborindexes = [neighbor.GetIdx() for neighbor in molecule.GetAtomWithIdx(atomId).GetNeighbors()]
        for aromaticneighborindex in aromaticneighborindexes:
            if aromaticneighborindex not in aromaticatomindexes:
                substituents = substituents + 1
               #print("aromaticindexes:",aromaticatomindexes, "atomId:",atomId,"typeatomId:",type(atomId),"aromaticneighborindexes:",aromaticneighborindexes,"aromaticneighborindex:",aromaticneighborindex,"substituents:",substituents)



df.at[x,'functionalgroupssubstituents'] = functionalgroupssubstituents
df.at[x,'substituents'] = substituents
df.to_csv('numberofaromaticsubstituents.csv', sep =';',index = False)
