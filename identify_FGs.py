#This script calculates the of aromatic substituents and unique aromatic substituents
#You have to execute this script where your casnumbersmiles.csv file is saved

from identifyfunctionalgroups import *
from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem import Mol
import pandas
import numpy as np
import re

# Read csv files containing casnumber and smiles into a dataframe
df = pandas.read_csv('corrections.csv', sep =',')

# Add a column numberoffunctionalgroups functionalgroup to the dataframe
df['unique_substituents'] = ''
df['nonfunctionalsubstituents'] = np.nan

#df['functionalgroups'] = np.nan
# 1. Add column nitrogens
df['numberofnitrogens'] = np.nan

# 2. Add column halogens
df['numberofhalogens'] = np.nan

# 3. Add column numberofrings
df['numberofrings'] = np.nan

# 4. Add column quarttertcarbons
df['tertsquartscarbonatoms'] = np.nan

# 5. Add column heteroatomsinrings
df['numberofheteroatomsinrings'] = np.nan

# 6. Add a column aromaticsubstituents to the dataframe
#df['aromaticsubstituents'] = np.nan

# 7. Add a column functionalgroups to the dataframe
df['numberfunctionalgroups'] = np.nan

# 8. Add a column carbonyloxygens
#df['carbonylgroupnocarboxylic'] = np.nan

# 9. Add a column noncarbonyloxygens
df['numberofnoncarbonyloxygens'] = np.nan



# Get the number of substituents (functional groups)
for (x,row) in df.iterrows():

   nonfunctionalsubstituents = 0
   #read the smiles and save it as molecule
   molecule = Chem.MolFromSmiles(row['smiles'])
   allatomindexes = [atom.GetIdx() for atom in molecule.GetAtoms()]
   #perform GetAromaticAatoms func. on the molcule to see which atoms are in a ring
   aromaticatoms = Chem.Mol.GetAromaticAtoms(molecule)
   #Get each atom's index, index is given only to aromatic atoms
   #The code searches through the variable: aromaticatoms, since the aromatic atoms are saved as this
   #aromaticatomindexes = [atom.GetIdx() for atom in aromaticatoms]
   aromaticatomindexes = [atom.GetIdx() for atom in aromaticatoms]
   nonaromaticatomindexes = [z for z in allatomindexes if z not in aromaticatomindexes]
   functionalgroups = identify_functional_groups(molecule)
   functionalgroupindexes = [group.atomIds for group in functionalgroups]
   flattenedfunctionalgroupindexes = [item for sublist in functionalgroupindexes for item in sublist]
   uniquesubstituents_types = set()

   for idx in aromaticatomindexes:
       neighbourindexes = [neighbor.GetIdx() for neighbor in molecule.GetAtomWithIdx(idx).GetNeighbors()]
       nonaromaticneighbourindexes = set(neighbourindexes).intersection(nonaromaticatomindexes)
       if len(nonaromaticneighbourindexes)>0:
           nonaromaticfunctionalneighbourindexes = set(nonaromaticneighbourindexes).intersection(flattenedfunctionalgroupindexes)
           if len(nonaromaticfunctionalneighbourindexes)>0:
               for ind in nonaromaticfunctionalneighbourindexes:
                   for i, groupindexes in enumerate(functionalgroupindexes):
                       if ind in groupindexes:
                           uniquesubstituents_types.add(functionalgroups[i][2])
           else:
               nonfunctionalsubstituents+=1

   df.at[x,'unique_substituents'] = len(uniquesubstituents_types)+nonfunctionalsubstituents
   df.at[x,'nonfunctionalsubstituents'] = nonfunctionalsubstituents

   numberofaromaticindexes = len(aromaticatomindexes)
   if (numberofaromaticindexes > 0):
       for z in range(0, numberofaromaticindexes):
           atomId = aromaticatomindexes[z]
           aromaticneighborindexes = [neighbor.GetIdx() for neighbor in molecule.GetAtomWithIdx(atomId).GetNeighbors()]
           for aromaticneighborindex in aromaticneighborindexes:
               if aromaticneighborindex not in aromaticatomindexes:
                   nonfunctionalsubstituents = nonfunctionalsubstituents + 1
                   #print("aromaticindexes:",aromaticatomindexes, "atomId:",atomId,"typeatomId:",type(atomId),"aromaticneighborindexes:",aromaticneighborindexes,"aromaticneighborindex:",aromaticneighborindex,"substituents:",substituents)

   def label(a):
       return 100*int(a.GetHybridization() + a.GetAtomicNum())

   nms = [Chem.Mol(f) for f in molecule]
   for nm in nms:
       for at in nm.GetAtoms():
           at.SetIsotope(label(at))
   mcs=rdFMCS.FindMCS(nms,atomCompare=rdFMCS.AtomCompare.CompareIsotopes)
   print(mcs.smartsString)





   #Number of carbonyl oxygens, number of non-carbonyl oxygens
   count = 0
   molecule = Chem.MolFromSmiles(row['smiles'])
   functionalgroups = identify_functional_groups(molecule)

   for group in functionalgroups:

       if group[1] == 'C=O':
            count+=1


    #noncarbonyl oxygens
   oxygen = 0

   for atom in molecule.GetAtoms():
      atomname = atom.GetSymbol()
      if atomname == 'O' or atomname == 'o':
          oxygen +=1

    #print(functionalgroups[2])
    #df.at[x, 'functionalgroups'] = functionalgroups[2]
   df.at[x,'numberofcarbonyloxygens'] = count
   df.at[x,'numberofnoncarbonyloxygens'] = oxygen - count

 #1. Number of nitrogens
   numberofnitrogens = 0
   for atom in molecule.GetAtoms():
       element = atom.GetSymbol()
       if element == 'N':
           numberofnitrogens +=1
   df.at[x,'numberofnitrogens'] = numberofnitrogens

# 2. Add column halogens
   numberofhalogens = 0
   for atom in molecule.GetAtoms():
       element = atom.GetSymbol()
       if element == 'F' or element == 'Br' or element == 'Cl':
           numberofhalogens +=1
   df.at[x,'numberofhalogens'] = numberofhalogens

# 3. Calculate number of rings
   ssr = Chem.GetSSSR(molecule)
   df.at[x,'numberofrings'] = ssr


# 4. Count the tertiarycarbonatoms and quarternarycarbonatoms
   tertiarycarboncount = 0
   quarternarycarboncount = 0

	#Loop over the atoms int the molecule
   for atom in molecule.GetAtoms():

		#Initialize counters
        singlecount = 0
        aromaticcount = 0
        doublecount = 0

		#Get atom properties
        atomindex = atom.GetIdx()
        neighborindexes = [neighbor.GetIdx() for neighbor in atom.GetNeighbors()]
        neighbornames = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]

		#When atom is a carbon atom, count carbon neighbors that are single bonded, double bonded
        if(atom.GetSymbol() == 'C'):
            for neighbor in atom.GetNeighbors():
                bondname = str(molecule.GetBondBetweenAtoms(atomindex,neighbor.GetIdx()).GetBondType())
                name = neighbor.GetSymbol()
                if(name == 'C' and bondname == 'SINGLE'):
                    singlecount = singlecount + 1
                    #print("singlecount:",singlecount, "for neighbor", name)
                if(name == 'C' and bondname == 'AROMATIC'):
                    aromaticcount = aromaticcount + 1
                    #print("singlecount:",singlecount, "for neighbor", name)
                if(name == 'C' and bondname == 'DOUBLE'):
                    doublecount = doublecount + 1
                    #print("doublecount:",doublecount, "for neighbor", name)
        valence = molecule.GetAtomWithIdx(atomindex).GetExplicitValence()
        if(valence == 4 and singlecount == 4):
            quarternarycarboncount = quarternarycarboncount + 1
        if(valence == 4 and singlecount == 3):
            tertiarycarboncount = tertiarycarboncount + 1
        if(valence == 4 and singlecount == 2 and doublecount == 1):
            tertiarycarboncount = tertiarycarboncount + 1
        if(valence == 4 and singlecount == 1 and aromaticcount == 2):
            tertiarycarboncount = tertiarycarboncount + 1
        #print(row['casnumber'],atomindex,atom.GetSymbol(),neighborindexes,neighbornames,"singlecount:",singlecount,"aromaticcount:",aromaticcount,"doublecount:",doublecount)

    #print("tertiarycarboncount:",tertiarycarboncount,"quarternarycarboncount:",quarternarycarboncount)
   df.at[x,'tertsquartscarbonatoms'] = tertiarycarboncount + quarternarycarboncount

# 7. Number of heteroatomsinrings

   heterocount = 0

   #Looping through atoms in molecule
   for atom in molecule.GetAtoms():
       atomisinring = atom.IsInRing()
       atomname = atom.GetSymbol()
       if (atomisinring == True and atomname != 'C' and atomname != 'H'):
           heterocount += 1
        #print(Chem.MolToSmiles(molecule), "after atomname", atomname, "count is:", count)
   df.at[x,'numberofheteroatomsinrings'] = heterocount

# 7. Get the number of functional groups

   functionalgroups = identify_functional_groups(molecule)
   numberofalkenes = 0
   numberamides = 0
   numbercarbamates = 0
   alkenecount = 0
   alkene_adjustment = 0

   for group in functionalgroups:
       if 'C=C' in group[2]:
          numberofalkenes = group[2].count('=')
          alkene_adjustment = numberofalkenes - 1
           # alkene_adjustment accounts for: for example, the loss in one count of functional group when a diene is counted as one type.

       if 'CNC(C)=O' in group[2] or 'cN(C)C(C)=O' in group[2]:
           numberamides +=1

       if 'COC(N)=O' in group[2]:
           numbercarbamates +=1

       functionalhetero = 0
       aromaticatoms = Chem.Mol.GetAromaticAtoms(molecule)
       aromaticatomindexes = [atom.GetIdx() for atom in aromaticatoms]
       for atom in molecule.GetAtoms():
           atom_index = atom.GetIdx
           inring = atom.IsInRing()
           atomname = atom.GetSymbol()
           if (inring == True and atomname != 'C' and atomname != 'H' and atom_index in aromaticatomindexes):
               functionalhetero += 1
   #print('hetero:', functionalhetero)
   print(aromaticatomindexes)
   functionalgroupindexes = [group.atomIds for group in functionalgroups] #if 'CN' not in group[2]]

   numberfunctionalgroups = len(functionalgroupindexes) + numberofalkenes + alkene_adjustment + numberamides + numbercarbamates - functionalhetero
   functionalgroupindexes = [item for sublist in functionalgroupindexes for item in sublist]

   print(functionalgroups)
   #print('numberoffunctionalgroups:' , numberfunctionalgroups)
   #print( 'numberofalkenes:', numberofalkenes)

   df.at[x,'numberfunctionalgroups'] = numberfunctionalgroups


   df.to_csv('correctedoutput.csv', sep =';',index = False)
