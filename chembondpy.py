'''
ChemBondPy beta version 0.0
Dave Kuntz, UNT/MSU
https://github.com/davekuntz/ChemBondPy

Thank you for downloading ChemBondPy!
ChemBondPy allows you to determine whether or not two atoms are bonded to each other based on their identities and interatomic distances in angstroms.  It is a very simple Python 2.7 function that works by a chain of if statements.  Because of the variation in element electronegativity, size, etc., no firm rules exist for determining a bond's length.  This function takes several hundred atom combinations, including organic, inorganic, metalorganic compounds, sandwich complexes, biological compounds and more and empirically derived bond distances between a variety of elements in many chemical environments.  It can perfectly predict the bonding of roughly 99% of the 30,000 molecules in the NIST Chemistry Webbook.  It can also be easily modified, tweaked, and expanded.

METHODOLOGY:
The methodology for assigning bonds in chembonder is simple: start with the smallest element in a bond and work outward. Bond lengths are kept to a minimum to prevent falsely creating nonexistent bonds.  This allows exceptions to bonding patterns to be handled first, and is intuitive and easily modifiable.  For example:

Almost all hydrogen (H) bonds are under 1.2 angstroms.  However, exceptions can occur.  H-I compounds, for instance can bond lengths approaching 1.7 A.  Simply saying that all hydrogen-containing pairs of elements with a distance of less than 1.7 angstroms wouldn't work, because a lot of false positives would occur in other compounds, such as alkanes where a hydrogen could be mistaken for bonding to a carbon adjacent to the one it is truly bonded to.  Hence:

A small element is tested first
Exceptionally long bonds are tested first per that element
If no long bond is found but the element is present, then a generic bond length is used to determine the bonding
The generic bond length is always shorter or the same as the long bond exceptions, and represents the smallest distance safely used to estimate a bond.
If the element isn't present in the bond, move onto the next element
At the end, if there are no matches, a general distance of 1.6 angstroms or less is used to describe bonding
Certain organometallic compounds have both hydrogen and nonhydrogen ligands.  To ensure that the hydrogens aren't falsely bonded to the other ligands, the exceptions list is created.  See "Function Usage"


MAKING YOUR OWN MODIFICATIONS:
You are welcome to make any modifications to this script.  Some of the bonds currently defined in it are only defined at a very generic level.  If you'd like to share your bond length modifications, please send me your changes to the email address listed under "Getting Help.'  You may also request to update your modifications through GitHub directly.  Please keep this documentation with the script.

INSTALLATION:
To install, save chembonder.py into your working directory or copy the function directly into your python file.  Either way, please keep this documentation intact.  To call the bonder function from its own file, type in
	from chembondpy import bonder

FUNCTION USEAGE:
bonder takes three arguments, bonder(pair,l,unique_atoms):
    pair: pair is a list of two elements by element abbreviation, IE ['C','H']
    l: l is the length in question between the elements in the list
    unique_atoms: unique_atoms is a list containing exactly one of each element present in the molecule of interest.  For instance, water's unique list would be ['H','O'] because elements aren't repeated twice (although it won't harm the function).  unique_atoms is used because certain organometallic compounds have bond lengths that are significantly different due to the presence of certain transition metals.  This option allows custom bond lengths when these elements are or are not present in a molecule but not directly present in a bond of interest.  This argument can also be ignored by simply including an empty list.
       
Examples:
#Bond Distance is 0.9 A, which is a bond
>>>> print bonder(['H','O'],0.9,[])
1

#Bond distance is 2.9 A, too far for a H-O bond
>>>> print bonder(['H','O'],2.9,[]) 
0

#Bond distance is 1.3 A, which is in range for a H-Cl bond
>>>> print bonder(['H','Cl'],1.3,[]) 
1

#Bond distance is 1.3 A, which is in range for a H-Cl bond, BUT nickel is present in the unique atom list, so the hydrogen is not bonding to the chloride! (The molecule is nickel dichloride in this case)
>>>> print bonder(['H','Cl'],1.3,['Ni','H','Cl']) 
0

OUTPUT:
    The output is very simple.  The function returns 1 when there is a bond present and 0 when the two atoms of interest are not connected.
    
COMING SOON:
Future planned functionality includes:
    Converting xyz matrices to distance matrices
    Deriving adjacency matrices using the bond distance
    Counting bonds and bond types
    Automatically generating descriptors of molecules from xyz distances. MOL files

KNOWN ISSUES:
Certain bicyclo organic compounds may have very long carbon-carbon bonds which are ignored.  Most ionic bonds are still counted as bonds per my preference, but this can be changed by setting result to 0 for any bond in question.
C-N bonds are difficult and are set to 1.6 A.  Many ammine bonds are slightly longer, but increasing the length adds false positives.  Adjust this bond length to make more or fewer C-N bonds based on the type of system you are working with.

GETTING HELP:

Feel free to contact me at davidmicahkuntz ____at gmail.com.  I assume you know how to input the @ sign in its proper place.  This is a free program for both commercial and personal use, and I assume no responsibility for how you use it or any of its results.

CITE AS:
Kuntz, David M.  ChemBondPy V. 0.0, bonder function.  2016


'''


#Determine if two atoms bond from distance
def bonder(pair,l,unique_atoms):
    exceptions = ['Be', 'Na', 'Mg', 'Al', 'K', 'Ti', 'V', 'Cr', 'Fe', 'Ni', 'Ge', 'Zr', 'Nb', 'Mo', 'In', 'Sn', 'Sb', 'Cs', 'Ba', 'La', 'Gd', 'Lu', 'Ta', 'W', 'Ir', 'Hg', 'Pb', 'Bi']
    any_in = lambda a, b: any(i in b for i in a) #Function for determining if any atoms in molecule are in exceptions list
    
    result = -1 #Result not found
    
    if l == 0: result = 0 #Zero length between = non bonding (generally the same atom refering to itself)
    
    #Bonds with hydrogen
    elif 'H' in pair:
        if pair == ['H','H'] and l < 0.8 and any_in(unique_atoms,exceptions) == False: result = 1 #H-H diatomic bonding
        elif 'He' in pair and l < 0.8: result = 1 #H-He
        elif 'Li' in pair and l < 1.7: result = 1 #H-Li
        elif 'Be' in pair and l < 1.4: result = 1 #H-Be
        elif 'B' in pair and l < 1.5: result = 1 #H-B
        elif 'C' in pair and l < 1.2: result = 1 #H-C in Al or Hg complex
        elif 'N' in pair and l < 1.5: result = 1 #H-N
        elif 'O' in pair and l < 1.1 and (any_in(unique_atoms,exceptions) == False or 'Na' in unique_atoms or 'K' in unique_atoms or 'Ba' in unique_atoms): result = 1
        elif 'F' in pair and l < 1.1 and any_in(unique_atoms,exceptions) == False: result = 1 #H-F
        elif 'Mg' in pair and l < 1.8: result = 1 #H-Mg
        elif 'Al' in pair and l < 1.8: result = 1 #H-Al This one is questionable
        elif 'Si' in pair and l < 1.6: result = 1 #H-Si
        elif 'P' in pair and l < 1.5: result = 1 #H-P
        elif 'S' in pair and l < 1.5: result = 1 #H-S
        elif 'Cl' in pair and l < 1.45 and any_in(unique_atoms,exceptions) == False: result = 1 #H-Cl
        elif 'Ca' in pair and l < 2.2: result = 1 #H-Ca
        elif 'Sc' in pair and l < 1.9: result = 1 #H-Sc
        elif 'Ti' in pair and l < 1.8: result = 1 #H-Ti
        elif 'Mn' in pair and l < 1.6: result = 1 #H-Mn
        elif 'Co' in pair and l < 1.5: result = 1 #H-Co
        elif 'Ni' in pair and l < 1.5: result = 1 #H-Ni
        elif 'Ga' in pair and l < 1.7: result = 1 #H-Ga
        elif 'Ge' in pair and l < 1.6: result = 1 #H-Ge
        elif 'As' in pair and l < 1.6: result = 1 #H-As
        elif 'Se' in pair and l < 1.6: result = 1 #H-Se
        elif 'Br' in pair and l < 1.6 and any_in(unique_atoms,exceptions) == False: result = 1 #H-Br
        elif 'Te' in pair and l < 1.7: result = 1 #H-Te
        elif 'I' in pair and l < 1.7 and any_in(unique_atoms,exceptions) == False: result = 1 #H-I
        elif 'La' in pair and l < 2.05: result = 1 #H-La
        #elif l < 1.2 and any_in(unique_atoms,exceptions) == False: result = 1 #Generic hydrogen bonding length
        else: result = 0 #Results with long H distances are nonbonding
    
    #Bonds with helium
    elif 'He' in pair:
        if pair == ['He', 'He'] and l < 2.8: result = 1 #He-He diatomic bonding
        
    #Bonds with lithium
    elif 'Li' in pair:
        if pair == ['Li', 'Li'] and l < 2.2: result = 1 #Li-Li diatomic bonding
        elif pair == ['Li', 'Li'] and l < 2.8 and len(unique_atoms) == 1: result = 1 #Diatomic Li-Li has longer bonds
        elif 'B' in pair: result = 0 #B-Li
        elif 'C' in pair and l < 2.1: result = 1 #Li-C
        elif 'O' in pair and l < 1.7: result = 1 #Li-O
        elif 'F' in pair and l < 1.8: result = 1 #Li-F
        elif 'Si' in pair: result = 0 #Li-Si
        elif 'S' in pair: result = 0 #Li-S
        elif 'Br' in pair and l < 2.2: result = 1 #Li-Br
        elif l < 2.1: result = 1 #Generic lithium bonding template
    
    
    #Bonds with beryllium
    elif 'Be' in pair:
        if pair == ['Be', 'Be'] and l < 2.1: result = 1 #Be-Be diatomic bonding
        elif 'C' in pair and l < 1.7: result = 1 #Be-C
        elif 'O' in pair and l < 1.7: result = 1 #Be-O
        elif 'S' in pair and l < 2.0: result = 1 #Be-S
        elif 'Cl' in pair and l < 3.9: result = 1 #Be-Cl
        elif 'Br' in pair and l < 3.9: result = 1 #Be-Br
        elif 'I' in pair and l < 6.1: result = 1 #Be-I
        elif l < 1.7: result = 1 #Generic beryllium bonding length
    
    #Bonds with boron    
    elif 'B' in pair:
        if pair == ['B', 'B'] and l < 1.8: result = 1 #B-B diatomic bonding
        elif 'C' in pair and l < 1.7: result = 1 #B-C
        elif 'N' in pair and l < 1.7: result = 1 #B-N
        elif 'O' in pair and l < 1.7: result = 1 #B-O
        elif 'F' in pair and l < 1.6: result = 1 #B-F
        elif 'Si' in pair and l < 2.1: result = 1 #B-Si
        elif 'P' in pair and l < 1.9: result = 1 #B-P
        elif 'S' in pair and l < 2.0: result = 1 #B-S
        elif 'Cl' in pair and l < 2.0: result = 1 #B-Cl
        elif 'Se' in pair and l < 1.9: result = 1 #B-Se
        elif 'Br' in pair and l < 2.0: result = 1 #B-Br
        elif 'I' in pair and l < 2.3: result = 1 #B-I
        elif l < 1.7: result = 1 #Generic boron bonding length
        
    #Bonds with carbon    
    elif 'C' in pair:
        if pair == ['C', 'C'] and l < 2.0: result = 1 #C-C diatomic bonding
        elif 'N' in pair and l < 1.6: result = 1 #C-N 1.55
        elif 'O' in pair and l < 1.6: result = 1 #C-O
        elif 'F' in pair and l < 1.5: result = 1 #C-F
        elif 'Al' in pair and l < 2.2: result = 1 #C-Al
        elif 'Si' in pair and l < 2.0: result = 1 #C-Si
        elif 'P' in pair and l < 2.0: result = 1 #C-P
        elif 'S' in pair and l < 1.93: result = 1 #C-S
        elif 'Cl' in pair and l < 2.0: result = 1 #C-Cl
        elif 'Ti' in pair and l < 2.2: result = 1 #C-Ti
        elif 'V' in pair and l < 2.1: result = 1 #C-V
        elif 'Cr' in pair and l < 2.8: result = 1 #C-Cr
        elif 'Mn' in pair and l < 2.2: result = 1 #C-Mn
        elif 'Fe' in pair and l < 2.2: result = 1 #C-Fe
        elif 'Co' in pair and l < 2.2: result = 1 #C-Co
        elif 'Ni' in pair and l < 2.8: result = 1 #C-Ni
        elif 'Cu' in pair and l < 1.8: result = 1 #C-Cu
        elif 'Zn' in pair and l < 2.0: result = 1 #C-Zn
        elif 'Ge' in pair and l < 2.1: result = 1 #C-Ge
        elif 'As' in pair and l < 2.1: result = 1 #C-As
        elif 'Se' in pair and l < 2.1: result = 1 #C-Se
        elif 'Br' in pair and l < 2.3: result = 1 #C-Br
        elif 'Mo' in pair and l < 4.3: result = 1 #C-Mo
        elif 'Ru' in pair and l < 2.3: result = 1 #C-Ru
        elif 'Rh' in pair and l < 2.4: result = 1 #C-Rh
        elif 'Pd' in pair and l < 2.3: result = 1 #C-Pd
        elif 'In' in pair and l < 2.2: result = 1 #C-In
        elif 'Sn' in pair and l < 2.3: result = 1 #C-Sn
        elif 'Sb' in pair and l < 2.3: result = 1 #C-Sb
        elif 'Te' in pair and l < 2.3: result = 1 #C-Te
        elif 'I' in pair and l < 2.4: result = 1 #C-I
        elif 'Hf' in pair and l < 2.6: result = 1 #C-Hf
        elif 'W' in pair and l < 2.3: result = 1 #C-W
        elif 'Re' in pair and l < 2.4: result = 1 #C-Re
        elif 'Au' in pair and l < 2.1: result = 1 #C-Au
        elif 'Hg' in pair and l < 3.2: result = 1 #C-Hg
        elif 'Tl' in pair and l < 2.8: result = 1 #C-Tl
        elif 'Pb' in pair and l < 2.4: result = 1 #C-Pb
        elif 'Bi' in pair and l < 2.4: result = 1 #C-Bi
        elif l < 1.6: result = 1 #Generic carbon bonding length

    #Bonds with nitrogen
    elif 'N' in pair:
        if pair == ['N', 'N'] and l < 1.7: result = 1 #N-N diatomic bonding
        elif 'O' in pair and l < 1.6: result = 1 #N-O
        elif 'F' in pair and l < 1.5: result = 1 #N-F
        elif 'Al' in pair and l < 1.7: result = 1 #N-Al
        elif 'Si' in pair and l < 1.9: result = 1 #N-Si
        elif 'P' in pair and l < 1.9: result = 1 #N-P
        elif 'S' in pair and l < 1.9: result = 1 #N-S
        elif 'Cl' in pair and l < 2.1: result = 1 #N-Cl
        elif 'Ti' in pair and l < 2.0: result = 1 #N-Ti
        elif 'V' in pair and l < 1.9: result = 1 #N-V
        elif 'Cr' in pair and l < 1.6: result = 1 #N-Cr
        elif 'Mn' in pair and l < 2.2: result = 1 #N-Mn
        elif 'Fe' in pair and l < 2.1: result = 1 #N-Fe
        elif 'Co' in pair and l < 2.0: result = 1 #N-Co
        elif 'Ni' in pair and l < 2.0: result = 1 #N-Ni
        elif 'Cu' in pair and l < 2.0: result = 1 #N-Cu
        elif 'Zn' in pair and l < 2.0: result = 1 #N-Zn
        elif 'Ge' in pair and l < 1.9: result = 1 #N-Ge
        elif 'Se' in pair and l < 2.0: result = 1 #N-Se
        elif 'Br' in pair and l < 2.2: result = 1 #N-Br
        elif 'Nb' in pair and l < 1.7: result = 1 #N-Nb
        elif 'Pt' in pair and l < 2.2: result = 1 #N-Pt
        elif 'Sn' in pair and l < 2.1: result = 1 #N-Sn
        elif 'I' in pair and l < 3.05: result = 1 #N-I
        elif 'Hf' in pair and l < 18: result = 1 #N-Hf
        elif l < 1.5: result = 1 #Generic nitrogen bonding length
        
    #Bonds with oxygen
    elif 'O' in pair:
        if pair == ['O', 'O'] and l < 1.5 and any_in(unique_atoms,exceptions) == False: result = 1 #O-O diatomic bonding
        elif 'F' in pair and l < 1.5: result = 1 #O-F
        elif 'Na' in pair and l < 2.1: result = 1 #O-Na
        elif 'Mg' in pair and l < 2.0: result = 1 #O-Mg
        elif 'Al' in pair and l < 2.1: result = 1 #O-Al
        elif 'Si' in pair and l < 1.8: result = 1 #O-Si
        elif 'P' in pair and l < 1.8: result = 1 #O-P
        elif 'S' in pair and l < 1.8: result = 1 #O-S
        elif 'Cl' in pair and l < 1.8: result = 1 #O-Cl
        elif 'K' in pair and l < 2.4: result = 1 #O-K
        elif 'Ca' in pair and l < 2.4: result = 1 #O-Ca
        elif 'Ti' in pair and l < 2.1: result = 1 #O-Ti
        elif 'V' in pair and l < 1.8: result = 1 #O-V
        elif 'Cr' in pair and l < 2.0: result = 1 #O-Cr
        elif 'Mn' in pair and l < 1.6: result = 1 #O-Mn
        elif 'Fe' in pair and l < 2.0: result = 1 #O-Fe
        elif 'Co' in pair and l < 1.9: result = 1 #O-Co
        elif 'Cu' in pair and l < 2.0: result = 1 #O-Cu
        elif 'Zn' in pair and l < 2.0: result = 1 #O-Zn
        elif 'Ge' in pair and l < 1.9: result = 1 #O-Ge
        elif 'As' in pair and l < 1.9: result = 1 #O-As
        elif 'Se' in pair and l < 1.7: result = 1 #O-Se
        elif 'Br' in pair and l < 1.9: result = 1 #O-Br
        elif 'Sr' in pair and l < 2.2: result = 1 #O-Sr
        elif 'Y' in pair and l < 2.1: result = 1 #O-Y
        elif 'Zr' in pair and l < 1.8: result = 1 #O-Zr
        elif 'Nb' in pair and l < 1.7: result = 1 #O-Nb
        elif 'Mo' in pair and l < 2.0: result = 1 #O-Mo
        elif 'Rh' in pair and l < 2.4: result = 1 #O-Rh
        elif 'Pd' in pair and l < 1.8: result = 1 #O-Pd
        elif 'Cd' in pair and l < 2.0: result = 1 #O-Cd
        elif 'In' in pair and l < 2.2: result = 1 #O-In
        elif 'Sn' in pair and l < 2.2: result = 1 #O-Sn
        elif 'Sb' in pair and l < 2.0: result = 1 #O-Sb
        elif 'Te' in pair and l < 1.9: result = 1 #O-Te
        elif 'I' in pair and l < 2.1: result = 1 #O-I
        elif 'Ba' in pair and l < 2.4: result = 1 #Ba-V
        elif 'W' in pair and l < 3.2: result = 1 #O-W
        elif 'Re' in pair and l < 2.0: result = 1 #O-Re
        elif 'Ir' in pair and l < 1.8: result = 1 #O-Ir
        elif 'Pt' in pair and l < 1.8: result = 1 #O-Pt
        elif 'Hg' in pair and l < 2.0: result = 1 #O-Hg
        elif 'Tl' in pair and l < 2.5: result = 1 #O-Tl
        elif l < 1.8: result = 1 #Generic oxygen bonding length

    #Bonds with flourine
    elif 'F' in pair:
        if pair == ['F', 'F'] and l < 1.5: result = 1 #F-F diatomic bonding
        elif 'Na' in pair and l < 2.3: result = 1 #F-Na
        elif 'Mg' in pair and l < 3.2: result = 1 #F-Mg
        elif 'Al' in pair and l < 1.8: result = 1 #F-Al
        elif 'Si' in pair and l < 1.7: result = 1 #F-Si
        elif 'P' in pair and l < 1.7: result = 1 #F-P
        elif 'S' in pair and l < 1.8: result = 1 #F-S
        elif 'Cl' in pair and l < 1.7: result = 1 #F-Cl
        elif 'K' in pair and l < 2.2: result = 1 #F-K
        elif 'Ca' in pair and l < 2.1: result = 1 #F-Ca
        elif 'Sc' in pair and l < 1.9: result = 1 #F-Sc
        elif 'Ti' in pair and l < 1.8: result = 1 #F-Ti
        elif 'V' in pair and l < 1.8: result = 1 #F-V
        elif 'Cr' in pair and l < 2.4: result = 1 #F-Cr
        elif 'Mn' in pair and l < 1.7: result = 1 #F-Mn
        elif 'Fe' in pair and l < 1.9: result = 1 #F-Fe
        elif 'Cu' in pair and l < 1.7: result = 1 #F-Cu
        elif 'Ge' in pair and l < 3.8: result = 1 #F-Ge
        elif 'As' in pair and l < 1.9: result = 1 #F-As
        elif 'Se' in pair and l < 1.8: result = 1 #F-Se
        elif 'Br' in pair and l < 1.9: result = 1 #F-Br
        elif 'Sr' in pair and l < 2.1: result = 1 #F-Sr
        elif 'Mo' in pair and l < 3.9: result = 1 #F-Mo
        elif 'Sn' in pair and l < 3.3: result = 1 #F-Sn
        elif 'Sb' in pair and l < 1.9: result = 1 #F-Sb
        elif 'I' in pair and l < 2.0: result = 1 #F-I
        elif 'Cs' in pair and l < 3.0: result = 1 #F-Cs
        elif 'Xe' in pair and l < 3.6: result = 1 #F-Xe
        elif 'Ba' in pair and l < 2.2: result = 1 #F-Ba
        elif 'La' in pair and l < 2.6: result = 1 #F-La
        elif 'Gd' in pair and l < 1.8: result = 1 #F-Gd
        elif 'W' in pair and l < 3.4: result = 1 #F-W
        elif 'Re' in pair and l < 1.9: result = 1 #F-Re
        elif 'Ir' in pair and l < 3.7: result = 1 #F-Ir
        elif 'Hg' in pair and l < 3.3: result = 1 #F-Hg
        elif 'Tl' in pair and l < 2.0: result = 1 #F-Tl
        elif 'Pb' in pair and l < 3.3: result = 1 #F-Pb
        elif 'Bi' in pair and l < 3.9: result = 1 #F-Bi
        elif l < 1.7: result = 1 #Generic flouride bonding length
    
    #Bonds with sodium
    elif 'Na' in pair:
        if pair == ['Na', 'Na'] and l < 3.1: result = 1 #Na-Na diatomic bonding
        elif 'P' in pair and l < 2.7: result = 1 #Na-P
        elif 'Cl' in pair and l < 4.1: result = 1 #Na-P
        elif l < 2.7: result = 1 #Generic sodium bonding length
        
    #Bonds with magnesium
    elif 'Mg' in pair:
        if pair == ['Mg', 'Mg'] and l < 7.4: result = 1 #Mg-Mg diatomic bonding
        elif 'Mg' in pair and l < 4.2: result = 1 #Mg-Cl
        elif 'Br' in pair and l < 3.3: result = 1 #Mg-Br
        elif 'I' in pair and l < 6.2: result = 1 #Mg-I
        elif l < 3.2: result = 1 #Generic magnesium bonding length
    
    #Bonds with aluminum
    elif 'Al' in pair:
        if pair == ['Al', 'Al'] and l < 2.6: result = 1 #Al-Al diatomic bonding
        elif 'P' in pair and l < 2.1: result = 1 #Al-P
        elif 'S' in pair and l < 2.1: result = 1 #Al-S
        elif 'Cl' in pair and l < 2.8: result = 1 #Al-Cl
        elif 'Br' in pair and l < 2.3: result = 1 #Al-Br
        elif 'I' in pair and l < 5.1: result = 1 #Al-I
        elif l < 2.0: result = 1 #Generic aluminum bonding length
        
    #Bonds with silicon
    elif 'Si' in pair:
        if pair == ['Si', 'Si'] and l < 2.6: result = 1 #Si-Si diatomic bonding
        elif 'P' in pair and l < 2.4: result = 1 #Si-P
        elif 'S' in pair and l < 2.2: result = 1 #Si-S
        elif 'Cl' in pair and l < 2.2: result = 1 #Si-Cl
        elif 'Co' in pair and l < 2.4: result = 1 #Si-Co
        elif 'Ge' in pair and l < 2.4: result = 1 #Si-Ge
        elif 'As' in pair and l < 2.4: result = 1 #Si-As
        elif 'Se' in pair and l < 2.3: result = 1 #Si-Se
        elif 'Br' in pair and l < 2.3: result = 1 #Si-Br
        elif 'Sn' in pair and l < 2.7: result = 1 #Si-Sn
        elif 'Te' in pair and l < 2.6: result = 1 #Si-Te
        elif 'I' in pair and l < 2.6: result = 1 #Si-I
        elif 'Hg' in pair and l < 2.6: result = 1 #Si-Hg
        elif l < 2.2: result = 1 #Generic silicon bonding length
    
    #Bonds with phosphorous
    elif 'P' in pair: 
        if pair == ['P', 'P'] and l < 2.3: result = 1 #P-P diatomic bonding
        elif 'S' in pair and l < 2.4: result = 1 #P-S
        elif 'Cl' in pair and l < 2.3: result = 1 #P-Cl
        elif 'Fe' in pair and l < 2.2: result = 1 #P-Fe
        elif 'Ge' in pair and l < 2.4: result = 1 #P-Ge
        elif 'As' in pair and l < 2.4: result = 1 #P-As
        elif 'Se' in pair and l < 2.3: result = 1 #P-Se
        elif 'Br' in pair and l < 2.3: result = 1 #P-Br
        elif 'In' in pair and l < 2.5: result = 1 #P-In
        elif 'Te' in pair and l < 2.6: result = 1 #P-Te
        elif 'I' in pair and l < 2.6: result = 1 #P-I
        elif l < 2.2: result = 1 #Generic phosphorous bonding length
        
    #Bonds with sulfur
    elif 'S' in pair:
        if pair == ['S', 'S'] and l < 2.5: result = 1 #S-S diatomic bonding
        elif 'Cl' in pair and l < 2.3: result = 1 #S-Cl
        elif 'Fe' in pair and l < 2.0: result = 1 #S-Fe
        elif 'Ni' in pair and l < 2.1: result = 1 #S-Ni
        elif 'Cu' in pair and l < 2.0: result = 1 #S-Cu
        elif 'Zn' in pair and l < 2.4: result = 1 #S-Zn
        elif 'Ge' in pair and l < 2.3: result = 1 #S-Ge
        elif 'As' in pair and l < 2.3: result = 1 #S-As
        elif 'Se' in pair and l < 2.1: result = 1 #S-Se
        elif 'Br' in pair and l < 2.4: result = 1 #S-Br
        elif 'Mo' in pair and l < 2.4: result = 1 #S-Mo
        elif 'Sn' in pair and l < 2.6: result = 1 #S-Sn
        elif 'Sb' in pair and l < 2.6: result = 1 #S-Sb
        elif 'I' in pair and l < 2.6: result = 1 #S-I
        elif 'Pt' in pair and l < 2.2: result = 1 #S-Pt
        elif 'Pb' in pair and l < 2.6: result = 1 #S-Pb
        elif l < 2.0: result = 1 #Generic sulfur bonding length
        
    #Bonds with chlorine
    elif 'Cl' in pair:
        if pair == ['Cl', 'Cl'] and l < 2.1: result = 1 #Cl-Cl diatomic bonding
        elif 'K' in pair and l < 3.25: result = 1 #Cl-K
        elif 'Ca' in pair and l < 2.6: result = 1 #Cl-Ca
        elif 'V' in pair and l < 2.2: result = 1 #Cl-V
        elif 'Cr' in pair and l < 2.2: result = 1 #Cl-Cr
        elif 'Fe' in pair and l < 2.3: result = 1 #Cl-Fe
        elif 'Co' in pair and l < 2.3: result = 1 #Cl-Co
        elif 'Ni' in pair and l < 2.1: result = 1 #Cl-Ni
        elif 'Cu' in pair and l < 2.1: result = 1 #Cl-Cu
        elif 'Zn' in pair and l < 2.1: result = 1 #Cl-Zn
        elif 'Ge' in pair and l < 2.2: result = 1 #Cl-Ge
        elif 'As' in pair and l < 2.3: result = 1 #Cl-As
        elif 'Se' in pair and l < 2.4: result = 1 #Cl-Se
        elif 'Br' in pair and l < 2.2: result = 1 #Cl-Br
        elif 'Rb' in pair and l < 3.6: result = 1 #Cl-Rb
        elif 'Sr' in pair and l < 2.7: result = 1 #Cl-Sr
        elif 'Nb' in pair and l < 4.9: result = 1 #Cl-Nb
        elif 'Ti' in pair and l < 5.2: result = 1 #Cl-Ti
        elif 'Mo' in pair and l < 4.9: result = 1 #Cl-Mo
        elif 'Ag' in pair and l < 2.4: result = 1 #Cl-Ag
        elif 'In' in pair and l < 3.8: result = 1 #Cl-In
        elif 'Sn' in pair and l < 4.3: result = 1 #Cl-Sn
        elif 'Sb' in pair and l < 7.9: result = 1 #Cl-Sb
        elif 'Te' in pair and l < 2.6: result = 1 #Cl-Te
        elif 'Cs' in pair and l < 4.2: result = 1 #Cl-Cs
        elif 'Ba' in pair and l < 2.8: result = 1 #Cl-Ba
        elif 'La' in pair and l < 5.7: result = 1 #Cl-La
        elif 'Gd' in pair and l < 3.2: result = 1 #Cl-Gd
        elif 'Lu' in pair and l < 3.0: result = 1 #Cl-Lu
        elif 'Hf' in pair and l < 2.5: result = 1 #Cl-Hf
        elif 'Ta' in pair and l < 4.9: result = 1 #Cl-Ta
        elif 'W' in pair and l < 7.2: result = 1 #Cl-W
        elif 'Re' in pair and l < 2.6: result = 1 #Cl-Re
        elif 'Pt' in pair and l < 2.5: result = 1 #Cl-Pt
        elif 'Hg' in pair and l < 8.2: result = 1 #Cl-Hg
        elif 'Pb' in pair and l < 3.2: result = 1 #Cl-Pb
        elif 'Bi' in pair and l < 4.5: result = 1 #Cl-Bi
        elif 'Th' in pair and l < 2.6: result = 1 #Cl-Th
        elif l < 2.1: result = 1 #Generic chlorine bonding length
    
    #Bonds with potassium
    elif 'K' in pair:
        if pair == ['K', 'K'] and l < 4.0: result = 1 #K-K diatomic bonding
        elif 'Br' in pair and l < 3.4: result = 1 #K-Br
        elif 'I' in pair and l < 3.6: result = 1 #K-I
        if l < 3.3: result = 1 #Generic potassium bonding length
    
    
    #Bonds with scandium
    elif 'Sc' in pair:
        if 'Se' in pair and l < 2.5: result = 1 #Sc-Se
       
    #Bonds with manganese
    elif 'Mn' in pair:
        if 'As' in pair and l < 2.5: result = 1 #Mn-As
        elif 'Se' in pair and l < 2.1: result = 1 #Mn-Se
        elif 'I' in pair and l < 2.9: result = 1 #Mn-I
        elif 'Re' in pair and l < 3.1: result = 1 #Mn-Re
        elif l < 2.1: result = 1 #Generic manganese bonding length
        
    #Bonds with iron
    elif 'Fe' in pair:
        if 'Br' in pair and l < 2.4: result = 1 #Fe-Br
        elif 'I' in pair and l < 2.8: result = 1 #Fe-I
        elif l < 2.4: result = 1 #Generic iron bonding length    
    
    #Bonds with cobalt
    elif 'Co' in pair:
        if 'Zn' in pair and l < 2.4: result = 1 #Co-Zn
        if 'Cd' in pair and l < 2.6: result = 1 #Co-Cd
        elif l < 2.4: result = 1 #Generic cobalt bonding length

    #Bonds with copper
    elif 'Cu' in pair:
        if pair == ['Cu', 'Cu'] and l < 2.1: result = 1 #Cu-Cu diatomic bonding
    
    #Bonds with gallium
    elif 'Ga' in pair:
        if 'Sb' in pair and l < 2.7: result = 1 #Ga-Sb
        elif 'Te' in pair and l < 1.8: result = 1 #Ga-Te
        elif l < 1.7: result = 1 #Generic gallium bonding length

    #Bonds with germanium
    elif 'Ge' in pair:
        if pair == ['Ge', 'Ge'] and l < 2.5: result = 1 #Ge-Ge diatomic bonding
        elif 'Se' in pair and l < 3.8: result = 1 #Ge-Se
        elif 'Br' in pair and l < 2.4: result = 1 #Ge-Br
        elif 'Sn' in pair and l < 2.7: result = 1 #Ge-Sn
        elif 'Te' in pair and l < 2.6: result = 1 #Ge-Te
        elif 'I' in pair and l < 2.7: result = 1 #Ge-I
        elif l < 2.4: result = 1 #Generic germanium bonding length
    
    #Bonds with arsenic
    elif 'As' in pair:
        if pair == ['As', 'As'] and l < 2.5: result = 1 #As-As diatomic bonding
        elif 'Se' in pair and l < 2.4: result = 1 #As-Se
        elif 'Br' in pair and l < 2.4: result = 1 #As-Br
        elif 'I' in pair and l < 2.7: result = 1 #As-I
        elif l < 2.4: result = 1 #Generic arsenic bonding length
        
    #Bonds with selenium
    elif 'Se' in pair:
        if pair == ['Se', 'Se'] and l < 2.4: result = 1 #Se-Se diatomic bonding
        elif 'Br' in pair and l < 2.4: result = 1 #Se-Br
        elif 'Y' in pair and l < 2.7: result = 1 #Se-Y
        elif 'In' in pair and l < 2.9: result = 1 #Se-In
        elif 'Sn' in pair and l < 2.6: result = 1 #Se-Sn
        elif 'Te' in pair and l < 2.4: result = 1 #Se-Te
        elif 'I' in pair and l < 2.7: result = 1 #Se-I
        elif 'Pb' in pair and l < 2.8: result = 1 #Se-Pb
        elif l < 2.4: result = 1 #Generic selenium bonding length   
    
    #Bonds with bromine
    elif 'Br' in pair:
        if pair == ['Br', 'Br'] and l < 2.4: result = 1 #Br-Br diatomic bonding
        elif 'Zr' in pair and l < 3.3: result = 1 #Br-Zr
        elif 'Mo' in pair and l < 2.6: result = 1 #Br-Mo
        elif 'In' in pair and l < 3.2: result = 1 #Br-In
        elif 'Sn' in pair and l < 3.2: result = 1 #Br-Sn
        elif 'Sb' in pair and l < 2.7: result = 1 #Br-Sb
        elif 'Te' in pair and l < 3.1: result = 1 #Br-Te
        elif 'Cs' in pair and l < 3.3: result = 1 #Br-Cs
        elif 'Ba' in pair and l < 6.2: result = 1 #Br-Ba
        elif 'W' in pair and l < 3.5: result = 1 #Br-W
        elif 'Hg' in pair and l < 3.9: result = 1 #Br-Hg
        elif 'Pb' in pair and l < 6.1: result = 1 #Br-Pb
        elif 'Bi' in pair and l < 4.1: result = 1 #Br-Bi
        elif l < 3.2: result = 1 #Generic bromine bonding length
    
    #Bonds with rubidium
    elif 'Kr' in pair:
        if pair == ['Kr', 'Kr'] and l < 3.9: result = 1 #Kr-Kr diatomic bonding
        
    #Bonds with rubidium
    elif 'Rb' in pair:
        if pair == ['Rb', 'Rb'] and l < 4.3: result = 1 #Rb-Rb diatomic bonding
    
    #Bonds with strontium
    #elif 'Sr' in pair:
        #if l < 2.2: result = 1 #Generic strontium bonding length  
            
    #Bonds with zirconium
    elif 'Zr' in pair:
        if 'I' in pair and l < 3.5: result = 1 #Zr-I 
               
    #Bonds with molybdenum
    elif 'Mo' in pair:
        if pair == ['Mo', 'Mo'] and l < 2.2: result = 1 #Mo-Mo diatomic bonding
        elif 'I' in pair and l < 6.2: result = 1 #Mo-I
        #if l < 2.1: result = 1 #Generic molybdenum bonding length 
        
    #Bonds with silver
    elif 'Ag' in pair:
        if pair == ['Ag', 'Ag'] and l < 2.6: result = 1 #Ag-Ag diatomic bonding
    
    #Bonds with indium
    elif 'In' in pair:
        if 'Te' in pair and l < 3.4: result = 1 #In-Te
        elif 'I' in pair and l < 6.5: result = 1 #In-I
        
    #Bonds with tin
    elif 'Sn' in pair:
        if pair == ['Sn', 'Sn'] and l < 2.9: result = 1 #Sn-Sn diatomic bonding
        if 'Te' in pair and l < 2.8: result = 1 #Sn-Te
        elif 'I' in pair and l < 3.6: result = 1 #Sn-I
        elif l < 2.8: result = 1 #Generic tin bonding length
        
    #Bonds with antimony
    elif 'Sb' in pair:
        if pair == ['Sb', 'Sb'] and l < 2.6: result = 1 #Sb-Sb diatomic bonding
        elif 'I' in pair and l < 2.9: result = 1 #Sb-I
        elif l < 2.6: result = 1 #Generic antimony bonding length
    
    #Bonds with tellurium
    elif 'Te' in pair:
        if pair == ['Te', 'Te'] and l < 2.9: result = 1 #Te-Te diatomic bonding
        elif 'I' in pair and l < 3.5: result = 1 #Te-I
        
    #Bonds with iodine
    elif 'I' in pair:
        if pair == ['I', 'I'] and l < 3.0: result = 1 #I-I diatomic bonding
        elif 'Cs' in pair and l < 4.4: result = 1 #I-Cs
        elif 'Ba' in pair and l < 3.3: result = 1 #I-Ba
        elif 'Lu' in pair and l < 10.6: result = 1 #I-Lu
        elif 'Re' in pair and l < 2.9: result = 1 #I-Re
        elif 'Hg' in pair and l < 4.2: result = 1 #I-Hg
        elif 'Pb' in pair and l < 3.9: result = 1 #I-Pb
        elif 'Bi' in pair and l < 6.4: result = 1 #I-Bi
        elif l < 3.3: result = 1 #Generic iodine bonding length    
    
    #Bonds with Xenon
    elif 'Xe' in pair:
        if pair == ['Xe', 'Xe'] and l < 7.1: result = 1 #Xe-Xe diatomic bonding
        
    #Bonds with cesium
    elif 'Cs' in pair:
        if pair == ['Cs', 'Cs'] and l < 4.8: result = 1 #Cs-Cs diatomic bonding
                
    #Bonds with rhenium
    elif 'Re' in pair:
        if pair == ['Re', 'Re'] and l < 2.6: result = 1 #Re-Re diatomic bonding    
  
    #Generic Bond for anything not covered, 1.6 A or less by default
    if result == -1:
        if l < 1.4 and 'H' not in pair: result = 1
        else: result = 0
    return result
