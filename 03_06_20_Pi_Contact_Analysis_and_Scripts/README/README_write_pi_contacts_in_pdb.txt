README written by Sara Y. Cheng on 03/06/20

Contacts_in_PDB.py originally written by Robert Vernon 

You need to following python scripts to run Contacts_in_PDB.py 

AnnotatePiPlanesBASE.py and AnnotateArcPiPlanes.py

See UsageNotes: /global/home/users/sycheng/software/scripts/PI_PI_SCRIPTS/Usage_Notes_Pi_Contact_scripts.txt

1. Run script to generate output of residues involved in pi-contacts for a single PDB

> python Contacts_in_PDB.py PDB.pdb | grep Contacts

Output format

SCSC: A.17:ARG*A.39:PHE   (sidechain:sidechain contact between residue 17 and residue 39, which are both chain A)
BBSC: A.1:SERASP*A.5:ARG (backbone:sidechain contact involving the peptide bond between residue 1 and residue 2 (which are serine and aspartate) and the sidechain of residue 5)
BBBB: A.1:SERASP*A.148:SERALA (backbone:backbone contact)

2. Run script to generated output of pi-contact types

> python Contacts_in_PDB.py PDB.pdb | grep Sums

Output format

Sums     ModelNum        0   #Atoms       78          5       16   100.00 #PiPi       16   #SCSC        5   #SCBB       11   #BBBB        0

The four numbers after #Atoms are the relevant ones:
1) Number of atoms in contacts
2) Number of native PDB sidechain-sidechain contacts present
3) Number of native PDB contacts present
4) Percentage of native PDB contacts present