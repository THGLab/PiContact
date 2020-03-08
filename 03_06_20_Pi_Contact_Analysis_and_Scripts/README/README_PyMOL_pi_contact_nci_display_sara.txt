README written by Sara Y. Cheng on 03/06/20

MakePi.py, MakePiPlanes.py originally written by Robert Vernon

nci.py originally written by Julia Contreras-Garcia’s group

https://github.com/rmera/ncipy

1. Download PyMOL and Install on your OS

Educational PyMOL for Mac https://pymol.org/edu/?q=educational/

Opensource PyMOL for Mac https://pymolwiki.org/index.php/MAC_Install

2. Create ~/.pymolrc in home directory

add the following lines to this file and save

“import nci”
import MakePi”

3. Open your PyMOL installation and execute the following line on the 

“print sys.path”

PyMOL>print sys.path ['/', '/Applications/PyMOL.app/Contents/lib/python37.zip', '/Applications/PyMOL.app/Contents/lib/python3.7', '/Applications/PyMOL.app/Contents/lib/python3.7/lib-dynload', '/Applications/PyMOL.app/Contents/lib/python3.7/site-packages']

4. Copy MakePi.py, MakePiPlanes.py into one of the directories in the system path

cp nci.py MakePi.py, MakePiPlanes.py /Applications/PyMOL.app/Contents/lib/python3.7/site-packages

5. Open PyMOL, if the scripts were successfully incorporated you should see the following lines at the top of your PyMOL GUI with no errors


PyMOL>import nci PyMOL>import MakePi

6. Change directories into the local directory where PDB.pdb PDB-dens.cube and PDB-grad.cube files are located

7. Load PDB.pdb, PDB-dens.cube, and PDB-grad.cube into PyMOL

8. To draw pi-contacts call the following

PyMOL>makepi PDB.pdb

9. To display NCIplot isosurfaces call the following

PyMOL>nci PDB

10. To save PyMOL session got to File-> Save Session As 
