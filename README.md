This is the find long path algorithm (flp) written by Paul Zierep 
used as core element in the SeMPI pipeline
http://www.pharmaceutical-bioinformatics.de/SeMPI 
in order to create the sub-paths for a series of at least 3800 molecules.

The object-oriented approch allows to modify the code and apply it to different situations.

Requirements:
python 2.7,
numpy (at least version 1.11.0),
rdkit (at least version 2015.03.1)
PIL (version 2.1.0, there seems to be a problem with rdkit and newer versions of PIL, or PILLOW)

This is a python package, which shoud be usable anywhere on a unix or linux related system, just clone or download it 
and try it out. System wide usage can be allowed by including the path to the PYTHONPATH.
An installation script will also be available in future versions.



