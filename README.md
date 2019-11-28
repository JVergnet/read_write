# read_write
Read wrtite vasp inputs 
convenient proxies for pymatgen & vasp on slurm architecture
Contains 4 main functionalties that eanble the following workflow :

*Write* : parse Cif / poscar, modify the structure and write corresponding vasp input files 

*Launch* : launch one or several vasp jobs from specified folders 

*Read* : parse the output of the run (vasprun.xml) and post-run analysis 

*Rerun* : parse either vasprun or poscar and write a modified copy vasp input files elsewhere



Additionally : *Read_run_app* is a visual interface which allows to perform this workflow from gui

