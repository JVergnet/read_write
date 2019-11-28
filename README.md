# read_write
Read wrtite vasp inputs 
convenient proxies for pymatgen & vasp on slurm architecture

----------------------------------------------------------

Contains 4 main functionalties that eanble the following workflow :

*Write_job* : parse Cif / poscar, modify the structure and write corresponding vasp input files 

*Launch_job* : launch one or several vasp jobs from specified folders 

*Read_job* : parse the output of the run (vasprun.xml) and post-run analysis 

*Rerun_job* : parse either vasprun or poscar and write a modified copy vasp input files elsewhere


Additionally : *Read_run_app* is a visual interface which allows to perform this workflow from gui

----------------------------------------------------------

These high-level scripts rely on 2 main subpackages : 

 -- pre-run analysis, involving only the structure (subpackage structural_analysis)

 -- post-run analysis, which involves also the electronic properties (subpackage electronic _analysis)

 -- utils & filtering subpackages help for all these functions 

