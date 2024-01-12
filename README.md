# openmm_personal


May 17 2021


Hi

Those are the OpenMM scripts I wrote for OpenMM +v7.6 to make our lives a bit easier. Here you should find:

1- EquNVT_c36.py
2- EquNPT_c36.py
3- EquNPT_c36_1step.py
4- Production_c36.py

which were written using the proper MD parameters for CHARMM36 FF. 


Also, the make_restraint.py is a script written to read a .pdb file and parses Backbone and Sidechain atoms in order to write a restraint_atoms.dat file to be used during equilibration. For more about it, run 

$python make_restraint.py -h


All these files are currently being tested and so far they worked pretty well for my systems. I would really appreciate any comments or suggestions to make them more failproof =)


######################################################
Things to bear in mind:

I   -> These scripts are argument parsed and expect input files to be set in your command line. For example:
    
    $ python Production_c36.py -psf 1aki_full.psf -pdb 1aki.NPT.pdb -toppar toppar_c36.str -state 1aki.NPT.rst -runtime 1000 -nstride 10 -dt 2 -savefreq 10 -printfreq 10 -temp 298 -outname 1aki.prod

the line above executes Production_c36.py using the files set in -psf, -pdb and -rst, using the FF parameters set in -toppar and askes the script to write output files with the prefix "1aki.prod". This line also defines the simulation time to 1000 ps (-runtime) for each of the 10 strides (-nstride), with a total of 10 ns. The time step is defined to 2 fs (-dt), the temperature is defined at 298 K and the coordinates are written every 10 ps. The .log file is also written every 10 ps.

The complete list of flags available within those scripts can be accessed by using -h. For instance, "$python Production_C36.py -h" produces the help below:

usage: Production_c36.py [-h] -pdb PDB -psf PSF -toppar TOPPAR -state STATE
                         [-outname OUTNAME] [-runtime RUNTIME]
                         [-nstride NSTRIDE] [-dt DT] [-savefreq SAVEFREQ]
                         [-printfreq PRINTFREQ] [-temp TEMP]
                         [-pressure PRESSURE]

optional arguments:
  -h, --help            show this help message and exit
  -pdb PDB              Input coordinate file (.pdb)
  -psf PSF              Topology file in XPLOR format (.psf)
  -toppar TOPPAR        Force field stream file (ex. "toppar.str").
  -state STATE          XML file to read positions/velocities from (.rst).
  -outname OUTNAME      Default name for output files. Default is "output".
  -runtime RUNTIME      Simulation length of each stride (in ps). Default 100.
  -nstride NSTRIDE      Number of strides/chunks in which the simulation will
                        be splitted. Default 5.
  -dt DT                Integration step (in fs). Default 1.
  -savefreq SAVEFREQ    Frequency (in ps) to save coordinates, checkpoints and
                        trajectory.
  -printfreq PRINTFREQ  Frequency (in ps) to print and write in .log file.
  -temp TEMP            Target temperature, in Kelvin. Default is 298.
  -pressure PRESSURE    Target pressure, in bar. Default is 1.0.


II  -> I realize there are a lot of flags to be set and you might want to have a backup of the flags you used on a recent run, so the script will automatically write a .dat file containing the most recent command line you executed when running the script.


III   -> Checkpoints (.chk files) are now implemented and it is written at the same frequency as the coordinates (-savefreq). If your run crashes, the script will look for a .chk file with the prefix set in -outname, calculate the last written frame and run the remaining number of steps. Note that if you modify the -outname prefix used on your crashed run, the script will start from scratch, so double check your -outname usage.


IV    -> When using the checkpoint to restart a simulation, the script will automatically backup your .log file to #outname.log.1#. If there is already a #outname.log.1#, then it writes a #outname.log.2# and so on.


V     -> Equilibration scripts DO NOT have -nstride options since they are designed to be short. However, they DO HAVE checkpoint implementations.

VI    -> Boxsize is automatically read from the provided -state file, EXCEPT for EquNVT_c36.py, which contains a -boxsize flag so the user can provide the initial length of box vectors.


######################################################
Example of equilibrating a system on c36 to run a production run on Drude:


$ python EquNVT_c36.py -psf 1aki_full.psf -pdb 1aki.min.pdb -toppar toppar_C36.str -boxsize 7.0 7.0 7.0 -runtime 1000 -dt 2 -savefreq 10 -printfreq 10 -temp 298 -outname 1aki.nvt
[simulation time = 1 ns; time step = 2 fs; temperature 298 K; pressure = 1 bar (default)]


$ python EquNPT_c36.py -psf 1aki_full.psf -pdb 1aki.nvt.pdb -toppar toppar_C36.str -state 1aki.nvt.rst -runtime 1000 -dt 2 -savefreq 10 -printfreq 10 -temp 298 -outname 1aki.npt
[simulation time = 1 ns; time step = 2 fs; temperature 298 K; pressure = 1 bar (default)]




