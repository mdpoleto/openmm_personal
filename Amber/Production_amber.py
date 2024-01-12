#!/usr/bin/python
###########################################
# Runs a NPT equilibration with CHARMM36m FF.
#
# Outputs DCD files, last frame (CRD) and RST files with positions and
# velocities. Requires a PSF and a CRD. It allows restarting a
# simulation at any given point (using .chk). It also  makes a backup of
# restarted .log file in order to save previous saved energies.
#
# USAGE: python EquNPT_1equ.py -h
#
# mdpoleto@vt.edu -> version 2022
###########################################
import openmm as mm
from openmm import *
from openmm.app import *
from openmm.unit import *

from sys import stdout, exit, stderr
import os, math, fnmatch
import argparse



# Parse user input and options
ap = argparse.ArgumentParser(description=__doc__)

# Mandatory
ap.add_argument('-crd', type=str, default=None, required=True,
                help='Input coordinate file (.crd)')
ap.add_argument('-prmtop', type=str, default=None, required=True,
                help='Topology file in XPLOR format (.psf)')
ap.add_argument('-pdb', type=str, default=None, required=True,
                help='Input coordinate file (.pdb)')
ap.add_argument('-state', type=str, required=True,
                help='XML file to read positions/velocities from (.rst).')

# Options
ap.add_argument('-outname', type=str, default="output",
                help='Default name for output files. Default is "output".')

ap.add_argument('-runtime', default=100, type=float,
                help='Simulation length of each stride (in ps). Default 100.')
ap.add_argument('-nstride', default=1, type=int,
                help='Number of strides/chunks in which the simulation will be splitted. Default 1.')
ap.add_argument('-dt', default=2, type=int,
                help='Integration step (in fs). Default 1.')

ap.add_argument('-savefreq', type=float, default=5,
                help='Frequency (in ps) to save coordinates, checkpoints and trajectory.')
ap.add_argument('-printfreq', type=float, default=5,
                help='Frequency (in ps) to print and write in .log file.')

ap.add_argument('-firststride', type=int, default=1,
                help='First stride number. Default 1.')

ap.add_argument('-temp', default=298, type=float,
                help='Target temperature, in Kelvin. Default is 298.')
ap.add_argument('-pressure', default=1.0, type=float,
                help='Target pressure, in bar. Default is 1.0.')

ap.add_argument('-fc_dih', default=4, type=int,
                help='Force constant for dihedral restraint (kJ/mol/nm^2). Default: 4')
ap.add_argument('-dih_restraint_file', default=None, type=str,
                help='File containing dihedrals to restraint.')

cmd = ap.parse_args()

with open('Production.dat', 'w') as out:
	out.write(' '.join(sys.argv[1:]))

#############################################

jobname = cmd.outname

stride_time = cmd.runtime*picosecond		# in ps
dt = cmd.dt*femtosecond						#fs
nstride = cmd.nstride

print_freq  = cmd.printfreq*picosecond
savcrd_freq = cmd.savefreq*picosecond

temperature = cmd.temp*kelvin
pressure	= cmd.pressure*bar

nsteps = int(stride_time.value_in_unit(picosecond)/dt.value_in_unit(picosecond))
nprint  = int(print_freq.value_in_unit(picosecond)/dt.value_in_unit(picosecond))
nsavcrd = int(savcrd_freq.value_in_unit(picosecond)/dt.value_in_unit(picosecond))

#############################################
# Defining functions to use below:
class ForceReporter(object):
	def __init__(self, file, reportInterval, forcegroup=1, append=False):
		self._out = open(file, 'w')
		self._reportInterval = reportInterval
		self._fgroup = forcegroup
		self._append = False
		self._out = open(file, 'w')

	def __del__(self):
		self._out.close()

	def describeNextReport(self, simulation):
		steps = self._reportInterval - simulation.currentStep%self._reportInterval
		return (steps, False, False, True, False, None)

	def report(self, simulation, state):
		if not self._append:
			self._out.write("#Step\tForce energy\n")
			self._append = True
			try:
				self._out.flush()
			except AttributeError:
				pass

		fenergy = simulation.context.getState(getEnergy=True, groups={self._fgroup})
		step    = simulation.currentStep

		self._out.write('%d\t%f\n' % (step, fenergy.getPotentialEnergy().value_in_unit(kilojoules/mole)))
		try:
			self._out.flush()
		except AttributeError:
			pass

def backup_old_log(pattern, string):
	result = []
	for root, dirs, files in os.walk("./"):
		for name in files:
			if fnmatch.fnmatch(name, pattern):

				try:
					number = int(name[-2])
					avail = isinstance(number, int)
					#print(name,avail)
					if avail == True:
						result.append(number)
				except:
					pass

	if len(result) > 0:
		maxnumber = max(result)
	else:
		maxnumber = 0

	backup_file = "\#" + string + "." + str(maxnumber + 1) + "#"
	os.system("mv " + string + " " + backup_file)
	return backup_file

def dih_restraints(system, crd, fc_dih, restraint_file):

	# dihedral restraints
	dihres = CustomTorsionForce('fc_dih*max(0, abs(diff+wrap) - rwidth)^2; \
	                                 wrap = 2*pi*(step(-diff-pi)-step(diff-pi)); \
	                                 diff = theta - rtheta0; \
	                                 rtheta0 = theta0*pi/180; \
	                                 rwidth = width*pi/180;')
	dihres.addGlobalParameter('fc_dih', fc_dih)
	dihres.addGlobalParameter('pi', 3.141592653589793)
	dihres.addPerTorsionParameter('width')
	dihres.addPerTorsionParameter('theta0')
	if os.path.isfile(restraint_file):
		for line in open(restraint_file, 'r'):
			if line[0] != "#" and line[0] != "!" and line[0] != "@":
				segments = line.strip().split()
				atom1  = int(segments[0])
				atom2  = int(segments[1])
				atom3  = int(segments[2])
				atom4  = int(segments[3])
				theta0 = float(segments[4])
				width  = float(segments[5])
				dihres.addTorsion(atom1,atom2,atom3,atom4,[width,theta0])
	system.addForce(dihres)

	return system

##############################################

#############################################
print("\n> Simulation details:\n")
print("\tJob name = " + jobname)
print("\tCoordinate file = " + str(cmd.crd))
print("\tPDB file = " + str(cmd.pdb))
print("\tTopology file = " + str(cmd.prmtop))

print("\n\tSimulation_time = " + str(stride_time*nstride))
print("\tIntegration timestep = " + str(dt))
print("\tTotal number of steps = " +  str(nsteps*nstride))
print("\tNumber of strides = " + str(cmd.nstride) + " (" + str(stride_time) + " in each stride)")

print("\n\tSave coordinates each " + str(savcrd_freq))
print("\tSave checkpoint each " + str(savcrd_freq))
print("\tPrint in log file each " + str(print_freq))

print("\n\tTemperature = " + str(temperature))
print("\tPressure = " + str(pressure))
#############################################

print("\n> Setting the system:\n")

print("\t- Reading topology and structure file...")
prmtop = AmberPrmtopFile(cmd.prmtop)
inpcrd = AmberInpcrdFile(cmd.crd)

print("\t- Creating system and setting parameters...")
nonbondedMethod = PME
nonbondedCutoff = 1.0*nanometers
ewaldErrorTolerance = 0.0005
rigidWater = True
constraintTolerance = 0.000001
friction = 1.0
system = prmtop.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,
                           constraints=HBonds, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)

if dih_restraint_file != None:
	print("\t- Applying dihedral restraints... (using " + str(dih_restraint_file) + ")")
	system = dih_restraints(system, inpcrd, fc_dih, dih_restraint_file)
else:
	pass


print("\t- Setting barostat...")
system.addForce(MonteCarloBarostat(pressure, temperature))

##############################################################################
######################## FLAT-BOTTOM RESTRAINT SECTION #######################
##############################################################################
flat_constant = 1000 # kJ/mol/nm^2
print("\t- Setting Flat-bottom potential for ligand distances (constant = " + str(flat_constant) + " kJ/mol/nm^2)...")
# 0-based list
BHET_O4  = 3827 # 3828
BHET_CX2 = 3819 # 3820
Y87_HN   = 797  # 798-res59
M161_HN  = 1915 # 1916-res133
S160_OG  = 1910 # 1911-res132
H237_NE2 = 2990 # 2991-res209
S160_HG  = 1911 # 1912-res132
H237_HD1 = 2986 # 2988-res209
D206_CG  = 2534 # 2535-res178

flat_bond_force = CustomBondForce('step(r-r0)*(k/2)*(r-r0)^2')
flat_bond_force.addPerBondParameter('r0')
flat_bond_force.addPerBondParameter('k')
flat_bond_force.setName("Flat-bond-force")
flat_bond_force.setUsesPeriodicBoundaryConditions(True)

flat_bond_force.addBond(BHET_O4, Y87_HN,   [2.00*angstrom, flat_constant*kilojoule_per_mole/nanometer**2])
flat_bond_force.addBond(BHET_O4, M161_HN,  [2.00*angstrom, flat_constant*kilojoule_per_mole/nanometer**2])
flat_bond_force.addBond(S160_OG, BHET_CX2, [2.50*angstrom, flat_constant*kilojoule_per_mole/nanometer**2])
flat_bond_force.addBond(S160_HG, H237_NE2, [2.00*angstrom, flat_constant*kilojoule_per_mole/nanometer**2])
flat_bond_force.addBond(H237_HD1, D206_CG, [2.25*angstrom, flat_constant*kilojoule_per_mole/nanometer**2])
system.addForce(flat_bond_force)

fbout = jobname + "_restraints.dat"

for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)
        if f.getName() == "Flat-bond-force":
                flatbond_index = i

##############################################################################
##############################################################################
##############################################################################


print("\t- Setting integrator...")
integrator = LangevinIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

print('\t- Using platform:', simulation.context.getPlatform().getName())

# Opening a loop of extension NSTRIDE to simulate the entire STRIDE_TIME*NSTRIDE
for n in range(cmd.firststride, nstride + 1):

	print("\n\n>>> Simulating Stride #" + str(n) + " <<<")

	dcd_file = jobname + "_" + str(n) + ".dcd"
	chk_file = jobname + "_" + str(n) + ".chk"
	log_file = jobname + "_" + str(n) + ".log"
	rst_file = jobname + "_" + str(n) + ".rst"
	prv_rst_file = jobname + "_" + str(n-1) + ".rst"
	pdb_file = jobname + "_" + str(n) + ".pdb"
	fbout = jobname + "_" + str(n) + "_restraints.dat"

	if os.path.exists(rst_file):
		print("> Stride #" + str(n) + " finished (" + rst_file + " present). Moving to next stride... <")
		continue


	if os.path.exists(chk_file):
		simulation.loadCheckpoint(chk_file)

		print("> Restarting from checkpoint > " + chk_file + " <")

		# Calculating elapsed time
		chk_time = simulation.context.getState().getTime()
		chk_time_val = chk_time.value_in_unit(picosecond)
		chk_step = math.ceil(chk_time.in_units_of(picosecond)/dt.in_units_of(picosecond))
		print("\t- Elapsed simulation time = " + str(round(chk_time_val,2)*picosecond) + " (step = " + str(chk_step) + ")")

		# Calculating remaining running time
		remaining_time = (stride_time*n - chk_time).in_units_of(picosecond)
		remaining_nsteps = int(math.ceil(remaining_time.in_units_of(picosecond)/dt.in_units_of(picosecond)))

		if remaining_nsteps < 0:
			sys.exit("\n>>>> WARNING: checkpoint last step saved is beyond the last nstep requested. Maybe increase simulation time?<<<<<<\n")

		# Adjust remaining running time
		simulation.currentStep = chk_step
		simulation.context.setTime(chk_time.value_in_unit(picosecond))
		print("\t- Restarting from step " + str(chk_step) + " ...")

		dcd=DCDReporter(dcd_file, nsavcrd, append=True)
		print("\t- Appending to file " + jobname + ".dcd ...")

		backup_file = backup_old_log("*" + log_file + "*", log_file)
		print("\t- Backuping old log to " + backup_file + " ...")

		simulation.reporters.append(dcd)
		simulation.reporters.append(StateDataReporter(stdout, nprint, step=True, speed=True, progress=True, totalSteps=(nsteps*nstride), remainingTime=True, separator='\t\t'))
		simulation.reporters.append(StateDataReporter(log_file, nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))
		simulation.reporters.append(CheckpointReporter(chk_file, nsavcrd))
		simulation.reporters.append(ForceReporter(fbout, nsavcrd,forcegroup=flatbond_index, append=True))

		print("\n> Simulating " + str(remaining_nsteps) + " steps... (Stride #" + str(n) + ")")
		simulation.step(remaining_nsteps)

	else:
		if n == cmd.firststride:
			equil_rst_file = cmd.state
			print("\n> Loading previous state from equilibration > " + equil_rst_file + " <")
			with open(equil_rst_file, 'r') as f:
				simulation.context.setState(XmlSerializer.deserialize(f.read()))
				currstep = int((n-1)*nsteps)
				currtime = currstep*dt.in_units_of(picosecond)
				simulation.currentStep = currstep
				simulation.context.setTime(currtime)
				print("> Current time: " + str(currtime) + " (Step = " + str(currstep) + ")")

		else:
			print("> Loading previous state from > " + prv_rst_file + " <")
			with open(prv_rst_file, 'r') as f:
				simulation.context.setState(XmlSerializer.deserialize(f.read()))
				currstep = int((n-1)*nsteps)
				currtime = currstep*dt.in_units_of(picosecond)
				simulation.currentStep = currstep
				simulation.context.setTime(currtime)
				print("> Current time: " + str(currtime) + " (Step = " + str(currstep) + ")")


		dcd = DCDReporter(dcd_file, nsavcrd)
		firstdcdstep = (currstep) + nsavcrd
		dcd._dcd = DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(), firstdcdstep, nsavcrd) # charmm doesn't like first step to be 0

		simulation.reporters.append(dcd)
		simulation.reporters.append(StateDataReporter(stdout, nprint, step=True, speed=True, progress=True, totalSteps=(nsteps*nstride), remainingTime=True, separator='\t\t'))
		simulation.reporters.append(StateDataReporter(log_file, nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))
		simulation.reporters.append(CheckpointReporter(chk_file, nsavcrd))
		simulation.reporters.append(ForceReporter(fbout, nsavcrd,forcegroup=flatbond_index, append=False))

		print("\n> Simulating " + str(nsteps) + " steps... (Stride #" + str(n) + ")")
		simulation.step(nsteps)

	simulation.reporters.clear() # remove all reporters so the next iteration don't trigger them.


	##################################
	# Writing last frame information of stride
	print("\n> Writing state file (" + str(rst_file) + ")...")
	state = simulation.context.getState( getPositions=True, getVelocities=True )
	with open(rst_file, 'w') as f:
		f.write(XmlSerializer.serialize(state))

	last_frame = int(nsteps/nsavcrd)
	print("> Writing coordinate file (" + str(pdb_file) + ", frame = " + str(last_frame) + ")...")
	positions = simulation.context.getState(getPositions=True).getPositions()
	PDBFile.writeFile(simulation.topology, positions, open(pdb_file, 'w'))

try:
	quote = Quotes.getquote()
	print(quote)
except:
	pass
print("\n> Finished!\n")
