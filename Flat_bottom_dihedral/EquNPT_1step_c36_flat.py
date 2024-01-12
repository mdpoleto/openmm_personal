#!/usr/bin/python
###########################################
# Runs a NPT equilibration with CHARMM36m FF.
#
# Outputs DCD files, last frame (CRD) and RST files with positions and
# velocities. Requires a PSF and a CRD. It allows restarting a
# simulation at any given point (using .chk). It also  makes a backup of
# restarted .log file in order to save previous saved energies.
#
# USAGE: python EquNPT_c36_1step.py -h
#
# mdpoleto@vt.edu -> version 2023
###########################################
import simtk.openmm as mm
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
import MDAnalysis as mda
import warnings
import parmed as pmd
import time

from sys import stdout, exit, stderr
import os, math, fnmatch
import argparse

warnings.filterwarnings("ignore", message="Found no information for attr:")
warnings.filterwarnings("ignore", message="Found missing chainIDs")
warnings.filterwarnings("ignore", message="Supplied AtomGroup was missing the following attributes")


# Parse user input and options
ap = argparse.ArgumentParser(description=__doc__)

# Mandatory
ap.add_argument('-crd', type=str, default=None, required=True,
                help='Input coordinate file (.crd)')
ap.add_argument('-psf', type=str, default=None, required=True,
                help='Topology file in XPLOR format (.psf)')
ap.add_argument('-toppar', type=str, default='toppar.str', required=True,
                help='Force field stream file (ex. "toppar.str").')
ap.add_argument('-restraint_file', default="restraint.dat", type=str, required=True,
                help='File containing heavy atoms to restraint separeted by BB and SC.')
ap.add_argument('-boxsize', type=float, required=True, nargs=3,
                help='Box vector lengths (a,b,c) in nanometers')

# Options
ap.add_argument('-outname', type=str, default="output.c36.npt",
                help='Default name for output files. Default is "output".')

ap.add_argument('-runtime', default=100, type=float,
                help='Simulation length (in ps). Default = 100.')
ap.add_argument('-dt', default=2, type=int,
                help='Integration step (in fs). Default = 1.')
ap.add_argument('-savefreq', type=float, default=5,
                help='Frequency (in ps) to save coordinates, checkpoints and trajectory. Default = 5')
ap.add_argument('-printfreq', type=float, default=5,
                help='Frequency (in ps) to print and write in .log file. Default = 5')

ap.add_argument('-temp', default=298, type=float,
                help='Target temperature, in Kelvin. Default = 298.')
ap.add_argument('-pressure', default=1.0, type=float,
                help='Target pressure, in bar. Default is 1.0.')

ap.add_argument('-fc_bb', default=500, type=int,
                help='Force constant for backbone heavy atoms (kJ/mol/nm^2). Default: 500')
ap.add_argument('-fc_sc', default=500, type=int,
                help='Force constant for sidechain heavy atoms (kJ/mol/nm^2). Default: 500')
ap.add_argument('-fc_dih', default=4, type=int,
                help='Force constant for dihedral restraint (kJ/mol/nm^2). Default: 4')
ap.add_argument('-dih_restraint_file', default=None, type=str,
                help='File containing dihedrals to restraint.')


ap.add_argument('-state_pequil', type=str, required=False, default=None,
                help='XML file to read positions/velocities from previous equilibration (.rst).')

cmd = ap.parse_args()

with open('Equilibration_NPT.dat', 'w') as out:
	out.write(' '.join(sys.argv[1:]))

#############################################

jobname = cmd.outname

simulation_time = cmd.runtime*picosecond		# in ps
dt = cmd.dt*femtosecond						#fs

print_freq  = cmd.printfreq*picosecond
savcrd_freq = cmd.savefreq*picosecond

temperature = cmd.temp*kelvin
pressure	= cmd.pressure*bar

fc_bb = float(cmd.fc_bb)
fc_sc = float(cmd.fc_sc)
fc_dih = float(cmd.fc_dih)
restraint_file = cmd.restraint_file
dih_restraint_file = cmd.dih_restraint_file

nsteps = int(simulation_time.value_in_unit(picosecond)/dt.value_in_unit(picosecond))
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
			self._out.write("#Step\tForce energy (kJ/mol)\n")
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

def gen_box(psf, boxsize):

	boxvectors = boxsize

	boxlx = boxvectors[0]*nanometer
	boxly = boxvectors[1]*nanometer
	boxlz = boxvectors[2]*nanometer

	psf.setBox(boxlx, boxly, boxlz)
	return psf

def read_toppar(filename):
	extlist = ['rtf', 'prm', 'str']

	parFiles = ()
	for line in open(filename, 'r'):
		if '!' in line: line = line.split('!')[0]
		parfile = line.strip()
		if len(parfile) != 0:
			ext = parfile.lower().split('.')[-1]
			if not ext in extlist: continue
			parFiles += ( parfile, )

	params = CharmmParameterSet( *parFiles )
	return params, parFiles

def restraints(system, crd, fc_bb, fc_sc, restraint_file):

	boxlx = system.getDefaultPeriodicBoxVectors()[0][0].value_in_unit(nanometers)
	boxly = system.getDefaultPeriodicBoxVectors()[1][1].value_in_unit(nanometers)
	boxlz = system.getDefaultPeriodicBoxVectors()[2][2].value_in_unit(nanometers)

	if fc_bb > 0 or fc_sc > 0:
		# positional restraints for protein
		posresPROT = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2;')
		posresPROT.addPerParticleParameter('k')
		posresPROT.addPerParticleParameter('x0')
		posresPROT.addPerParticleParameter('y0')
		posresPROT.addPerParticleParameter('z0')
		for line in open(restraint_file, 'r'):
			segments = line.strip().split()
			atom1 = int(segments[0])
			state = segments[1]
			xpos  = crd.positions[atom1].value_in_unit(nanometers)[0]
			ypos  = crd.positions[atom1].value_in_unit(nanometers)[1]
			zpos  = crd.positions[atom1].value_in_unit(nanometers)[2]
			if state == 'BB' and fc_bb > 0:
				fc_ppos = fc_bb
				posresPROT.addParticle(atom1, [fc_ppos, xpos, ypos, zpos])
			if state == 'SC' and fc_sc > 0:
				fc_ppos = fc_sc
				posresPROT.addParticle(atom1, [fc_ppos, xpos, ypos, zpos])
		system.addForce(posresPROT)

	return system

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
	dihres.setName("Dihedral-restraint-force")
	for line in open(restraint_file, 'r'):
		if line.strip(): #skip empty lines
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

def vfswitch(system, psf, r_on, r_off):

	# custom nonbonded force for force-switch
	chknbfix = False
	for force in system.getForces():
		if isinstance(force, NonbondedForce):
			nonbonded = force
		if isinstance(force, CustomNonbondedForce):
			nbfix     = force
			chknbfix  = True

	# vfswitch
	vfswitch = CustomNonbondedForce('step(Ron-r)*(ccnba*tr6*tr6-ccnbb*tr6+ccnbb*onoff3-ccnba*onoff6) \
	                                 +step(r-Ron)*step(Roff-r)*(cr12*rjunk6 - cr6*rjunk3) \
	                                 -step(r-Ron)*step(Ron-r)*(cr12*rjunk6 - cr6*rjunk3); \
	                                 cr6  = ccnbb*ofdif3*rjunk3; \
	                                 cr12 = ccnba*ofdif6*rjunk6; \
	                                 rjunk3 = r3-recof3; \
	                                 rjunk6 = tr6-recof6; \
	                                 r3 = r1*tr2; \
	                                 r1 = sqrt(tr2); \
	                                 tr6 = tr2 * tr2 * tr2; \
	                                 tr2 = 1.0/s2; \
	                                 s2 = r*r; \
	                                 ccnbb = 4.0*epsilon*sigma^6; \
	                                 ccnba = 4.0*epsilon*sigma^12; \
	                                 sigma = sigma1+sigma2; \
	                                 epsilon = epsilon1*epsilon2; \
	                                 onoff3 = recof3/on3; \
	                                 onoff6 = recof6/on6; \
	                                 ofdif3 = off3/(off3 - on3); \
	                                 ofdif6 = off6/(off6 - on6); \
	                                 recof3 = 1.0/off3; \
	                                 on6 = on3*on3; \
	                                 on3 = c2onnb*Ron; \
	                                 recof6 = 1.0/off6; \
	                                 off6 = off3*off3; \
	                                 off3 = c2ofnb*Roff; \
	                                 c2ofnb = Roff*Roff; \
	                                 c2onnb = Ron*Ron; \
	                                 Ron  = %f; \
	                                 Roff = %f;' % (r_on, r_off) )
	vfswitch.addPerParticleParameter('sigma')
	vfswitch.addPerParticleParameter('epsilon')
	vfswitch.setNonbondedMethod(vfswitch.CutoffPeriodic)
	vfswitch.setCutoffDistance(nonbonded.getCutoffDistance())
	for i in range(nonbonded.getNumParticles()):
		chg, sig, eps = nonbonded.getParticleParameters(i)
		nonbonded.setParticleParameters(i, chg, 0.0, 0.0) # zero-out LJ
		sig = sig*0.5
		eps = eps**0.5
		vfswitch.addParticle([sig, eps])
	for i in range(nonbonded.getNumExceptions()):
		atom1, atom2 = nonbonded.getExceptionParameters(i)[:2]
		vfswitch.addExclusion(atom1, atom2)
	vfswitch.setForceGroup(psf.NONBONDED_FORCE_GROUP)
	system.addForce(vfswitch)

	# vfswitch14
	vfswitch14 = CustomBondForce('step(Ron-r)*(ccnba*tr6*tr6-ccnbb*tr6+ccnbb*onoff3-ccnba*onoff6) \
	                              +step(r-Ron)*step(Roff-r)*(cr12*rjunk6 - cr6*rjunk3) \
	                              -step(r-Ron)*step(Ron-r)*(cr12*rjunk6 - cr6*rjunk3); \
	                              cr6  = ccnbb*ofdif3*rjunk3; \
	                              cr12 = ccnba*ofdif6*rjunk6; \
	                              rjunk3 = r3-recof3; \
	                              rjunk6 = tr6-recof6; \
	                              r3 = r1*tr2; \
	                              r1 = sqrt(tr2); \
	                              tr6 = tr2 * tr2 * tr2; \
	                              tr2 = 1.0/s2; \
	                              s2 = r*r; \
	                              ccnbb = 4.0*epsilon*sigma^6; \
	                              ccnba = 4.0*epsilon*sigma^12; \
	                              onoff3 = recof3/on3; \
	                              onoff6 = recof6/on6; \
	                              ofdif3 = off3/(off3 - on3); \
	                              ofdif6 = off6/(off6 - on6); \
	                              recof3 = 1.0/off3; \
	                              on6 = on3*on3; \
	                              on3 = c2onnb*Ron; \
	                              recof6 = 1.0/off6; \
	                              off6 = off3*off3; \
	                              off3 = c2ofnb*Roff; \
	                              c2ofnb = Roff*Roff; \
	                              c2onnb = Ron*Ron; \
	                              Ron  = %f; \
	                              Roff = %f;' % (r_on, r_off) )
	vfswitch14.addPerBondParameter('sigma')
	vfswitch14.addPerBondParameter('epsilon')
	for i in range(nonbonded.getNumExceptions()):
		atom1, atom2, chg, sig, eps = nonbonded.getExceptionParameters(i)
		nonbonded.setExceptionParameters(i, atom1, atom2, chg, 0.0, 0.0) # zero-out LJ14
		vfswitch14.addBond(atom1, atom2, [sig, eps])
	system.addForce(vfswitch14)

	# vfswitch_NBFIX
	if chknbfix:
		nbfix.setEnergyFunction('step(Ron-r)*(ccnba*tr6*tr6-ccnbb*tr6+ccnbb*onoff3-ccnba*onoff6) \
		                         +step(r-Ron)*step(Roff-r)*(cr12*rjunk6 - cr6*rjunk3) \
		                         -step(r-Ron)*step(Ron-r)*(cr12*rjunk6 - cr6*rjunk3); \
		                         cr6  = ccnbb*ofdif3*rjunk3; \
		                         cr12 = ccnba*ofdif6*rjunk6; \
		                         rjunk3 = r3-recof3; \
		                         rjunk6 = tr6-recof6; \
		                         r3 = r1*tr2; \
		                         r1 = sqrt(tr2); \
		                         tr6 = tr2 * tr2 * tr2; \
		                         tr2 = 1.0/s2; \
		                         s2 = r*r; \
		                         ccnbb = bcoef(type1, type2); \
		                         ccnba = acoef(type1, type2)^2; \
		                         onoff3 = recof3/on3; \
		                         onoff6 = recof6/on6; \
		                         ofdif3 = off3/(off3 - on3); \
		                         ofdif6 = off6/(off6 - on6); \
		                         recof3 = 1.0/off3; \
		                         on6 = on3*on3; \
		                         on3 = c2onnb*Ron; \
		                         recof6 = 1.0/off6; \
		                         off6 = off3*off3; \
		                         off3 = c2ofnb*Roff; \
		                         c2ofnb = Roff*Roff; \
		                         c2onnb = Ron*Ron; \
		                         Ron  = %f; \
		                         Roff = %f;' % (r_on, r_off) )

		# turn off long range correction (OpenMM Issues: #2353)
		nbfix.setUseLongRangeCorrection(False)

	return system
##############################################

#############################################
print("\n> NPT Equilibration! Let's do this:")
print("\n> Simulation details:\n")
print("\tJob name = " + jobname)
print("\tCRD file = " + str(cmd.crd))
print("\tPSF file = " + str(cmd.psf))
print("\tToppar stream file = " + str(read_toppar(cmd.toppar)[1]))

print("\n\tSimulation_time = " + str(simulation_time))
print("\tIntegration timestep = " + str(dt))
print("\tTotal number of steps = " +  str(nsteps))

print("\n\tSave coordinates each " + str(savcrd_freq))
print("\tSave checkpoint each " + str(savcrd_freq))
print("\tPrint in log file each " + str(print_freq))

print("\n\tTemperature = " + str(temperature))
print("\tPressure = " + str(pressure) + " (NPT ensemble)")
print("\tBB Restraint constant = " + str(fc_bb) + " (in kJ/mol/nm^2)")
print("\tSC Restraint constant = " + str(fc_sc) + " (in kJ/mol/nm^2)")
if dih_restraint_file != None:
	print("\tDihedral Restraint constant = " + str(fc_dih) + " (kJ/mol/rad^2)")

#############################################

print("\n> NPT Equilibration! Let's do this:\n")
print("\n> Setting the system:\n")
print("\t- Reading force field directory...")
charmm_params = read_toppar(cmd.toppar)[0]

print("\t- Reading topology and structure file...")
psf = CharmmPsfFile(cmd.psf)
crd = CharmmCrdFile(cmd.crd)

print("\t- Setting box (using user's information)...")
psf = gen_box(psf, cmd.boxsize)

print("\t- Creating system and setting parameters...")
system = psf.createSystem(charmm_params, nonbondedMethod=PME, nonbondedCutoff=1.2*nanometers, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0005, constraints=HBonds)
#system = vfswitch(system, psf, 1.0, 1.2)


if dih_restraint_file != None:
	print("\t- Applying dihedral restraints... (using " + str(dih_restraint_file) + ")")
	system = dih_restraints(system, crd, fc_dih, dih_restraint_file)
else:
	print("\t- Applying restraints... (using " + str(restraint_file) + ")")
	system = restraints(system, crd, fc_bb, fc_sc,restraint_file)

for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)
        if f.getName() == "Dihedral-restraint-force":
                flatdih_index = i

print("\t- Setting barostat...")
system.addForce(MonteCarloBarostat(pressure, temperature))


print("\t- Setting integrator...")
integrator = LangevinIntegrator(temperature, 1/picosecond, dt)
simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(crd.positions)

# check if starting from scratch or from previous equilibration
if cmd.state_pequil is None:
	simulation.context.setTime(0.0)
	simulation.currentStep = 0

	print("\t- Energy minimization: 1000 steps")
	simulation.minimizeEnergy(tolerance=10*kilojoule/mole, maxIterations=1000)
	print("\t-> Potential Energy = " + str(simulation.context.getState(getEnergy=True).getPotentialEnergy()))

	print("\t- Setting initial velocities...")
	simulation.context.setVelocitiesToTemperature(temperature)

else:
	print("\n> Loading previous state from equilibration > " + cmd.state_pequil + " <")
	with open(cmd.state_pequil, 'r') as f:
		simulation.context.setState(XmlSerializer.deserialize(f.read()))
	simulation.context.setTime(0.0)
	simulation.currentStep = 0


print('\t- Using platform:', simulation.context.getPlatform().getName())

##########################################

dcd_file = jobname + ".dcd"
chk_file = jobname + ".chk"
log_file = jobname + ".log"
rst_file = jobname + ".rst"
prv_rst_file = jobname + ".rst"
crd_file = jobname + ".crd"
fbout = jobname + "_restraints.dat"

if os.path.exists(chk_file):
	simulation.loadCheckpoint(chk_file)

	print("> Restarting from checkpoint > " + chk_file + " <")

	# Calculating elapsed time
	chk_time = simulation.context.getState().getTime()
	chk_time_val = round(chk_time.value_in_unit(picosecond),4)
	chk_step = math.ceil(chk_time_val/dt.value_in_unit(picosecond))
	print("\t- Elapsed simulation time = " + str(round(chk_time_val,2)*picosecond) + " (step = " + str(chk_step) + ")")

	# Calculating remaining running time
	remaining_time = (simulation_time - chk_time).in_units_of(picosecond)
	remaining_nsteps = int(math.ceil(remaining_time.in_units_of(picosecond)/dt.in_units_of(picosecond)))

	# Adjust remaining running time
	simulation.currentStep = chk_step
	print("\t- Restarting from step " + str(chk_step) + " ...")

	dcd=DCDReporter(dcd_file, nsavcrd, append=True)
	print("\t- Appending to file " + jobname + ".dcd ...")

	backup_file = backup_old_log("*" + log_file + "*", log_file)
	print("\t- Backuping old log to " + backup_file + " ...")

	simulation.reporters.append(dcd)
	simulation.reporters.append(StateDataReporter(stdout, nprint, step=True, speed=True, progress=True, totalSteps=nsteps, remainingTime=True, separator='\t\t'))
	simulation.reporters.append(StateDataReporter(log_file, nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))
	simulation.reporters.append(CheckpointReporter(chk_file, nsavcrd))
	simulation.reporters.append(ForceReporter(fbout, nsavcrd,forcegroup=flatdih_index, append=True))

	print("\n> Simulating " + str(remaining_nsteps) + " steps...")
	simulation.step(remaining_nsteps)

else:
	dcd = DCDReporter(dcd_file, nsavcrd)
	firstdcdstep = (nsteps) + nsavcrd
	dcd._dcd = DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(), firstdcdstep, nsavcrd) # charmm doesn't like first step to be 0

	simulation.reporters.append(dcd)
	simulation.reporters.append(StateDataReporter(stdout, nprint, step=True, speed=True, progress=True, totalSteps=nsteps, remainingTime=True, separator='\t\t'))
	simulation.reporters.append(StateDataReporter(log_file, nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))
	simulation.reporters.append(CheckpointReporter(chk_file, nsavcrd))
	simulation.reporters.append(ForceReporter(fbout, nsavcrd,forcegroup=flatdih_index, append=False))

	print("\n> Simulating " + str(nsteps) + " steps...")
	simulation.step(nsteps)

simulation.reporters.clear() # remove all reporters so the next iteration don't trigger them.


##################################
# Writing last frame information of stride
print("\n> Writing state file (" + str(rst_file) + ")...")
state = simulation.context.getState( getPositions=True, getVelocities=True )
with open(rst_file, 'w') as f:
	f.write(XmlSerializer.serialize(state))

last_frame = int(nsteps/nsavcrd)
print("> Writing last coordinate (" + str(crd_file) + ", frame = " + str(last_frame) + ")...")
u = mda.Universe(cmd.psf, dcd_file)
system = u.select_atoms('all')
for ts in u.trajectory[(len(u.trajectory)-1):len(u.trajectory)]:
	with mda.Writer(str(crd_file), system.n_atoms, extended = True) as W:
		W.write(system)

try:
	quote = Quotes.getquote()
	print(quote)
except:
	pass

print("\n> Finished!\n")
