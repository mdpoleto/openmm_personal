#!/usr/bin/python
###########################################
# Runs a MD simulation with C36 FF.
#
# Outputs DCD files, last frame (CRD) and RST files with positions and
# velocities. Requires a PSF, a CRD and RST from equilibration. It allows
# restarting a simulation at any given point (using both .chk and .rst files).
# It also makes a backup of restarted .log file in order to save previous saved
# energies.
#
# USAGE: python Production_c36.py -h
#
# mdpoleto@vt.edu -> version 2023
###########################################
import openmm as mm
from openmm import *
from openmm.app import *
from openmm.unit import *
from openmmplumed import *
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

ap.add_argument('-plumed', default=None, type=str,
                help='Plumed script to use Plumed with OpenMM')


cmd = ap.parse_args()

with open('Production_c36.dat', 'w') as out:
	out.write(' '.join(sys.argv[1:]))

#############################################
# Although we use checkpoints to restart simulations, an unexpected crash may
# harm the dcd integrity beyond repair. It is rare, but that may happenself.
# Therefore, we use 10 simulation strides over a loop, creating 10 dcd files that
# are concatenated afterwards.

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

plumedscript = cmd.plumed

#############################################
# Defining functions to use below:
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

def get_cubic_box(psf, rstfile):

	f = open(rstfile, 'r')

	box = {}

	while True:
		line = f.readline()
		if not line: break

		if line.split()[0] == "<A":
			size = line.split()[1].strip('x="')
			box['A'] = float(size)
		elif line.split()[0] == "<B":
			size = line.split()[2].strip('y="')
			box['B'] = float(size)
		elif line.split()[0] == "<C":
			size = line.split()[3].strip('z="').strip('"/>')
			box['C'] = float(size)
		else:
			pass

	boxX = box['A']*nanometer
	boxY = box['B']*nanometer
	boxZ = box['C']*nanometer

	psf.setBox(boxX, boxY, boxZ)

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
print("\n> Simulation details:\n")
print("\tJob name = " + jobname)
print("\tCRD file = " + str(cmd.crd))
print("\tPSF file = " + str(cmd.psf))
print("\tToppar stream file = " + str(read_toppar(cmd.toppar)[1]))

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
print("\t- Reading force field directory...")
charmm_params = read_toppar(cmd.toppar)[0]

print("\t- Reading topology and structure file...")
psf = CharmmPsfFile(cmd.psf)
crd = CharmmCrdFile(cmd.crd)

print("\t- Setting box (using information on -state file)...")
psf = get_cubic_box(psf, cmd.state)

print("\t- Creating system and setting parameters...")
system = psf.createSystem(charmm_params, nonbondedMethod=PME, nonbondedCutoff=1.2*nanometers, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0005, constraints=HBonds)
#system = vfswitch(system, psf, 1.0, 1.2)  # if turned on, remove switchDistance on the line above

print("\t- Setting barostat...")
system.addForce(MonteCarloBarostat(pressure, temperature))

if plumedscript != None:
	print("\t- Loading Plumed script...")
	with open(plumedscript, 'r') as plumedfile:
		plumed_param = plumedfile.read()
	system.addForce(PlumedForce(plumed_param))
else:
	pass

print('\t- Setting simulation context...')
platform = Platform.getPlatformByName("CUDA")
cuda_devices=os.getenv('CUDA_VISIBLE_DEVICES')
properties = dict(CudaPrecision='single', CudaDeviceIndex=cuda_devices) if platform.getName() == 'CUDA' else dict()

simulation = Simulation(psf.topology, system, integrator, platform, properties)
simulation.context.setPositions(crd.positions)
simulation.context.computeVirtualSites()

print('\t- Using platform:', simulation.context.getPlatform().getName())

# Opening a loop of extension NSTRIDE to simulate the entire STRIDE_TIME*NSTRIDE
for n in range(cmd.firststride, nstride + 1):

	print("\n\n>>> Simulating Stride #" + str(n) + " <<<")

	dcd_file = jobname + "_" + str(n) + ".dcd"
	chk_file = jobname + "_" + str(n) + ".chk"
	log_file = jobname + "_" + str(n) + ".log"
	rst_file = jobname + "_" + str(n) + ".rst"
	prv_rst_file = jobname + "_" + str(n-1) + ".rst"
	crd_file = jobname + "_" + str(n) + ".crd"

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

		print("\n> Simulating " + str(nsteps) + " steps... (Stride #" + str(n) + ")")
		simulation.step(nsteps)

	simulation.reporters.clear() # remove all reporters so the next iteration don't trigger them.


	##################################
	# Writing last frame information of stride
	print("\n> Writing stride state file (" + str(rst_file) + ")...")
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
