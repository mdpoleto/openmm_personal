#!/usr/bin/python
###########################################
# Creates a restraint file to be used in positional restraints
# Atom naming is CHARMM compatible!!! Bear in mind your atom naming!
#
# USAGE: python make_restraint.py -crd filename.crd -pdb filename.pdb -addsegments HETA
#
# By default, all Drude particles and lone pairs are ignored.
#
# Optionally, you can add segments to the restraint list!
#
# mdpoleto@vt.edu -> version 2023
###########################################
import os
import sys
import argparse

############################################################
bb_ref_list = ["N", "CA", "C", "O", "O5'", "C5'", "O3'", "O1P", "O2P",
               "P", "C4'", "C3'", "C2'", "C1'", "O4'", "O3'", "O2'", "O1'"]

sc_ref_list = ["CB", "CG", "CG1", "CG2", "OG", "OG1", "CD", "CD1", "CD2", "SD",
               "OD1", "OD2", "ND1", "ND2", "CE", "CE1", "CE2", "CE3", "NE",
               "NE1", "NE2", "OE1", "OE2", "CZ", "CZ2", "CZ3", "NZ", "OH","CH2",
               "NH2", "NH1", "SG", "OT1", "OT2", "N1", "C6", "C2", "O2", "N3",
               "C4", "N4", "C5", "O4", "C5M", "N9", "N2", "O6", "N7", "C8", "N6"]

res_ref_list = ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "SER",
                "THR", "ASN", "GLN", "CYS", "GLY", "PRO", "ARG", "HIS", "LYS",
                "ASP", "GLU", "HSD", "HIE", "HSE", "HSP", "CYX", "HID", "GUA",
                "ADE", "URA", "THY", "CYT"]

ions_ref_list = ["POT", "CLA", "MG", "CAL", "SOD"]
exclude_list  = ["SWM4", "TIP3", "TIP", "WAT"]
############################################################

# Parse user input and options
ap = argparse.ArgumentParser(description=__doc__)

# Mandatory
ap.add_argument('-pdb', type=str, default=None, required=False,
                help='Input coordinate file (.pdb)')

ap.add_argument('-crd', type=str, default=None, required=False,
                help='Input coordinate file (.crd)')

ap.add_argument('-addsegments', type=str, default=None, required=False,
                help='Other segments to restraint as BB (heavy atoms only). Use comma for multiple segments.')

cmd = ap.parse_args()

file = cmd.pdb
filecrd = cmd.crd

if cmd.addsegments is None:
	segments = []
else:
	segments = cmd.addsegments.strip(" ").split(",")

################################################
def use_pdb(file,segments):
	f = open(file,"r")

	line = ""

	bb_list = []
	sc_list = []
	addsegment_list = []

	print("\n>>> Reading " + file + " file...")
	while True:
		line = f.readline()
		if not line: break

		if line.split()[0] == "ATOM" or line.split()[0] == "HETATM":
			atom_nr = int(line[6:11].strip(" "))
			atom_name = line[12:16].strip(" ")
			resname = line[17:21].strip(" ")
			resnumber = line[22:26].strip(" ")
			segid = line[72:76]

			atom_nr = str(atom_nr - 1)

			# ignoring Drude and lone pairs
			if atom_name[0] == "D" or atom_name[0:2] == "LP":
				continue
			else:

				if segid not in segments:
					# Ignoring water molecules based on those residue names
					if resname in exclude_list:
						continue

					else:
						if resname in res_ref_list:
							if atom_name in bb_ref_list:
								bb_list.append(atom_nr)
							if atom_name in sc_ref_list:
								sc_list.append(atom_nr)
							if atom_name in ions_ref_list and atom_name[0] != "D" and atom_name[0:2] != "LP":
								bb_list.append(atom_nr)
						else:
							if resname not in exclude_list:
								if segid not in segments:
									#print(atom_nr, atom_name,resnumber, segid)
									print(">>>>>> Residue " + str(resname+"-"+resnumber) + " not found in the RESIDUES/IONS list! Check your input...")
				else:

					# Check for additional segments
					if atom_name[0] != "H" and atom_name[0] != "D" and atom_name[0:2] != "LP":
						# excluding OM atoms from SWM4
						if resname == "SWM4" and atom_name == "OM":
							pass
						else:
							addsegment_list.append(atom_nr)

	return bb_list, sc_list, addsegment_list

def use_crd(file,segments):
	f = open(file,"r")

	line = ""

	bb_list = []
	sc_list = []
	addsegment_list = []

	flag = False

	print("\n>>> Reading " + file + " file...")
	while True:
		line = f.readline()
		if not line: break
		if len(line.split()) > 1:

			if "*" in line:
				continue
			elif "EXT" in line.split():
				readflag = True
				continue

			if readflag == True:
				atom_nr   = int(line[0:10])
				globalres = int(line[11:20])
				resname   = line[21:30].strip()
				atom_name = line[31:40].strip()
				xcoor     = float(line[41:60])
				ycoor     = float(line[61:80])
				zcoor     = float(line[81:100])
				segid     = line[102:107].strip()
				resnumber = line[108:120].strip()
				bfac      = 0.0

				#print(atom_nr, resnumber, resname, atom_name, segid)

				atom_nr = str(atom_nr - 1)

				# ignoring Drude and lone pairs
				if atom_name[0] == "D" or atom_name[0:2] == "LP":
					continue
				else:

					if segid not in segments:
						# Ignoring water molecules based on those residue names
						if resname == "TIP" or resname == "TIP3" or resname == "SWM4":
							continue

						else:
							if resname in res_ref_list:
								if atom_name in bb_ref_list:
									bb_list.append(atom_nr)
								if atom_name in sc_ref_list:
									sc_list.append(atom_nr)
								if atom_name in ions_ref_list and atom_name[0] != "D" and atom_name[0:2] != "LP":
									bb_list.append(atom_nr)
							else:
								if resname not in exclude_list:
									if segid not in segments:
										#print(atom_nr, atom_name,resnumber, segid)
										print(">>>>>> Residue " + str(resname+"-"+resnumber) + " not found in the RESIDUES/IONS list! Check your input...")
					else:

						# Check for additional segments
						if atom_name[0] != "H" and atom_name[0] != "D" and atom_name[0:2] != "LP":
							# excluding OM atoms from SWM4
							if resname == "SWM4" and atom_name == "OM":
								pass
							else:
								addsegment_list.append(atom_nr)

		else:
			continue

	return bb_list, sc_list, addsegment_list

################################################
if filecrd != None:
	bb_list, sc_list, addsegment_list = use_crd(filecrd,segments)
else:
	bb_list, sc_list, addsegment_list = use_pdb(file,segments)


o = open("restraint_atoms.dat", "w")

for index in bb_list:
	string = " " + index + " BB\n"
	o.write(string)

for index in sc_list:
	string = " " + index + " SC\n"
	o.write(string)

for index in addsegment_list:
	string = " " + index + " BB\n"
	o.write(string)

o.close()
print(">>> Done!")
