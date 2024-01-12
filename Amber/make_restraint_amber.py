#!/usr/bin/python
###########################################
# Creates a restraint file to be used in run_EquilibrationNVT.py
# Atom naming is CHARMM compatible!!! Bear in mind your atom naming!
#
# USAGE: python make_restraint.py input.pdb
#
# Optionally, you can define if ions/waters should be restrained or not!
#
# mdpoleto@vt.edu -> version 2021
###########################################
import os
import sys
import argparse

############################################################
bb_ref_list = ["N", "CA", "C", "O"]
sc_ref_list = ["CB", "CG", "CG1", "CG2", "OG", "OG1", "CD", "CD1", "CD2", "SD",
               "OD1", "OD2", "ND1", "ND2", "CE", "CE1", "CE2", "CE3", "NE",
               "NE1", "NE2", "OE1", "OE2", "CZ", "CZ2", "CZ3", "NZ", "OH",
               "CH2", "NH2", "NH1", "SG", "OT1", "OT2"]

res_ref_list = ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "SER",
                "THR", "ASN", "GLN", "CYS", "GLY", "PRO", "ARG", "HIS", "LYS",
                "ASP", "GLU", "HSD", "HIE", "HSE", "HSP", "CYX", "HID"]

ions_ref_list = ["POT", "CLA", "MG", "CAL"]
exclude_list = ["SWM4", "TIP3", "TIP", "WAT"]
############################################################

# Parse user input and options
ap = argparse.ArgumentParser(description=__doc__)

# Mandatory
ap.add_argument('-pdb', type=str, default=None, required=False,
                help='Input coordinate file (.pdb)')

ap.add_argument('-addresname', type=str, default=None, required=False,
                help='Other segments to restraint as BB (heavy atoms only). Use comma for multiple segments.')

cmd = ap.parse_args()

file = cmd.pdb

if cmd.addresname is None:
	segments = []
else:
	segments = cmd.addresname.strip(" ").split(",")

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

				if resname not in segments:
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


################################################
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
