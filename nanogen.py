### Written by Matthew Kwiecien Jul, 2015
### Modified by CGG, spring 2016

from __future__ import print_function
from __future__ import division

import numpy as np
from scipy import *
import subprocess
import time
import os
import re
from itertools import chain
import shutil

global basePath
basePath = os.path.abspath('..') + "/Results/"

# Maybe adding later... need to think more about the implementation
# class Nanotube:

# 	def __init__(self, rings, n, m):
# 		self.rings = rings
# 		self.n = n
# 		self.m = m

# Simple function to create a quoted string to protect file/path names with spaces
def quoted(theString):

	return '"' + theString + '"'

# VMD generation of nanotube and periodic boundary conditions
def tubeGen(inFile, pbcFile, N_0, n, m):
	""" tubeGen creates a periodic nanotube with the following input parameters:  inFile 
	is the name of the initial nanotube with number of rings N_0 and dimensions n x m.  pbcFile is 
	the name of the same nanotube but now with periodic boundary conditions applied to it. """

	#Bonds lengths of different armchair nanotubes in nanometers
	s0 = 0.1418
	#calculates the length of the nanotube based on bond lengths
	l = float((N_0-0.75))*s0*np.sqrt(3)

	tubePath = basePath + "cnt" + str(N_0) + "_" + str(n) + "x" + str(m) + "/"
	pbcPath = tubePath + "PBC/" + inFile + "_run/"
	if not os.path.exists(pbcPath):
		os.makedirs(pbcPath)

	# Opening a pipe to VMD in the shell
	VMDin=subprocess.Popen(["vmd","-dispdev", "none"], stdin=subprocess.PIPE)

	# runs CNTtools.tcl script to generate nanotube and generate PBCs
	sourceCNT = "source CNTtools.tcl\n"
	CNTtools = "package require CNTtools 1.0\n"

	genNT = "genNT " + quoted(inFile) + " " + quoted(tubePath) + " " + str(l) + " " + str(n) + " " + str(m) + "\n"
	pbcNT = "pbcNT " + quoted(tubePath + inFile) + " " + quoted(pbcPath + pbcFile) + " default\n"
	fixNT = "fixNT " + quoted(pbcPath + pbcFile) + "\n"
	removeLangevinWater = "removeLangevinWater " + quoted(pbcPath + pbcFile) + "\n"

	# run commands through pipe and saves to file
	VMDin.stdin.write(sourceCNT)
	VMDin.stdin.write(CNTtools)
	VMDin.stdin.write(genNT)
	VMDin.stdin.write(pbcNT)
	VMDin.stdin.write(fixNT)
	VMDin.stdin.write(removeLangevinWater)

	# finished creating periodic nanotubes in VMD
	VMDin.stdin.flush()
	VMDin.stdin.close
	VMDin.communicate()
	if VMDin.returncode == 0:
		return tubePath, pbcPath


def solvate(inFile, N_0, S, n, m, force):
	""" The solvate module will create an armchair swcnt with N_0 rings and add N_0+S water molecules inside
	the nanotube, then write out new psf and pdb files for the nanotube. """
	tPath, pPath = tubeGen(inFile,inFile,N_0,n,m)

	os.remove(tPath + inFile + ".pdb")
	os.remove(tPath + inFile + ".psf")
	# Opens input nanotube psf and pdb files, and reads all the lines of each file into lists
	with open(pPath + inFile + ".psf") as psfFile:
		psfLines = psfFile.readlines()
	with open(pPath + inFile + ".pdb") as pdbFile:
		pdbLines = pdbFile.readlines()

	# Grabs the lengths of each of the lists
	lenPsf = len(psfLines)
	lenPdb = len(pdbLines)

	# String formats for the PSF and PDB file writing
	# Pdb 
	dampingCoeff = 0.00
	oxygen = "ATOM{0:>7}  OH2 TIP3            0.000   0.000{1:>8.3f}  0.00  0.00      TUB  O\n"
	hydro1 = "ATOM{0:>7}  H1  TIP3            0.000   0.766{1:>8.3f}  0.00  0.00      TUB  H\n"
	hydro2 = "ATOM{0:>7}  H2  TIP3            0.000  -0.766{1:>8.3f}  0.00  0.00      TUB  H\n"

	# Psf 
	opsf = "    {0:>5} TUB  {1:>5}  TIP3 OH2  OT    -0.834000       15.9994           0\n"
	h1psf = "    {0:>5} TUB  {1:>5}  TIP3 H1   HT     0.417000        1.0080           0\n"
	h2psf = "    {0:>5} TUB  {1:>5}  TIP3 H2   HT     0.417000        1.0080           0\n"

	# String format for the bonds and angles in the psf file
	sBondFormat = " {0: >8}{1: >8}{2: >8}{3: >8}{4: >8}{5: >8}{6: >8}{7: >8}\n"
	sAngleFormat = " {0: >8}{1: >8}{2: >8}{3: >8}{4: >8}{5: >8}{6: >8}{7: >8}{8: >8}\n"


	# Bonds lengths for different armchair nanotubes, now in Angstroms
	s0 = 1.418
	# Calculates the apothem of each regular hexagon in the nanotube and then calculates the distance
	## between the center of each ring in the tube
	#apothem = s*np.sqrt(3)/2
	l = float((N_0-0.75))*s0*np.sqrt(3)
	dist = s0*np.sqrt(3)

	# Initializing lists used below
	atoms = []
	bonds = []
	angles = []
	preAtoms = []	
	postAngles = []

	intBonds = []
	bondsFinal = []
	intAngles = []
	anglesFinal = []

	# Finds the original number of atoms in the Pdb file
	nAtoms = lenPdb-2
	# Calculates the new number of atoms after solvating
	newAtoms = nAtoms + (3*(N_0+S)) 
	# Calculates the new number of bonds and angles after solvating
	newBonds = int(nAtoms*(3./2)) + 2*(N_0+S)
	newAngles = nAtoms*3 + (N_0+S)

	# Iterates through all of the lines of the input Psf file, and records the index of the lines that contain
	## the string !NATOM, !NBOND, and !NTHETA, as wella as changes the line to update the new number of each
	for i in range(0,lenPsf):
		if "!NATOM" in psfLines[i]:
			psfLines[i] = "     {:3d} !NATOM\n".format( newAtoms )
			atomIndex = i
		elif "!NBOND" in psfLines[i]:
			psfLines[i] = "     {:3d} !NBOND: bonds\n".format( newBonds )
			bondIndex = i
		elif "!NTHETA" in psfLines[i]:
			psfLines[i] = "     {:3d} !NTHETA: angles\n".format( newAngles )
			angleIndex = i

	# Stores all of the original text lines that come before the atom section into a list
	for i in range(0,atomIndex):
		preAtoms.append(psfLines[i])

	# Stores the atoms into a list
	count = 1
	while psfLines[atomIndex+count].strip():
		atoms.append( psfLines[atomIndex+count] )
		count+=1
	# Stores the bonds into a list
	count = 1
	while psfLines[bondIndex+count].strip():
		bonds.append( psfLines[bondIndex+count] )
		count+=1
	# Stores the angles into a list
	count = 1
	while psfLines[angleIndex+count].strip():
		angles.append( psfLines[angleIndex+count] )
		count+=1
	# Stores all the text lines after the angles into a list
	for i in range(angleIndex+count,lenPsf):
		postAngles.append(psfLines[i])

	# Takes the bonds and angles in the original file and splits each line into individual numbers
	for bond in bonds:
		intBonds.append( bond.strip("\n").split() )
	for angle in angles:
		intAngles.append( angle.strip("\n").split() )

	# Compresses the list of lists into a single list of all of the angles and bonds in the original file
	intBonds = list(chain.from_iterable(intBonds))
	intAngles = list(chain.from_iterable(intAngles))	

	# Adds the new atoms to the original list of atoms
	for i in range(nAtoms+1, (3*(N_0+S)) + nAtoms+1, 3):
		atoms.append(opsf.format(i,i))
		atoms.append(h1psf.format(i+1,i+1))
		atoms.append(h2psf.format(i+2,i+2))

		intAngles.append(str(i+1))
		intAngles.append(str(i))
		intAngles.append(str(i+2))

		intBonds.append(str(i))
		intBonds.append(str(i+1))
		intBonds.append(str(i))
		intBonds.append(str(i+2))

	# Formats the list of bonds into the psf format with 8 columns
	for i in range(0,len(intBonds),8):
		try:
			bondsFinal.append( sBondFormat.format(intBonds[i],intBonds[i+1],intBonds[i+2],intBonds[i+3],
				intBonds[i+4], intBonds[i+5], intBonds[i+6], intBonds[i+7]) )

		except:
			diff = len(intBonds) - i
			tempStr = ""
			for j in range(i,i+diff):
				tempStr = tempStr + "{:>8}".format(intBonds[j])

			bondsFinal.append( " " + tempStr + "\n" )

	# Formates the list of angles into the psf format with 9 columns
	for i in range(0,len(intAngles),9):
		try:
			anglesFinal.append( sAngleFormat.format(intAngles[i],intAngles[i+1],intAngles[i+2],intAngles[i+3],
				intAngles[i+4], intAngles[i+5], intAngles[i+6], intAngles[i+7], intAngles[i+8]) )
	
		except:
			diff = len(intAngles) - i
			tempStr = ""
			for j in range(i,i+diff):
				tempStr = tempStr + "{:>8}".format(intAngles[j])

			anglesFinal.append( " " + tempStr + "\n" )

	sFactorAdjust = float(N_0) / (N_0 + S)

	for i in range(0,3*(N_0+S),3):
		if i==0:
			pdbLines[lenPdb-1] = oxygen.format(nAtoms+1, 0.000)
			pdbLines.append( hydro1.format(nAtoms+(i+2), 0.570) )
			pdbLines.append( hydro2.format(nAtoms+(i+3), 0.570) )
		
		else:
			pdbLines.append( oxygen.format(nAtoms+i+1, ((i/3)*dist*sFactorAdjust)) )
			pdbLines.append( hydro1.format(nAtoms+(i+2), 0.570+((i/3)*dist*sFactorAdjust)) )
			pdbLines.append( hydro2.format(nAtoms+(i+3), 0.570+((i/3)*dist*sFactorAdjust)) )

	# Writes the new pdb lines to a new pdb file
	pdbLines.append("END\n")
	pdbOut = open(pPath + inFile + "-solv.pdb",'w')
	pdbOut.writelines(pdbLines)
	pdbOut.close()

	# Writes the new psf lines to a new psf file
	psfOut = open(pPath + inFile + "-solv.psf",'w')
	psfOut.writelines(preAtoms)
	psfOut.writelines(psfLines[atomIndex])
	psfOut.writelines(atoms)
	psfOut.writelines("\n" + psfLines[bondIndex])
	psfOut.writelines(bondsFinal)
	psfOut.writelines("\n" + psfLines[angleIndex])
	psfOut.writelines(anglesFinal)
	psfOut.writelines(postAngles)
	psfOut.close()

	# Generates the file that indicates how strong the force is on each atom by modifying the occupancy column

	forceVal = "{0:.2f}".format(force)
	forceVal = forceVal[0:4]

	with open (pPath + inFile + "-solv.pdb") as f:
		flines = f.readlines()

	for i in range(1,len(flines)-1):
		if i < lenPdb-1:
			flines[i] = flines[i][0:56] + "0.00" + flines[i][60::]
		else:
			if "OH2" in flines[i]:
				flines[i] = flines[i][0:33] + "0.000   0.000   1.000  " + forceVal + flines[i][60::]
			else:
				flines[i] = flines[i][0:56] + "0.00" + flines[i][60::]


	outFile = open(pPath + inFile + "-force.pdb",'w')
	outFile.writelines(flines)
	outFile.close()

	# Generate a file that restrains each carbon atom to it's initial position with force constant k

	kVal = 3.0
	with open (pPath + inFile + "-solv.pdb") as kFile:
		kfLines = kFile.readlines()

	for i in range(1,len(kfLines)-1):
		if i < lenPdb-1:
			kfLines[i] = kfLines[i][0:56] + "{:.2f}".format(kVal) + kfLines[i][60::]
		else:
			kfLines[i] = kfLines[i][0:56] + "{:.2f}".format(0.00) + kfLines[i][60::]

	outkFile = open(pPath + inFile + "-restraint.pdb",'w')
	outkFile.writelines(kfLines)
	outkFile.close()

	return pPath


def simWrite(pbcFile, CNTpath, minimize, temp = 300, length = 20000, output = "waterSim"):
	""" simWrite generates a .conf file to use as input to namd2. To organize simulations
	with different parameters, simWrite will create a directory for the simulation using 
	the following template: ~/Simulations/cntl_nxm/temp/length/sim.dcd, sim.conf.  """

	simPath = CNTpath + str(temp) + "/" + str(length) + "/"
	if not os.path.exists(simPath):
		os.makedirs(simPath)

	# Grabs the CNT basis vectors
	x,y,z = getCNTBasis(CNTpath + pbcFile)

	# Read in lines of simulation file
	with open(basePath + "templates/sim_template.conf") as inFile:
		simLines = inFile.readlines()

	simLines[12] = "structure          " + quoted(CNTpath + pbcFile + "-solv.psf") + "\n"
	simLines[13] = "coordinates        " + quoted(CNTpath + pbcFile + "-solv.pdb") + "\n"

	simLines[15] = "set temperature    {:3d}\n".format(temp)
	simLines[30] = "cellBasisVector1    {0:<10.3f}{1:<10}{2:}\n".format(x,0.,0.)
	simLines[31] = "cellBasisVector2    {0:<10}{1:<10.3f}{2:}\n".format(0.,y,0.)
	simLines[32] = "cellBasisVector3    {0:<10}{1:<10}{2:.3f}\n".format(0.,0.,z)
	simLines[33] = "cellOrigin          {0:<10}{1:<10}{2:.3f}\n\n".format( 0,0, float(z)/2 )

	simLines[16] = "set outputname     " + quoted(simPath + output) + "\n"
	simLines[71] = "fixedAtomsFile      " + CNTpath + pbcFile + "-solv.pdb\n"
	simLines[75] = "consref 			" + CNTpath + pbcFile + "-restraint.pdb\n"
	simLines[76] = "conskfile			" + CNTpath + pbcFile + "-restraint.pdb\n"
	simLines[78] = "constraintScaling 	" + "{:.2f}".format(1.00)

	simLines[94] = "consforcefile 	    " + CNTpath + pbcFile + "-force.pdb\n"
	simLines[100] = "minimize " + str(minimize) + "\n"
	simLines[102] = "run {:5d} \n".format(length)

	# Write contents out to original file
	if os.path.exists(simPath):
		outFile = open(simPath + output + ".conf", "w")
		outFile.writelines(simLines)
		outFile.close()
		paramFile = basePath + "templates/par_all27_prot_lipid.prm"
		shutil.copy(paramFile, simPath)
		return simPath + output + ".conf"


def runSim(simPath):
	""" given an input path to a simulation file, runSim will call namd2 to run the simulation """
	Namd2in=subprocess.Popen(["namd2",simPath], stdin=subprocess.PIPE)
	Namd2in.stdin.flush()
	Namd2in.stdin.close
	Namd2in.communicate()
	if Namd2in.returncode==0:
		print("################################################################\nSimulation finished\nSimulation file saved into " + simPath.replace(".conf",".dcd") + ".\n################################################################\n")


# find cell basis
def getCNTBasis(CNT):
	""" getCNTBasis finds the basis of a nanotube with filename outFile. """
	# Opens the CNT prebond file and reads the header of the file
	with open(CNT + "-prebond.pdb") as basisFile:
		header = basisFile.next()

	# Splits the first line of the CNT-prebond file, and finds the x,y,z basis vectors of the CNT 
	basis = re.split('\s+',header)
	xVec = eval(basis[1])
	yVec = eval(basis[2])
	zVec = eval(basis[3])
	
	return xVec, yVec, zVec


def main(FNAME, N_0, S, n, m, TEMP, LENGTH, FORCESTRENGTH, minimize):
	pbcPath = solvate(FNAME, N_0, S, n, m, FORCESTRENGTH)
	simPath = simWrite(FNAME, pbcPath, minimize, TEMP, LENGTH, FNAME)
	runSim(simPath)

