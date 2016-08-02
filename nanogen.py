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
import datetime
import socket

import myutil as my

# Maybe adding later... need to think more about the implementation
# class Nanotube:

#     def __init__(self, rings, n, m):
#         self.rings = rings
#         self.n = n
#         self.m = m

def cntWrite(paths, fileName, N0, n, m):
    """ nanoWrite creates a periodic nanotube with the following input parameters:
    fileName is the name of the initial nanotube with number of rings N0 and
    dimensions n x m."""

    # Bonds lengths of different armchair nanotubes in nanometers
    s0 = 0.1418
    # calculates the length of the nanotube based on bond lengths
    l = float(N0-0.75) * s0 * np.sqrt(3)

    if os.path.isfile(my.configFile(paths['pbc'], fileName + '.psf')) \
        and os.path.isfile(my.configFile(paths['pbc'], fileName + '.pdb')):

        
        print("################################################################\n" \
            + "Configuration files already exist in:\n" + my.config(paths['pbc']) + ".\n" \
            + "################################################################\n")

        return

    print("################################################################\n" \
        + "Running VMD to generate nanotube configuration files.\n" \
        + "################################################################\n")

    if not os.path.exists(my.config(paths['N0'])):
        os.makedirs(my.config(paths['N0']))
    if not os.path.exists(my.config(paths['pbc'])):
        os.makedirs(my.config(paths['pbc']))
        
    # Create the nanotube in VMD

    logFile = open(my.config(paths['N0']) + fileName + ".log", "w")

    # Opening a pipe to VMD in the shell
    VMDin=subprocess.Popen(["vmd", "-dispdev", "none"], stdin=subprocess.PIPE, \
                       stdout=logFile)

    sourceCNT = "source CNTtools.tcl\n"
    CNTtools = "package require CNTtools 1.0\n"

    genNT = "genNT " + my.quoted(fileName) + " " + my.quoted(my.config(paths['N0'])) + " " \
        + str(l) + " " + str(n) + " " + str(m) + "\n"

    # run commands through pipe and saves to file
    VMDin.stdin.write(sourceCNT)
    VMDin.stdin.write(CNTtools)
    VMDin.stdin.write(genNT)
    VMDin.stdin.flush()
    VMDin.stdin.close
    VMDin.communicate()

    # Make the nanotube periodic in VMD

    logFile = open(my.config(paths['pbc']) + fileName + ".log", "w")

    # Opening a pipe to VMD in the shell
    VMDin=subprocess.Popen(["vmd", "-dispdev", "none"], stdin=subprocess.PIPE, \
                       stdout=logFile)

    sourceCNT = "source CNTtools.tcl\n"
    CNTtools = "package require CNTtools 1.0\n"
    pbcNT = "pbcNT " + my.quoted(my.config(paths['N0']) + fileName) + " " \
        + my.quoted(my.config(paths['pbc']) + fileName) + " default\n"

    # run commands through pipe and saves to file
    VMDin.stdin.write(sourceCNT)
    VMDin.stdin.write(CNTtools)
    VMDin.stdin.write(pbcNT)
    VMDin.stdin.flush()
    VMDin.stdin.close
    VMDin.communicate()
    
    # Delete the now-uneeded -prebond pdb and psf files
    os.remove(my.config(paths['pbc']) + fileName + "-prebond.pdb")
    os.remove(my.config(paths['pbc']) + fileName + "-prebond.psf")
    
    # Read back in the pdb file and make sure the z-coordinate spacing is regular
    
    alignCNT(my.config(paths['pbc']) + fileName + ".pdb")
    
def alignCNT(pdbFileName, z0 = 0.614):

    with open(pdbFileName) as pdbFile:
        pdbLines = pdbFile.readlines()

    lenPDB = len(pdbLines)
    
    for i in range(0, lenPDB):
        
        if "CRYST" in pdbLines[i]:
        
            L = float(pdbLines[i][24:33])
            L = z0 * np.round(L/z0)
            pdbLines[i] = pdbLines[i][:24] + "{:9.3f}".format(L) + pdbLines[i][33:]
            
        elif "ATOM" in pdbLines[i]:

            z = float(pdbLines[i][46:54])
            z = z0 * np.round(z/z0)
            pdbLines[i] = pdbLines[i][:46] + "{:8.3f}".format(z) + pdbLines[i][54:]
            
    pdbFile.close()
    
    pdbOut = open(pdbFileName, 'w')
    pdbOut.writelines(pdbLines)
    pdbOut.close()
    
def waterWrite(paths, fileName, N0, S):
    """ waterWrite adds N0 + S water molecules to the inside of the nanotube,
        then write out new psf and pdb files for the nanotube. """

    if os.path.isfile(my.configFile(paths['solvate'], fileName + '.psf')) \
        and os.path.isfile(my.configFile(paths['solvate'], fileName + '.pdb')):
    
        print("################################################################\n" \
            + "Configuration files already exist in:\n" + my.config(paths['solvate']) + ".\n" \
            + "################################################################\n")
        
        return
        
    print("################################################################\n" \
        + "Adding water to the CNT and writing configuration files.\n" \
        + "################################################################\n")

    if not os.path.exists(my.config(paths['solvate'])):
        os.makedirs(my.config(paths['solvate']))

    # Opens input nanotube psf and pdb files, and reads all the lines of each file into lists
    with open(my.configFile(paths['pbc'], fileName + '.psf')) as psfFile:
        psfLines = psfFile.readlines()
    with open(my.configFile(paths['pbc'], fileName + '.pdb')) as pdbFile:
        pdbLines = pdbFile.readlines()

    # Grabs the lengths of each of the lists
    lenPsf = len(psfLines)
    lenPdb = len(pdbLines)

    # String formats for the PSF and PDB file writing
    # Pdb
    dampingCoeff = 0.00
    oxygen = "ATOM{0:>7}  OH2 TIP3 {1:4.0f}       0.000   0.000{2:>8.3f}  0.00  0.00      WTR  O\n"
    hydro1 = "ATOM{0:>7}  H1  TIP3 {1:4.0f}       0.000   0.766{2:>8.3f}  0.00  0.00      WTR  H\n"
    hydro2 = "ATOM{0:>7}  H2  TIP3 {1:4.0f}       0.000  -0.766{2:>8.3f}  0.00  0.00      WTR  H\n"

    # Psf
    # opsf  = "   {0:>5} TUB {1:>5} TIP3 OH2  OT    -0.834000       15.9994           0\n"
    # h1psf = "   {0:>5} TUB {1:>5} TIP3 H1   HT     0.417000        1.0080           0\n"
    # h2psf = "   {0:>5} TUB {1:>5} TIP3 H2   HT     0.417000        1.0080           0\n"

    opsf  = "   {0:>5} WTR  {1:<4} TIP3 OH2  OT    -0.834000       15.9994           0\n"
    h1psf = "   {0:>5} WTR  {1:<4} TIP3 H1   HT     0.417000        1.0080           0\n"
    h2psf = "   {0:>5} WTR  {1:<4} TIP3 H2   HT     0.417000        1.0080           0\n"

    # String format for the bonds and angles in the psf file
    sBondFormat  = " {0: >8}{1: >8}{2: >8}{3: >8}{4: >8}{5: >8}{6: >8}{7: >8}\n"
    sAngleFormat = " {0: >8}{1: >8}{2: >8}{3: >8}{4: >8}{5: >8}{6: >8}{7: >8}{8: >8}\n"

    # Bonds lengths for different armchair nanotubes, now in Angstroms
    s0 = 1.418
    
    # Calculates the apothem of each regular hexagon in the nanotube and then
    # calculates the distance between the center of each ring in the tube
    #apothem = s*np.sqrt(3)/2
    l = float(N0-0.75) * s0 * np.sqrt(3)
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
    newAtoms = nAtoms + (3*(N0+S)) 
    # Calculates the new number of bonds and angles after solvating
    newBonds = int(nAtoms*(3./2)) + 2*(N0+S)
    newAngles = nAtoms*3 + (N0+S)

    # Iterates through all of the lines of the input Psf file, and records the
    # index of the lines that contain the string !NATOM, !NBOND, and !NTHETA,
    # as well as changes the line to update the new number of each
    for i in range(0, lenPsf):
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
    for i in range(0, atomIndex):
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
    for i in range(angleIndex+count, lenPsf):
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
    for i in range(nAtoms+1, (3*(N0+S)) + nAtoms+1, 3):
        atoms.append(opsf.format(i, int((i-nAtoms)/3)+2))
        atoms.append(h1psf.format(i+1, int((i-nAtoms)/3)+2))
        atoms.append(h2psf.format(i+2, int((i-nAtoms)/3)+2))

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
            bondsFinal.append( sBondFormat.format(intBonds[i], intBonds[i+1],
                intBonds[i+2], intBonds[i+3], intBonds[i+4], intBonds[i+5],
                intBonds[i+6], intBonds[i+7]) )

        except:
            diff = len(intBonds) - i
            tempStr = ""
            for j in range(i, i+diff):
                tempStr = tempStr + "{:>8}".format(intBonds[j])

            bondsFinal.append( " " + tempStr + "\n" )

    # Formates the list of angles into the psf format with 9 columns
    for i in range(0, len(intAngles), 9):
        try:
            anglesFinal.append( sAngleFormat.format(intAngles[i], intAngles[i+1],
                intAngles[i+2], intAngles[i+3], intAngles[i+4], intAngles[i+5],
                intAngles[i+6], intAngles[i+7], intAngles[i+8]) )
    
        except:
            diff = len(intAngles) - i
            tempStr = ""
            for j in range(i, i+diff):
                tempStr = tempStr + "{:>8}".format(intAngles[j])

            anglesFinal.append( " " + tempStr + "\n" )

    sFactorAdjust = float(N0) / (N0 + S)

    for i in range(0,3*(N0+S), 3):
        oxZ = dist/4 + ((i/3)*dist*sFactorAdjust)
        hyZ = oxZ + 0.570
        if i==0:
            pdbLines[lenPdb-1] = oxygen.format(nAtoms+1, i/3+2, oxZ)
            pdbLines.append( hydro1.format(nAtoms+(i+2), i/3+2, hyZ) )
            pdbLines.append( hydro2.format(nAtoms+(i+3), i/3+2, hyZ) )
        
        else:
            pdbLines.append( oxygen.format(nAtoms+i+1, i/3+2, oxZ) )
            pdbLines.append( hydro1.format(nAtoms+(i+2), i/3+2, hyZ) )
            pdbLines.append( hydro2.format(nAtoms+(i+3), i/3+2, hyZ) )

    # Writes the new pdb lines to a new pdb file
    pdbLines.append("END\n")
    pdbOut = open(my.configFile(paths['solvate'], fileName + '.pdb'), 'w')
    pdbOut.writelines(pdbLines)
    pdbOut.close()

    # Writes the new psf lines to a new psf file
    psfOut = open(my.configFile(paths['solvate'], fileName + '.psf'), 'w')
    psfOut.writelines(preAtoms)
    psfOut.writelines(psfLines[atomIndex])
    psfOut.writelines(atoms)
    psfOut.writelines("\n" + psfLines[bondIndex])
    psfOut.writelines(bondsFinal)
    psfOut.writelines("\n" + psfLines[angleIndex])
    psfOut.writelines(anglesFinal)
    psfOut.writelines(postAngles)
    psfOut.close()

    print("################################################################\n" \
        + "Post-processing the CNT and rewriting configuration files.\n" \
        + "################################################################\n")

    logFile = open(my.config(paths['solvate']) + fileName + ".log", "w")

    # Opening a pipe to VMD in the shell
    VMDin=subprocess.Popen(["vmd", "-dispdev", "none"], stdin=subprocess.PIPE, \
                           stdout=logFile)

    sourceCNT = "source CNTtools.tcl\n"
    CNTtools = "package require CNTtools 1.0\n"

    centerNT = "centerNT " + my.quoted(my.config(paths['solvate']) + fileName) + "\n"

    VMDin.stdin.write(sourceCNT)
    VMDin.stdin.write(CNTtools)
    VMDin.stdin.write(centerNT)
    VMDin.stdin.flush()
    VMDin.stdin.close
    VMDin.communicate()
    
def pdbWrite(paths, fileName, name, value):
    
    if os.path.isfile(my.configFile(paths[name], fileName + '.pdb')):
        
        print("################################################################\n" \
              + "Configuration files already exist in:\n" + my.config(paths[name]) + ".\n" \
              + "################################################################\n")
              
        return

    print("################################################################\n" \
          + "Writing configuration file in:\n" + my.config(paths[name]) + ".\n" \
          + "################################################################\n")
    
    
    if not os.path.exists(my.config(paths[name])):
        os.makedirs(my.config(paths[name]))

    # Opening a pipe to VMD in the shell
    logFile = open(my.config(paths[name]) + fileName + ".log", "w")
    VMDin=subprocess.Popen(["vmd", "-dispdev", "none"], stdin=subprocess.PIPE, \
                           stdout=logFile)

    sourceCNT = "source CNTtools.tcl\n"
    CNTtools = "package require CNTtools 1.0\n"

    tclName = "NT" + name + " " + my.quoted(my.config(paths['solvate']) + fileName) + " " \
        + my.quoted(my.config(paths[name]) + fileName) + " " + str(value) + "\n"

    VMDin.stdin.write(sourceCNT)
    VMDin.stdin.write(CNTtools)
    VMDin.stdin.write(tclName)
    VMDin.stdin.flush()
    VMDin.stdin.close
    VMDin.communicate()

def confWrite(paths, fileName, force, restraint, temp, thermostat, duration, minDuration, dt, \
        outputFreq, PME, hostname, extend):
    """ confWrite generates a .conf file to use as input to namd2. To organize simulations
    with different parameters, confWrite will create a directory for the simulation."""

    confFile = paths['run'] + fileName + '.conf'
    wisdomFile = "FFTW_NAMD_2.11_" + hostname + ".txt"
    paramFile = "par_all27_prot_lipid.prm"

    print("################################################################\n" \
          + "Writing NAMD configuration file " + confFile + ".\n" \
          + "################################################################\n")
      
    if not os.path.exists(paths['run']):
        os.makedirs(paths['run'])

    # Grabs the CNT basis vectors
    x, y, z = getCNTBasis(my.config(paths['pbc']) + fileName + ".pdb")

    outFile = open(confFile, "w")
    
    outFile.write("# Configuration file written on host " + hostname + "\n" \
        + "# on " + str(datetime.datetime.now().date()) \
        + " at " + str(datetime.datetime.now().time()).partition('.')[0] + ".\n\n")
        
    outFile.write("# Limit the length of the log file\n\n")
    
    outFile.write("{0:<20}".format("outputEnergies") + str(100*outputFreq) + "\n\n")
    
    outFile.write("# Set up periodic boundary conditions\n\n")
    
    outFile.write("cellBasisVector1    {0:>10.3f}{1:>10.3f}{2:>10.3f}\n".format(x, 0., 0.))
    outFile.write("cellBasisVector2    {0:>10.3f}{1:>10.3f}{2:>10.3f}\n".format(0., y, 0.))
    outFile.write("cellBasisVector3    {0:>10.3f}{1:>10.3f}{2:>10.3f}\n".format(0., 0., z))
    outFile.write("cellOrigin          {0:>10}{1:>10}{2:>10}\n\n".format(0, 0, 0))
    
    outFile.write("{0:<20}".format("wrapWater") + "off" + "\n")
    outFile.write("{0:<20}".format("wrapAll") + "off" + "\n")
    outFile.write("\n")

    outFile.write("# Set the input and output files\n\n")
    
    outFile.write("{0:<20}".format("structure") + fileName + ".psf" + "\n")
    outFile.write("{0:<20}".format("coordinates") + fileName + ".pdb" + "\n")
    outFile.write("{0:<20}".format("outputname") + fileName + "\n")
    outFile.write("\n")
    
    outFile.write("# Set the force field parameters\n\n")
    
    outFile.write("{0:<20}".format("paraTypeCharmm") + "on" + "\n")
    outFile.write("{0:<20}".format("parameters") + "par_all27_prot_lipid.prm" + "\n") 
    outFile.write("{0:<20}".format("exclude") + "scaled1-4" + "\n")
    outFile.write("{0:<20}".format("cutoff") + "12.0" + "\n")
    outFile.write("{0:<20}".format("pairlistdist") + "14.0" + "\n")
    outFile.write("{0:<20}".format("switching") + "on" + "\n")
    outFile.write("{0:<20}".format("switchdist") + "10.0" + "\n")
    if PME == 'on':
        outFile.write("{0:<20}".format("PME") + "yes" + "\n")
        outFile.write("{0:<20}".format("PMEGridSpacing") + "1.0" + "\n")
        outFile.write("{0:<20}".format("FFTWWisdomFile") + wisdomFile + "\n")
    else:
        outFile.write("{0:<20}".format("PME") + "no" + "\n")
    outFile.write("\n")

    outFile.write("# Set the integration parameters\n\n")
    outFile.write("{0:<20}".format("timestep") + str(dt) + "\n")
    outFile.write("{0:<20}".format("nonbondedFreq") + "2" + "\n")
    outFile.write("{0:<20}".format("fullElectFrequency") + "4" + "\n")
    outFile.write("{0:<20}".format("stepspercycle") + "20" + "\n")
    outFile.write("{0:<20}".format("rigidBonds") + "water" + "\n")

    outFile.write("\n")

    if not extend:
        outFile.write("{0:<20}".format("restartfreq") + str(outputFreq) + "\n")
        outFile.write("{0:<20}".format("dcdfreq") + str(outputFreq) + "\n")
        outFile.write("{0:<20}".format("veldcdfreq") + str(outputFreq) + "\n")
        outFile.write("{0:<20}".format("forcedcdfreq") + str(outputFreq) + "\n")
    else:
        outFile.write("{0:<20}".format("restartfreq") + str(100) + "\n")
        outFile.write("{0:<20}".format("dcdfreq") + str(100) + "\n")
        outFile.write("{0:<20}".format("veldcdfreq") + str(100) + "\n")
        outFile.write("{0:<20}".format("forcedcdfreq") + str(100) + "\n")
    
    outFile.write("\n")

    outFile.write("# Set the restraints on the carbon\n\n")
    
    if restraint == 0:
        outFile.write("{0:<20}".format("fixedAtoms") + "on" + "\n")
        outFile.write("{0:<20}".format("fixedAtomsFile") + fileName + "-restraint.pdb" + "\n")
        outFile.write("{0:<20}".format("fixedAtomsCol") + "B" + "\n")
    else:
        outFile.write("{0:<20}".format("constraints") + "on" + "\n")
        outFile.write("{0:<20}".format("consref") + fileName + "-restraint.pdb" + "\n")
        outFile.write("{0:<20}".format("conskfile") + fileName + "-restraint.pdb" + "\n")
        outFile.write("{0:<20}".format("conskcol") + "O" + "\n")

    outFile.write("\n")
    
    outFile.write("# Set the external forcing\n\n")
    
    if not force == 0:
        outFile.write("{0:<20}".format("constantForce") + "yes" + "\n")
    else:
        outFile.write("{0:<20}".format("constantForce") + "no" + "\n")
    outFile.write("{0:<20}".format("consForceFile") + fileName + "-forcing.pdb" + "\n")
    outFile.write("{0:<20}".format("consForceScaling") + "0.014392621" + "\n")
    outFile.write("\n")

    outFile.write("# Set the Langevin thermostat\n\n")
    
    if restraint == 0 or (not thermostat):
        outFile.write("{0:<20}".format("langevin") + "off" + "\n")
    else:
        outFile.write("{0:<20}".format("langevin") + "on" + "\n")
        outFile.write("{0:<20}".format("langevinFile") + fileName + "-langevin.pdb" + "\n")
        outFile.write("{0:<20}".format("langevinCol") + "O" + "\n")
        outFile.write("{0:<20}".format("langevinTemp") + str(temp) + "\n")
        
    outFile.write("\n")

    outFile.write("# Set the execution parameters\n\n")

    if not extend:
        outFile.write("{0:<20}".format("temperature") + str(temp) + "\n")
        outFile.write("{0:<20}".format("minimize") + str(minDuration) + "\n")
        outFile.write("{0:<20}".format("reinitvels") + str(temp) + "\n")
        outFile.write("{0:<20}".format("run") + str(duration) + "\n")
    else:
        outFile.write("{0:<20}".format("bincoordinates") + fileName + ".restart.coor" + "\n")
        outFile.write("{0:<20}".format("binvelocities") + fileName + ".restart.vel" + "\n")
        outFile.write("{0:<20}".format("run") + str(100000) + "\n")

    outFile.close()
    
    # Make a link to the parameter file in the simulation run folder
    if not os.path.exists(paths['run'] + paramFile):
        os.link("Parameter Files/" + paramFile, paths['run'] + paramFile)
        
    # Make a link to the wisdom file in the simulation run folder
    if (not os.path.exists(paths['run'] + wisdomFile)) \
        and os.path.exists("Parameter Files/" + wisdomFile):
        os.link("Parameter Files/" + wisdomFile, paths['run'] + wisdomFile)
    
    # Make links to all the other configuration files in the simulation run folder
    if not os.path.exists(paths['run'] + fileName + '.pdb'):
        os.link(my.configFile(paths['solvate'], fileName + '.pdb'), paths['run'] + fileName + '.pdb')
    if not os.path.exists(paths['run'] + fileName + '.psf'):
        os.link(my.configFile(paths['solvate'], fileName + '.psf'), paths['run'] + fileName + '.psf')
    if not os.path.exists(paths['run'] + fileName + '-restraint.pdb'):
        os.link(my.configFile(paths['restraint'], fileName + '.pdb'), paths['run'] + fileName + '-restraint.pdb')
    if not os.path.exists(paths['run'] + fileName + '-forcing.pdb'):
        os.link(my.configFile(paths['forcing'], fileName + '.pdb'), paths['run'] + fileName + '-forcing.pdb')
    if not os.path.exists(paths['run'] + fileName + '-langevin.pdb') and restraint != 0:
        os.link(my.configFile(paths['temperature'], fileName + '.pdb'), paths['run'] + fileName + '-langevin.pdb')
    
    if extend:
        if not os.path.exists(paths['run'] + fileName + ".restart.coor"):
            os.link(paths['parent data'] + fileName + ".restart.coor", paths['run'] + fileName + ".restart.coor")
        if not os.path.exists(paths['run'] + fileName + ".restart.vel"):
            os.link(paths['parent data'] + fileName + ".restart.vel", paths['run'] + fileName + ".restart.vel")

    return confFile
    
def runSim(simPath, simFile, output, hostname):
    """ given an input path to a simulation file, runSim will call namd2 to run the simulation """
    
    logFile = open(simPath + output + ".log", "w")
    
    print("################################################################\n" \
          + "Starting simulation " + "on " + str(datetime.datetime.now().date()) \
          + " at " + str(datetime.datetime.now().time()).partition('.')[0] + ".\n" \
          + "################################################################\n")
    
    startTime = time.time()
    
    numThreads = {}
    numThreads['cascade'] = '+p4'
    numThreads['Home-iMac'] = '+p4'
    numThreads['nanotube'] = '+p2'
    numThreads['CGG-MacBook-Air'] = '+p2'

    Namd2in=subprocess.Popen(['namd2', numThreads[hostname], simFile], stdin=subprocess.PIPE, stdout=logFile, stderr=logFile)
    Namd2in.stdin.flush()
    Namd2in.stdin.close
    Namd2in.communicate()
    
    elapsedTime = time.time() - startTime
    
    if Namd2in.returncode == 0:
        print("################################################################\n" \
              + "Simulation finished\nSimulation file saved into " \
              + simPath.replace(".conf", ".dcd") + ".\n" \
              + "Elapsed time: " + "{:.0f}".format(elapsedTime) + " s = " \
              + "{:.1f}".format(elapsedTime/60) + " m = " \
              + "{:.1f}".format(elapsedTime/3600) + " h\n" \
              + "################################################################\n")
        runTime = elapsedTime/3600
    else:
        print("################################################################\n" \
              + "Simulation FAIL.\n" \
              + "################################################################\n")
        runTime = 0
              
    return Namd2in.returncode, runTime

# find cell basis
def getCNTBasis(inFile):
    """ getCNTBasis finds the basis of a nanotube with filename outFile. """
    # Opens the CNT prebond file and reads the header of the file
    with open(inFile) as basisFile:
        header = basisFile.next()

    # Splits the first line of the CNT-prebond file, and finds the x, y, z basis vectors of the CNT 
    basis = re.split('\s+', header)
    xVec = eval(basis[1])
    yVec = eval(basis[2])
    zVec = eval(basis[3])
    
    return xVec, yVec, zVec

def run(paramDict):
    
    # The paths dictionary holds the path name to each tier of the data hierarchy
    paths = my.makePaths(paramDict)
    
    if os.path.exists(paths['data']) and paramDict['Run Type'] == 'New':
        print("################################################################\n" \
              + "Simulation data folder already exists ... aborting.\n" \
              + "################################################################\n")
        return
    elif paramDict['Run Type'] == 'Extend' and (not os.path.exists(paths['parent data'])):
        print("################################################################\n" \
              + "No simulation data folder to extend run from ... aborting.\n" \
              + "################################################################\n")
        return
            
    hostname = socket.gethostname().partition('.')[0]
    
    if paramDict['Run Type'] == 'New':
        paths['run'] = os.getenv('HOME') + '/' + paramDict['File Name'] + '/'
    else:
        paths['run'] = os.getenv('HOME') + '/' + 'Extend-1' + '/'
    
    if paramDict['Run Type'] == 'New':
    
        cntWrite(paths, paramDict['File Name'], paramDict['N0'], paramDict['n'], paramDict['m'])
        waterWrite(paths, paramDict['File Name'], paramDict['N0'], paramDict['S'])
        pdbWrite(paths, paramDict['File Name'], 'restraint', paramDict['Restraint'])
        pdbWrite(paths, paramDict['File Name'], 'forcing', paramDict['Force (pN)'])
        pdbWrite(paths, paramDict['File Name'], 'temperature', paramDict['Damping'])
        
        if paramDict['Thermostat'] == 'On':
            therm = True
        else:
            therm = False
        confFile = confWrite(paths, paramDict['File Name'], paramDict['Force (pN)'], \
            paramDict['Restraint'], paramDict['Temperature (K)'], therm, paramDict['Duration'], \
            paramDict['Min Duration'], paramDict['dt (fs)'], \
            paramDict['outputFreq'], paramDict['PME'], hostname, False)
            
    else:
    
        confFile = confWrite(paths, paramDict['File Name'], paramDict['Force (pN)'], \
            paramDict['Restraint'], paramDict['Temperature (K)'], paramDict['Duration'], \
            paramDict['Min Duration'], paramDict['dt (fs)'], \
            paramDict['outputFreq'], paramDict['PME'], hostname, True)
    
    result, runTime = runSim(paths['run'], confFile, paramDict['File Name'], hostname)
    
    if result == 0:
    
        if not os.path.isfile(hostname + '.csv'):
        
            outFile = open(hostname + '.csv', "w")
            outFile.write("Date,Run Time (h),File Name,Folder,Run Type,N0,S,n,m," \
                + "Temperature (K),Damping,Thermostat," \
                + "Force,PME,Restraint,Duration,Min Duration,dt (fs),outputFreq\n")
                
        outFile = open(hostname + '.csv', "a")
        outFile.write(str(datetime.date.today()) + ',' \
            + "{:.1f}".format(runTime) + ',' \
            + paramDict['File Name'] + ',' \
            + paramDict['Folder'] + ',' \
            + paramDict['Run Type'] + ',' \
            + str(paramDict['N0']) + ',' \
            + str(paramDict['S']) + ',' \
            + str(paramDict['n']) + ',' \
            + str(paramDict['m']) + ',' \
            + str(paramDict['Temperature (K)']) + ',' \
            + str(paramDict['Damping']) + ',' \
            + paramDict['Thermostat'] + ',' \
            + str(paramDict['Force (pN)']) + ',' \
            + paramDict['PME'] + ','
            + str(paramDict['Restraint']) + ',' \
            + str(paramDict['Duration']) + ',' \
            + str(paramDict['Min Duration']) + ',' \
            + str(paramDict['dt (fs)']) + ',' \
            + str(paramDict['outputFreq']) + '\n')
        
    shutil.move(paths['run'], paths['data'])
