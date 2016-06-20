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

# Maybe adding later... need to think more about the implementation
# class Nanotube:

#     def __init__(self, rings, n, m):
#         self.rings = rings
#         self.n = n
#         self.m = m

# Simple function to create a quoted string to protect file/path names with spaces
def quoted(theString):
    """quoted puts double quotes around any string"""

    return '"' + theString + '"'

def config(thePath):
    """config turns a path into the appropriate path for config files"""

    return thePath + 'Config Files/'

def configFile(thePath, theFile):
    """configFile takes a path and a file name and constructs the full path"""

    return config(thePath) + theFile

def cntWrite(paths, fileName, N0, n, m):
    """ nanoWrite creates a periodic nanotube with the following input parameters:
    fileName is the name of the initial nanotube with number of rings N0 and
    dimensions n x m."""

    # Bonds lengths of different armchair nanotubes in nanometers
    s0 = 0.1418
    # calculates the length of the nanotube based on bond lengths
    l = float(N0-0.75) * s0 * np.sqrt(3)

    if os.path.isfile(configFile(paths['pbc'], fileName + '.psf')) \
        and os.path.isfile(configFile(paths['pbc'], fileName + '.pdb')):

        
        print("################################################################\n" \
            + "Configuration files already exist in:\n" + config(paths['pbc']) + ".\n" \
            + "################################################################\n")

        return

    print("################################################################\n" \
        + "Running VMD to generate nanotube configuration files.\n" \
        + "################################################################\n")

    if not os.path.exists(config(paths['N0'])):
        os.makedirs(config(paths['N0']))
    if not os.path.exists(config(paths['pbc'])):
        os.makedirs(config(paths['pbc']))
        
    # Create the nanotube in VMD

    logFile = open(config(paths['N0']) + fileName + ".log", "w")

    # Opening a pipe to VMD in the shell
    VMDin=subprocess.Popen(["vmd", "-dispdev", "none"], stdin=subprocess.PIPE, \
                       stdout=logFile)

    sourceCNT = "source CNTtools.tcl\n"
    CNTtools = "package require CNTtools 1.0\n"

    genNT = "genNT " + quoted(fileName) + " " + quoted(config(paths['N0'])) + " " \
        + str(l) + " " + str(n) + " " + str(m) + "\n"

    # run commands through pipe and saves to file
    VMDin.stdin.write(sourceCNT)
    VMDin.stdin.write(CNTtools)
    VMDin.stdin.write(genNT)
    VMDin.stdin.flush()
    VMDin.stdin.close
    VMDin.communicate()

    # Make the nanotube periodic in VMD

    logFile = open(config(paths['pbc']) + fileName + ".log", "w")

    # Opening a pipe to VMD in the shell
    VMDin=subprocess.Popen(["vmd", "-dispdev", "none"], stdin=subprocess.PIPE, \
                       stdout=logFile)

    sourceCNT = "source CNTtools.tcl\n"
    CNTtools = "package require CNTtools 1.0\n"
    pbcNT = "pbcNT " + quoted(config(paths['N0']) + fileName) + " " \
        + quoted(config(paths['pbc']) + fileName) + " default\n"

    # run commands through pipe and saves to file
    VMDin.stdin.write(sourceCNT)
    VMDin.stdin.write(CNTtools)
    VMDin.stdin.write(pbcNT)
    VMDin.stdin.flush()
    VMDin.stdin.close
    VMDin.communicate()
        
def waterWrite(paths, fileName, N0, S):
    """ waterWrite adds N0 + S water molecules to the inside of the nanotube,
        then write out new psf and pdb files for the nanotube. """

    if os.path.isfile(configFile(paths['solvate'], fileName + '.psf')) \
        and os.path.isfile(configFile(paths['solvate'], fileName + '.pdb')):
    
        print("################################################################\n" \
            + "Configuration files already exist in:\n" + config(paths['solvate']) + ".\n" \
            + "################################################################\n")
        
        return
        
    print("################################################################\n" \
        + "Adding water to the CNT and writing configuration files.\n" \
        + "################################################################\n")

    if not os.path.exists(config(paths['solvate'])):
        os.makedirs(config(paths['solvate']))

    # Opens input nanotube psf and pdb files, and reads all the lines of each file into lists
    with open(configFile(paths['pbc'], fileName + '.psf')) as psfFile:
        psfLines = psfFile.readlines()
    with open(configFile(paths['pbc'], fileName + '.pdb')) as pdbFile:
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
    # as wella as changes the line to update the new number of each
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
        if i==0:
            pdbLines[lenPdb-1] = oxygen.format(nAtoms+1, i/3+2, 0.000)
            pdbLines.append( hydro1.format(nAtoms+(i+2), i/3+2, 0.570) )
            pdbLines.append( hydro2.format(nAtoms+(i+3), i/3+2, 0.570) )
        
        else:
            pdbLines.append( oxygen.format(nAtoms+i+1, i/3+2, ((i/3)*dist*sFactorAdjust)) )
            pdbLines.append( hydro1.format(nAtoms+(i+2), i/3+2, 0.570+((i/3)*dist*sFactorAdjust)) )
            pdbLines.append( hydro2.format(nAtoms+(i+3), i/3+2, 0.570+((i/3)*dist*sFactorAdjust)) )

    # Writes the new pdb lines to a new pdb file
    pdbLines.append("END\n")
    pdbOut = open(configFile(paths['solvate'], fileName + '.pdb'), 'w')
    pdbOut.writelines(pdbLines)
    pdbOut.close()

    # Writes the new psf lines to a new psf file
    psfOut = open(configFile(paths['solvate'], fileName + '.psf'), 'w')
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

    logFile = open(config(paths['solvate']) + fileName + ".log", "w")

    # Opening a pipe to VMD in the shell
    VMDin=subprocess.Popen(["vmd", "-dispdev", "none"], stdin=subprocess.PIPE, \
                           stdout=logFile)

    sourceCNT = "source CNTtools.tcl\n"
    CNTtools = "package require CNTtools 1.0\n"

    fixNT = "fixNT " + quoted(config(paths['pbc']) + fileName) + "\n"

    removeLangevinWater = "removeLangevinWater " + quoted(config(paths['pbc']) + fileName) + "\n"

    VMDin.stdin.write(sourceCNT)
    VMDin.stdin.write(CNTtools)
    # VMDin.stdin.write(fixNT)
    # VMDin.stdin.write(removeLangevinWater)

    # finished creating periodic nanotubes in VMD
    VMDin.stdin.flush()
    VMDin.stdin.close
    VMDin.communicate()
    
def pdbWrite(paths, fileName, name, value):
    
    if os.path.isfile(configFile(paths[name], fileName + '.pdb')):
        
        print("################################################################\n" \
              + "Configuration files already exist in:\n" + config(paths[name]) + ".\n" \
              + "################################################################\n")
              
        return

    print("################################################################\n" \
          + "Writing configuration file in:\n" + config(paths[name]) + ".\n" \
          + "################################################################\n")
    
    
    if not os.path.exists(config(paths[name])):
        os.makedirs(config(paths[name]))

    # Opening a pipe to VMD in the shell
    logFile = open(config(paths[name]) + fileName + ".log", "w")
    VMDin=subprocess.Popen(["vmd", "-dispdev", "none"], stdin=subprocess.PIPE, \
                           stdout=logFile)

    sourceCNT = "source CNTtools.tcl\n"
    CNTtools = "package require CNTtools 1.0\n"

    tclName = "NT" + name + " " + quoted(config(paths['solvate']) + fileName) + " " \
        + quoted(config(paths[name]) + fileName) + " " + str(value) + "\n"

    VMDin.stdin.write(sourceCNT)
    VMDin.stdin.write(CNTtools)
    VMDin.stdin.write(tclName)
    VMDin.stdin.flush()
    VMDin.stdin.close
    VMDin.communicate()

def confWrite(paths, fileName, restraint, temp, duration, minDuration):
    """ confWrite generates a .conf file to use as input to namd2. To organize simulations
    with different parameters, confWrite will create a directory for the simulation."""

    confFile = paths['tfinal'] + fileName + '.conf'

    if os.path.isfile(confFile):
    
        print("################################################################\n" \
              + "NAMD configuration file " + confFile + " already exists.\n" \
              + "################################################################\n")
        
        return confFile
          
    print("################################################################\n" \
          + "Writing NAMD configuration file " + confFile + ".\n" \
          + "################################################################\n")
      
    if not os.path.exists(paths['tfinal']):
        os.makedirs(paths['tfinal'])

    # Grabs the CNT basis vectors
    x, y, z = getCNTBasis(config(paths['pbc']) + fileName + '.pdb')

    outFile = open(confFile, "w")
    
    outFile.write("# Set up periodic boundary conditions\n\n")
    
    outFile.write("cellBasisVector1    {0:>10.3f}{1:>10.3f}{2:>10.3f}\n".format(x, 0., 0.))
    outFile.write("cellBasisVector2    {0:>10.3f}{1:>10.3f}{2:>10.3f}\n".format(0., y, 0.))
    outFile.write("cellBasisVector3    {0:>10.3f}{1:>10.3f}{2:>10.3f}\n".format(0., 0., z))
    outFile.write("cellOrigin          {0:>10}{1:>10}{2:>10}\n\n".format(0, 0, 0))
    
    outFile.write("{0:<20}".format("wrapWater") + "on" + "\n")
    outFile.write("{0:<20}".format("wrapAll") + "on" + "\n")
    outFile.write("\n")

    outFile.write("# Set the input and output files\n\n")
    
    outFile.write("{0:<20}".format("structure") + fileName + '.psf' + "\n")
    outFile.write("{0:<20}".format("coordinates") + fileName + '.pdb' + "\n")
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
    outFile.write("\n")

    outFile.write("# Set the integration parameters\n\n")
    outFile.write("{0:<20}".format("timestep") + "1" + "\n")
    outFile.write("{0:<20}".format("nonbondedFreq") + "2" + "\n")
    outFile.write("{0:<20}".format("fullElectFrequency") + "4" + "\n")
    outFile.write("{0:<20}".format("stepspercycle") + "20" + "\n")
    outFile.write("{0:<20}".format("rigidBonds") + "water" + "\n")

    outFile.write("\n")

    outFile.write("{0:<20}".format("restartfreq") + "100" + "\n")
    outFile.write("{0:<20}".format("dcdfreq") + "100" + "\n")
    outFile.write("{0:<20}".format("veldcdfreq") + "100" + "\n")
    outFile.write("{0:<20}".format("outputEnergies") + "100" + "\n")
    outFile.write("\n")

    outFile.write("# Set the restraints on the carbon\n\n")
    
    if restraint == 0:
        outFile.write("{0:<20}".format("fixedAtoms") + "on" + "\n")
        outFile.write("{0:<20}".format("fixedAtomsFile") + fileName + '-restraint.pdb' + "\n")
        outFile.write("{0:<20}".format("fixedAtomsCol") + "B" + "\n")
    else:
        outFile.write("{0:<20}".format("constraints") + "on" + "\n")
        outFile.write("{0:<20}".format("consref") + fileName + '-restraint.pdb' + "\n")
        outFile.write("{0:<20}".format("conskfile") + fileName + '-restraint.pdb' + "\n")
        outFile.write("{0:<20}".format("conskcol") + "O" + "\n")

    outFile.write("\n")
    
    outFile.write("# Set the external forcing\n\n")
    
    outFile.write("{0:<20}".format("constantForce") + "yes" + "\n")
    outFile.write("{0:<20}".format("consForceFile") + fileName + '-forcing.pdb' + "\n")
    outFile.write("{0:<20}".format("consForceScaling") + "0.014392631" + "\n")
    outFile.write("\n")

    outFile.write("# Set the Langevin thermostat\n\n")
    
    outFile.write("{0:<20}".format("langevin") + "on" + "\n")
    outFile.write("{0:<20}".format("langevinFile") + fileName + '-langevin.pdb' + "\n")
    outFile.write("{0:<20}".format("langevinCol") + "O" + "\n")
    outFile.write("{0:<20}".format("langevinTemp") + str(temp) + "\n")
    outFile.write("\n")

    outFile.write("# Set the execution parameters\n\n")

    outFile.write("{0:<20}".format("temperature") + str(temp) + "\n")
    outFile.write("{0:<20}".format("minimize") + str(minDuration) + "\n")
    outFile.write("{0:<20}".format("reinitvels") + str(temp) + "\n")
    outFile.write("{0:<20}".format("run") + str(duration) + "\n")
    
    outFile.close()
    paramFile = paths['home'] + "Templates/par_all27_prot_lipid.prm"
    shutil.copy(paramFile, paths['tfinal'])

    os.link(configFile(paths['solvate'], fileName + '.pdb'), paths['tfinal'] + fileName + '.pdb')
    os.link(configFile(paths['solvate'], fileName + '.psf'), paths['tfinal'] + fileName + '.psf')
    os.link(configFile(paths['restraint'], fileName + '.pdb'), paths['tfinal'] + fileName + '-restraint.pdb')
    os.link(configFile(paths['forcing'], fileName + '.pdb'), paths['tfinal'] + fileName + '-forcing.pdb')
    os.link(configFile(paths['temperature'], fileName + '.pdb'), paths['tfinal'] + fileName + '-langevin.pdb')
    
    return confFile
    
def tubeGen(inFile, pbcFile, N_0, n, m):
    """ tubeGen creates a periodic nanotube with the following input parameters:
    inFile is the name of the initial nanotube with number of rings N_0 and
    dimensions n x m.  pbcFile is the name of the same nanotube but now with 
    periodic boundary conditions applied to it. """

    # Bonds lengths of different armchair nanotubes in nanometers
    s0 = 0.1418
    # calculates the length of the nanotube based on bond lengths
    l = float((N_0-0.75))*s0*np.sqrt(3)
    
    paths = {}
    paths['home'] = os.path.abspath('..')
    paths['type'] = paths['home'] + '/Data/' + '(' + str(n) + ', ' + str(m) + ')/'
    paths['length'] = paths['type'] + 'N0 = ' + str(N_0) + '/'
    paths['pbc'] = paths['length'] + 'PBC/'
    
    files = {}
    # files['length-psf'] = config(paths['length']) + inFile + '.psf'
    # files['length-pdb'] = config(paths['length']) + inFile + '.pdb'
    # files['pbc-psf'] = config(paths['pbc']) + inFile + '.psf'
    # files['pbc-pdb'] = config(paths['pbc']) + inFile + '.pdb'
    # files['prebond-psf'] = config(paths['pbc']) + inFile + '-prebond.psf'
    # files['prebond-pdb'] = config(paths['pbc']) + inFile + '-prebond.pdb'
    files['length-psf'] = inFile + '.psf'
    files['length-pdb'] = inFile + '.pdb'
    files['pbc-psf'] = inFile + '.psf'
    files['pbc-pdb'] = inFile + '.pdb'
    files['prebond-psf'] = inFile + '-prebond.psf'
    files['prebond-pdb'] = inFile + '-prebond.pdb'
    
    # if os.path.isfile(files['pbc-psf']) \
    #     and os.path.isfile(files['pbc-pdb']) \
    #     and os.path.isfile(files['prebond-psf']) \
    #     and os.path.isfile(files['prebond-pdb']):
    if os.path.isfile(configFile(paths['pbc'], files['pbc-psf'])) \
        and os.path.isfile(configFile(paths['pbc'], files['pbc-pdb'])) \
        and os.path.isfile(configFile(paths['pbc'], files['prebond-psf'])) \
        and os.path.isfile(configFile(paths['pbc'], files['prebond-pdb'])):
        
        print("################################################################\n" \
              + "Configuration files already exist in:\n" + config(paths['pbc']) + ".\n" \
              + "################################################################\n")

        return paths, files

    else:
    
        print("################################################################\n" \
              + "Running VMD to generate nanotube configuration files.\n" \
              + "################################################################\n")

        if not os.path.exists(config(paths['pbc'])):
            os.makedirs(config(paths['pbc']))

        if not os.path.exists(config(paths['length'])):
            os.makedirs(config(paths['length']))

        logFile = open(config(paths['length']) + inFile + ".log", "w")

        # Opening a pipe to VMD in the shell
        VMDin=subprocess.Popen(["vmd", "-dispdev", "none"], stdin=subprocess.PIPE, \
                           stdout=logFile)

        # runs CNTtools.tcl script to generate nanotube and generate PBCs
        sourceCNT = "source CNTtools.tcl\n"
        CNTtools = "package require CNTtools 1.0\n"

        genNT = "genNT " + quoted(inFile) + " " + quoted(config(paths['length'])) + " " \
            + str(l) + " " + str(n) + " " + str(m) + "\n"

        pbcNT = "pbcNT " + quoted(config(paths['length']) + inFile) + " " \
            + quoted(config(paths['pbc']) + pbcFile) + " default\n"

        fixNT = "fixNT " + quoted(config(paths['pbc']) + pbcFile) + "\n"

        removeLangevinWater = "removeLangevinWater " + quoted(config(paths['pbc']) + pbcFile) + "\n"

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
            return paths, files


def solvate(inFile, N_0, S, n, m):
    """ The solvate module will create an armchair swcnt with N_0 rings and
    add N_0+S water molecules inside the nanotube, then write out new psf and
    pdb files for the nanotube. """
    paths, files = tubeGen(inFile, inFile, N_0, n, m)

    paths['solvate'] = paths['pbc'] + 'S = ' + str(S) + '/'
    
    # files['solvate-pdb'] = config(paths['solvate']) + inFile + '-solv.pdb'
    # files['solvate-psf'] = config(paths['solvate']) + inFile + '-solv.psf'
    files['solvate-pdb'] = inFile + '.pdb'
    files['solvate-psf'] = inFile + '.psf'
   
    if os.path.isfile(configFile(paths['solvate'], files['solvate-psf'])) \
        and os.path.isfile(configFile(paths['solvate'], files['solvate-pdb'])):
    
        print("################################################################\n" \
              + "Configuration files already exist in:\n" + config(paths['solvate']) + ".\n" \
              + "################################################################\n")
        
        return paths, files

    else:
    
        print("################################################################\n" \
          + "Adding water to the CNT and writing configuration files.\n" \
          + "################################################################\n")
    
        if not os.path.exists(config(paths['solvate'])):
            os.makedirs(config(paths['solvate']))
    
        # Opens input nanotube psf and pdb files, and reads all the lines of each file into lists
        with open(configFile(paths['pbc'], files['pbc-psf'])) as psfFile:
            psfLines = psfFile.readlines()
        with open(configFile(paths['pbc'], files['pbc-pdb'])) as pdbFile:
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

        # Iterates through all of the lines of the input Psf file, and records the
        # index of the lines that contain the string !NATOM, !NBOND, and !NTHETA,
        # as wella as changes the line to update the new number of each
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
        for i in range(nAtoms+1, (3*(N_0+S)) + nAtoms+1, 3):
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

        sFactorAdjust = float(N_0) / (N_0 + S)

        for i in range(0,3*(N_0+S), 3):
            if i==0:
                pdbLines[lenPdb-1] = oxygen.format(nAtoms+1, i/3+2, 0.000)
                pdbLines.append( hydro1.format(nAtoms+(i+2), i/3+2, 0.570) )
                pdbLines.append( hydro2.format(nAtoms+(i+3), i/3+2, 0.570) )
            
            else:
                pdbLines.append( oxygen.format(nAtoms+i+1, i/3+2, ((i/3)*dist*sFactorAdjust)) )
                pdbLines.append( hydro1.format(nAtoms+(i+2), i/3+2, 0.570+((i/3)*dist*sFactorAdjust)) )
                pdbLines.append( hydro2.format(nAtoms+(i+3), i/3+2, 0.570+((i/3)*dist*sFactorAdjust)) )

        # Writes the new pdb lines to a new pdb file
        pdbLines.append("END\n")
        pdbOut = open(configFile(paths['solvate'], files['solvate-pdb']), 'w')
        pdbOut.writelines(pdbLines)
        pdbOut.close()

        # Writes the new psf lines to a new psf file
        psfOut = open(configFile(paths['solvate'], files['solvate-psf']), 'w')
        psfOut.writelines(preAtoms)
        psfOut.writelines(psfLines[atomIndex])
        psfOut.writelines(atoms)
        psfOut.writelines("\n" + psfLines[bondIndex])
        psfOut.writelines(bondsFinal)
        psfOut.writelines("\n" + psfLines[angleIndex])
        psfOut.writelines(anglesFinal)
        psfOut.writelines(postAngles)
        psfOut.close()

        return paths, files


def forceWrite(inFile, paths, files, force):
    # Generates the file that indicates how strong the force is on each atom by modifying the occupancy column

    paths['forcing'] = paths['restraint'] + 'F = ' + str(force) + '/'
    
    # files['forcing-pdb'] = config(paths['forcing']) + inFile + '-forcing.pdb'
    files['forcing-pdb'] = inFile + '-forcing.pdb'
    
    if os.path.isfile(configFile(paths['forcing'], files['forcing-pdb'])):
        
        print("################################################################\n" \
              + "Configuration files already exist in:\n" + config(paths['forcing']) + ".\n" \
              + "################################################################\n")
              
        return paths, files

    else:
        
        print("################################################################\n" \
              + "Writing configuration file for external forcing.\n" \
              + "################################################################\n")
    
        if not os.path.exists(config(paths['forcing'])):
            os.makedirs(config(paths['forcing']))
        
        forceVal = "{0:.2f}".format(force)
        forceVal = forceVal[0:4]

        with open (configFile(paths['solvate'], files['solvate-pdb'])) as f:
            flines = f.readlines()

        for i in range(1, len(flines)-1):
            if "CNT" in flines[i]:
                flines[i] = flines[i][0:56] + "0.00" + flines[i][60::]
            elif "OH2" in flines[i]:
                    flines[i] = flines[i][0:33] + "0.000   0.000   1.000  " + forceVal + flines[i][60::]
            else:
                flines[i] = flines[i][0:56] + "0.00" + flines[i][60::]

        outFile = open(configFile(paths['forcing'], files['forcing-pdb']), 'w')
        outFile.writelines(flines)
        outFile.close()

        return paths, files

        
def restraintWrite(inFile, paths, files):
    # Generate a file that restrains each carbon atom to it's initial position with force constant k

    kVal = 3.0
    
    if kVal == 0.0:
        paths['restraint'] = paths['solvate'] + 'Fixed/'
    
    else:
        paths['restraint'] = paths['solvate'] + 'R = ' + str(kVal) + '/'
    
    # files['restraint-pdb'] = config(paths['restraint']) + inFile + '-restraint.pdb'
    files['restraint-pdb'] = inFile + '-restraint.pdb'
    
    if os.path.isfile(configFile(paths['restraint'], files['restraint-pdb'])):
        
        print("################################################################\n" \
              + "Configuration files already exist in:\n" + config(paths['restraint']) + ".\n" \
              + "################################################################\n")
              
        return paths, files

    else:
    
        if kVal == 0:
            print("################################################################\n" \
                  + "Writing configuration file for fixed nanotube.\n" \
                  + "################################################################\n")
        
        else:
            print("################################################################\n" \
                  + "Writing configuration file for flexible nanotube with R = " + str(kVal) + ".\n" \
                  + "################################################################\n")
        
        
        if not os.path.exists(config(paths['restraint'])):
            os.makedirs(config(paths['restraint']))

        with open (configFile(paths['solvate'], files['solvate-pdb'])) as kFile:
            kfLines = kFile.readlines()
        
        if kVal == 0.0:
            for i in range(1, len(kfLines)-1):
                if "CNT" in kfLines[i]:
                    kfLines[i] = kfLines[i][0:62] + "1.00" + kfLines[i][66::]
                        
        else:
            for i in range(1, len(kfLines)-1):
                if "CNT" in kfLines[i]:
                    kfLines[i] = kfLines[i][0:56] + "{:.2f}".format(kVal) + kfLines[i][60::]
                else:
                    kfLines[i] = kfLines[i][0:56] + "{:.2f}".format(0.00) + kfLines[i][60::]

        outkFile = open(configFile(paths['restraint'], files['restraint-pdb']), 'w')
        outkFile.writelines(kfLines)
        outkFile.close()

        return paths, files


def simWrite(inFile, paths, files, temp = 300, tf = 20000, minimize = 1000):
    """ simWrite generates a .conf file to use as input to namd2. To organize simulations
    with different parameters, simWrite will create a directory for the simulation."""

    paths['temperature'] = paths['forcing'] + 'Temp = ' + str(temp) + '/'
    paths['tfinal'] = paths['temperature'] + 'Tf = ' + str(tf) + '/'
    
    files['sim-conf'] = paths['tfinal'] + inFile + '.conf'

    if os.path.isfile(files['sim-conf']):
    
        print("################################################################\n" \
              + "NAMD configuration file already exists in:\n" + config(paths['tfinal']) + ".\n" \
              + "################################################################\n")
        
        return paths, files
          
    else:
        
        print("################################################################\n" \
              + "Writing NAMD configuration file for T = " + str(temp) + ", Tf = " + str(tf) + ".\n" \
              + "################################################################\n")
          
        if not os.path.exists(paths['tfinal']):
            os.makedirs(paths['tfinal'])

        # Grabs the CNT basis vectors
        x, y, z = getCNTBasis(config(paths['pbc']) + files['prebond-pdb'])

        # Read in lines of simulation file
        with open(paths['home'] + "/Templates/sim_template.conf") as tempFile:
            simLines = tempFile.readlines()

        simLines[12] = "structure          " + quoted(files['solvate-psf']) + "\n"
        simLines[13] = "coordinates        " + quoted(files['solvate-pdb']) + "\n"

        simLines[15] = "set temperature    {:3d}\n".format(temp)
        simLines[30] = "cellBasisVector1    {0:<10.3f}{1:<10}{2:}\n".format(x, 0., 0.)
        simLines[31] = "cellBasisVector2    {0:<10}{1:<10.3f}{2:}\n".format(0., y, 0.)
        simLines[32] = "cellBasisVector3    {0:<10}{1:<10}{2:.3f}\n".format(0., 0., z)
        # simLines[33] = "cellOrigin          {0:<10}{1:<10}{2:.3f}\n\n".format(0, 0, float(z)/2 )
        simLines[33] = "cellOrigin          {0:<10}{1:<10}{2:.3f}\n\n".format(0, 0, 0)

        # simLines[16] = "set outputname     " + quoted(paths['tfinal'] + inFile) + "\n"
        simLines[16] = "set outputname     " + quoted(inFile) + "\n"
        simLines[71] = "fixedAtomsFile      " + files['solvate-pdb'] + "\n"
        simLines[75] = "consref             " + files['restraint-pdb'] + "\n"
        simLines[76] = "conskfile            " + files['restraint-pdb'] + "\n"
        simLines[78] = "constraintScaling     " + "{:.2f}".format(1.00)

        simLines[94] = "consforcefile         " + files['forcing-pdb'] + "\n"
        simLines[100] = "minimize " + str(minimize) + "\n"
        simLines[102] = "run {:5d} \n".format(tf)

        # Write contents out to original file
        
        outFile = open(files['sim-conf'], "w")
        outFile.writelines(simLines)
        outFile.close()
        paramFile = paths['home'] + "/Templates/par_all27_prot_lipid.prm"
        shutil.copy(paramFile, paths['tfinal'])

        os.link(configFile(paths['solvate'], files['solvate-pdb']), paths['tfinal'] + files['solvate-pdb'])
        os.link(configFile(paths['solvate'], files['solvate-psf']), paths['tfinal'] + files['solvate-psf'])
        os.link(configFile(paths['restraint'], files['restraint-pdb']), paths['tfinal'] + files['restraint-pdb'])
        os.link(configFile(paths['forcing'], files['forcing-pdb']), paths['tfinal'] + files['forcing-pdb'])
        return paths, files


def runSim(simPath, simFile, output):
    """ given an input path to a simulation file, runSim will call namd2 to run the simulation """
    
    logFile = open(simPath + output + ".log", "w")
    
    print("################################################################\n" \
          + "Starting simulation.\n" \
          + "################################################################\n")
    
    startTime = time.time()

    Namd2in=subprocess.Popen(["namd2", simFile], stdin=subprocess.PIPE, stdout=logFile, stderr=logFile)
    Namd2in.stdin.flush()
    Namd2in.stdin.close
    Namd2in.communicate()
    
    elapsedTime = time.time() - startTime
    
    if Namd2in.returncode==0:
        print("################################################################\n" \
              + "Simulation finished\nSimulation file saved into " \
              + simPath.replace(".conf", ".dcd") + ".\n" \
              + "Elapsed time: " + "{:.0f}".format(elapsedTime) + " s = " \
              + "{:.1f}".format(elapsedTime/60) + " m = " \
              + "{:.1f}".format(elapsedTime/3600) + " h\n" \
              + "################################################################\n")
    else:
        print("################################################################\n" \
              + "Simulation FAIL.\n" \
              + "################################################################\n")


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


def main(FNAME, N_0, S, n, m, TEMP, LENGTH, FORCESTRENGTH, minimize):
    paths, files = solvate(FNAME, N_0, S, n, m)
    paths, files = restraintWrite(FNAME, paths, files)
    paths, files = forceWrite(FNAME, paths, files, FORCESTRENGTH)
    paths, files = simWrite(FNAME, paths, files, TEMP, LENGTH, minimize)
    runSim(paths['tfinal'], files['sim-conf'], FNAME)
    
def run(fileName = 'Test', N0 = 20, S = 0, n = 5, m = 5, temp = 300, damping = 1, \
        duration = 1000, force = 0, restraint = 0, minDuration = 1000):
    
    # The paths dictionary holds the path name to each tier of the data hierarchy
    paths = {}
    paths['home'] = os.path.abspath('..') + '/'
    paths['type'] = paths['home'] + 'Data/' + '(' + str(n) + ', ' + str(m) + ')/'
    paths['N0'] = paths['type'] + 'N0 = ' + str(N0) + '/'
    paths['pbc'] = paths['N0'] + 'PBC/'
    paths['solvate'] = paths['pbc'] + 'S = ' + str(S) + '/'
    paths['restraint'] = paths['solvate'] + 'R = ' + str(restraint) + '/'
    paths['forcing'] = paths['restraint'] + 'F = ' + str(force) + '/'
    paths['temperature'] = paths['forcing'] + 'Temp = ' + str(temp) + ', d = ' + str(damping) + '/'
    paths['tfinal'] = paths['temperature'] + 'Tf = ' + str(duration) + '/'
        
    cntWrite(paths, fileName, N0, n, m)
    waterWrite(paths, fileName, N0, S)
    pdbWrite(paths, fileName, 'restraint', restraint)
    pdbWrite(paths, fileName, 'forcing', force)
    pdbWrite(paths, fileName, 'temperature', damping)
    confFile = confWrite(paths, fileName, restraint, temp, duration, minDuration)
    runSim(paths['tfinal'], confFile, fileName)