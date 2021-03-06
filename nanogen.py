### Written by Matthew Kwiecien Jul, 2015
### Modified by CGG, spring 2016
### Additional modifications, and updating to python 3, summer 2017

# from __future__ import print_function
# from __future__ import division

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

def run(pD):
    """ This is the main entry point for performing simulations with a carbon nanotube
    and water. It takes a python dictionary of simulation parameters as its
    argument. Before calling run, you must call myutil.setWaterParamDefaults()
    to create the dictionary."""
    
    # The paths dictionary holds the path name to each tier of the data hierarchy.
    
    paths = my.makeWaterPaths(pD)
    
    # Find the hostname and set the name of the wisdom file used by PME.
    
    hostname = socket.gethostname().partition('.')[0]
    wisdomFile = "FFTW_NAMD_2.11_" + hostname + ".txt"

    # If we're doing a production run (not a minimization) on nanotube, we want to
    # use the GPU, so we run the CUDA version of namd2. The CUDA version of namd seems
    # to give errors occasionally when minimizing, so we don't use it for that.
    
    if pD['Run Type'] == 'Minimize':
        executable = 'namd2'
    else:
        if hostname == 'nanotube':
            executable = 'namd2-cuda'
        else:
            executable = 'namd2'
        
    # If the data folder for the simulation already exists, abort, under the assumption
    # that this is an accidental duplicate of an earlier run. New runs should go into
    # their own folders.
    
    if os.path.exists(paths['data']):
        print("################################################################\n" \
              + "Data folder for " + pD['Run Type'] + " run already exists ... aborting.\n" \
              + "################################################################\n")
        return

    # This is the base name for all the various pdb and psf files we will generate.
    
    fileName = pD['Config File Name']

    # Bonds lengths of armchair nanotubes in nanometers. From this we'll calculate the
    # actual length of the nanotube in nanometers to feed to VMD.
    
    s0 = 0.1418
    l = float(pD['N0']-0.75) * s0 * np.sqrt(3)

    # Create the desired nanotube. The files will be stored in the 'Config Files' folder
    # in the folder pointed to by paths['N0'].
    
    package = "package require CNTtools 1.0\n"
    commands = [ "genNT " \
                + my.quoted(my.configFile(paths['N0'], fileName)) + " " \
                + str(l) + " " + str(pD['n']) + " " + str(pD['m']) + "\n" ]
   
    my.pdbWrite(paths['N0'], fileName, package, commands)
    
    # Now we'll make the nanotube periodic in the z direction.
    
    commands = [ "pbcNT " + my.quoted(my.configFile(paths['N0'], fileName)) + " " \
                + my.quoted(my.configFile(paths['pbc'], fileName)) + " default\n" ]

    my.pdbWrite(paths['pbc'], fileName, package, commands)  
    
    # Delete the now-uneeded -prebond pdb and psf files.
    
    if os.path.exists(my.config(paths['pbc']) + fileName + "-prebond.pdb"):
        os.remove(my.config(paths['pbc']) + fileName + "-prebond.pdb")
    if os.path.exists(my.config(paths['pbc']) + fileName + "-prebond.psf"):
        os.remove(my.config(paths['pbc']) + fileName + "-prebond.psf")
    
    # Read back in the pdb file and make sure the z-coordinate spacing is regular.
    
    alignCNT(my.config(paths['pbc']) + fileName + ".pdb")

    # cntWrite(paths, pD['Config File Name'], pD['N0'], pD['n'], pD['m'])
    
    # Add the water to the carbon nanotube, and center the nanotube in the box.
    
    waterWrite(paths, pD['Config File Name'], pD['N0'], pD['S'])
    
    # Write the restraint pdb file.
    
    commands = [ "NTrestraint " \
                + my.quoted(my.configFile(paths['solvate'], fileName)) + " " \
                + my.quoted(my.configFile(paths['restraint'], fileName)) + " " \
                + str(pD['Restraint']) + "\n" ]
    my.pdbWrite(paths['restraint'], fileName, package, commands)
    
    # If the simulation is a New (not Minimization) run, we need to set the forcing
    # and the thermostat parameters, and preheat the water.
    
    if pD['Run Type'] == 'New':
    
        # Write the forcing pdb file.
    
        commands = [ "NTforcing " \
                    + my.quoted(my.configFile(paths['solvate'], fileName)) + " " \
                    + my.quoted(my.configFile(paths['forcing'], fileName)) + " " \
                    + str(pD['Force (pN)']) + "\n" ]
        my.pdbWrite(paths['forcing'], fileName, package, commands)

        # Write the Langevin pdb file.
    
        commands = [ "NTtemperature " \
                    + my.quoted(my.configFile(paths['solvate'], fileName)) + " " \
                    + my.quoted(my.configFile(paths['temperature'], fileName)) + " " \
                    + str(pD['Damping']) + "\n" ]
        my.pdbWrite(paths['temperature'], fileName, package, commands)

        # Set up the temporary folder for the simulation.
        
        paths['run'] = os.getenv('HOME') + '/' + pD['Data Folder'] + '/Preheat/'
        
        # Create a NAMD configuration file for the preheating run.
        
        preheatConfFile = confWrite(paths, pD, hostname, wisdomFile, preheat = True)
        
        # Link all the necessary files back to the temporary run folder.
        
        linkFiles(paths, pD, wisdomFile, preheat = True)
        
        # Do the preheating run.
        
        result, runTime = runSim(executable, paths['run'], preheatConfFile, \
            pD['Config File Name'], hostname)
    
    # Set up a temporary run folder. This is either for a production run or for a
    # minimization run.
    
    paths['run'] = os.getenv('HOME') + '/' + pD['Data Folder'] + '/'
    
    # Create the NAMD configuration file for the run. Again, this is a production run
    # or a minimization run.
    
    confFile = confWrite(paths, pD, hostname, wisdomFile)
    
    # Link all the necessary files back to the temporary run folder.

    linkFiles(paths, pD, wisdomFile)

    # Do the simulation!!
    
    result, runTime = runSim(executable, paths['run'], confFile, pD['Config File Name'], hostname)
    
    # If the run succeeded, log it to the appropriate spreadsheet.
    
    if result == 0:
    
        if not os.path.isfile(hostname + '.csv'):
        
            print("################################################################\n" \
                + "Creating file " + hostname + ".csv." + "\n" \
                + "################################################################\n")
        
            outFile = open(hostname + '.csv', "w")
            outFile.write("Date,Run Time (h),File Name,Folder,Run Type,N0,S,n,m," \
                + "Temperature (K),Damping,Thermostat," \
                + "Force (pN),PME,Restraint,Duration,Min Duration,dt (fs),outputFreq\n")
                
        print("################################################################\n" \
              + "Writing to file " + hostname + ".csv." + "\n" \
              + "################################################################\n")
              
        outFile = open(hostname + '.csv', "a")
        outFile.write(str(datetime.date.today()) + ',' \
            + "{:.5f}".format(runTime) + ',' \
            + pD['Data Folder'] + ',' \
            + pD['Top Folder'] + ',' \
            + pD['Run Type'] + ',' \
            + str(pD['N0']) + ',' \
            + str(pD['S']) + ',' \
            + str(pD['n']) + ',' \
            + str(pD['m']) + ',' \
            + str(pD['Temperature (K)']) + ',' \
            + str(pD['Damping']) + ',' \
            + pD['Thermostat'] + ',' \
            + str(pD['Force (pN)']) + ',' \
            + pD['PME'] + ','
            + str(pD['Restraint']) + ',' \
            + str(pD['Duration']) + ',' \
            + str(pD['Min Duration']) + ',' \
            + str(pD['dt (fs)']) + ',' \
            + str(pD['outputFreq']) + '\n')
    
    print("################################################################\n" \
          + "Finished writing to file " + hostname + ".csv." + "\n" \
          + "################################################################\n")
    print("################################################################\n" \
          + "Moving " + paths['run'] + " to " + paths['data'] + ".\n" \
          + "################################################################\n")
    
    # Move the temporary run file to the appropriate place in the data hierarchy.
    
    shutil.move(paths['run'], paths['data'])
    
def cntWrite(paths, fileName, N0, n, m):
    """ cntWrite creates a periodic nanotube with the following input parameters:
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

    CNTtools = "package require CNTtools 1.0\n"

    genNT = "genNT " + my.quoted(fileName) + " " + my.quoted(my.config(paths['N0'])) + " " \
        + str(l) + " " + str(n) + " " + str(m) + "\n"

    # run commands through pipe and saves to file
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

    CNTtools = "package require CNTtools\n"
    pbcNT = "pbcNT " + my.quoted(my.config(paths['N0']) + fileName) + " " \
        + my.quoted(my.config(paths['pbc']) + fileName) + " default\n"

    # run commands through pipe and saves to file
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

    if not my.needConfigFiles(paths['solvate'], fileName):
        return
    
    # Opens input nanotube psf and pdb files, and reads all the lines of each file into lists

    with open(my.configFile(paths['pbc'], fileName + '.psf')) as psfFile:
        psfLines = psfFile.readlines()
    with open(my.configFile(paths['pbc'], fileName + '.pdb')) as pdbFile:
        pdbLines = pdbFile.readlines()

    # Grabs the lengths of each of the lists

    lenPsf = len(psfLines)
    lenPdb = len(pdbLines)

    # String formats for the PSF and PDB file writing
    # PDB

    dampingCoeff = 0.00
    oxygen = "ATOM{0:>7}  OH2 TIP3 {1:4.0f}       0.000   0.000{2:>8.3f}  0.00  0.00      WTR  O\n"
    hydro1 = "ATOM{0:>7}  H1  TIP3 {1:4.0f}       0.000   0.766{2:>8.3f}  0.00  0.00      WTR  H\n"
    hydro2 = "ATOM{0:>7}  H2  TIP3 {1:4.0f}       0.000  -0.766{2:>8.3f}  0.00  0.00      WTR  H\n"

    # Psf
    
    opsf  = "   {0:>5} WTR  {1:<4} TIP3 OH2  OT    -0.834000       15.9994           0\n"
    h1psf = "   {0:>5} WTR  {1:<4} TIP3 H1   HT     0.417000        1.0080           0\n"
    h2psf = "   {0:>5} WTR  {1:<4} TIP3 H2   HT     0.417000        1.0080           0\n"

    # String format for the bonds and angles in the psf file

    sBondFormat  = " {0: >8}{1: >8}{2: >8}{3: >8}{4: >8}{5: >8}{6: >8}{7: >8}\n"
    sAngleFormat = " {0: >8}{1: >8}{2: >8}{3: >8}{4: >8}{5: >8}{6: >8}{7: >8}{8: >8}\n"

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

    oxZ = waterZ(N0, S)

    for i in range(0, N0+S):

        hyZ = oxZ[i] + 0.570
        if i == 0:
            pdbLines[lenPdb-1] = oxygen.format(nAtoms+1, 2, oxZ[i])
            pdbLines.append( hydro1.format(nAtoms+2, 2, hyZ) )
            pdbLines.append( hydro2.format(nAtoms+3, 2, hyZ) )
        else:
            pdbLines.append( oxygen.format(nAtoms + 3*i+1, i+2, oxZ[i]) )
            pdbLines.append( hydro1.format(nAtoms + (3*i+2), i+2, hyZ) )
            pdbLines.append( hydro2.format(nAtoms + (3*i+3), i+2, hyZ) )

    # Writes the new pdb lines to a new pdb file

    pdbLines.append("END\n")
    pdbOut = open(my.configFile(paths['solvate'], fileName + '-temp.pdb'), 'w')
    pdbOut.writelines(pdbLines)
    pdbOut.close()

    # Writes the new psf lines to a new psf file

    psfOut = open(my.configFile(paths['solvate'], fileName + '-temp.psf'), 'w')
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

    package = "package require CNTtools 1.0\n"
    commands = [ "centerNT " + my.quoted(my.configFile(paths['solvate'], fileName)) + "\n" ]
    my.pdbWrite(my.config(paths['solvate']), fileName, package, commands) 

    # Delete the now-uneeded -temp pdb and psf files.
    
    os.remove(my.config(paths['solvate']) + fileName + "-temp.pdb")
    os.remove(my.config(paths['solvate']) + fileName + "-temp.psf")

def waterZ(N0, S, equalSpacing = False):

    # Bonds lengths for different armchair nanotubes, in angstroms
    s0 = 1.418
    
    # Calculates the apothem of each regular hexagon in the nanotube and then
    # calculates the distance between the center of each ring in the tube
    # apothem = s*np.sqrt(3)/2
    dist = s0*np.sqrt(3)

    if equalSpacing:
    
        if N0 + S != 0:
            sFactorAdjust = float(N0) / (N0 + S)
        else:
            sFactorAdjust = 0
    
        oxZ = dist/4 + np.array(range(0, N0+S))*dist*sFactorAdjust
    
    else:
    
        x = np.array(range(0, N0+S))/float(N0+S)
        width = -np.sign(S)*N0/10

        x0 = np.array(range(0, 2*np.abs(S)))/(2.0*np.abs(S)) + 1.0/(4*np.abs(S))

        z = np.zeros(N0+S)

        for i in range(0,len(x0)):
            z = z + np.arctan(np.exp(width*(x-x0[i])))/np.pi
            
        if S > 0:
            z = z - 1

        oxZ = dist/4 + np.array(range(0, N0+S))*dist + dist*z
        
    return oxZ
    
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

    CNTtools = "package require CNTtools 1.0\n"

    tclName = "NT" + name + " " + my.quoted(my.config(paths['solvate']) + fileName) + " " \
        + my.quoted(my.config(paths[name]) + fileName) + " " + str(value) + "\n"

    VMDin.stdin.write(CNTtools)
    VMDin.stdin.write(tclName)
    VMDin.stdin.flush()
    VMDin.stdin.close
    VMDin.communicate()

def confWrite(paths, pD, hostname, wisdomFile, preheat = False):
    """ confWrite generates a .conf file to use as input to namd2. To organize simulations
    with different parameters, confWrite will create a directory for the simulation."""

    confFile = paths['run'] + pD['Config File Name'] + '.conf'
    if not os.path.exists(paths['run']):
        os.makedirs(paths['run'])

    print("################################################################\n" \
          + "Writing NAMD configuration file " + confFile + ".\n" \
          + "################################################################\n")
      
    # Grabs the CNT basis vectors
    x, y, z = getCNTBasis(my.config(paths['pbc']) + pD['Config File Name'] + ".pdb")

    outFile = open(confFile, "w")
    
    outFile.write("# Configuration file written on host " + hostname + "\n" \
        + "# on " + str(datetime.datetime.now().date()) \
        + " at " + str(datetime.datetime.now().time()).partition('.')[0] + ".\n\n")
        
    outFile.write("# Limit the length of the log file\n\n")
    
    outFile.write("{0:<20}".format("outputEnergies") + str(100*pD['outputFreq']) + "\n\n")
    
    outFile.write("# Set up periodic boundary conditions\n\n")
    
    outFile.write("cellBasisVector1    {0:>10.3f}{1:>10.3f}{2:>10.3f}\n".format(x, 0., 0.))
    outFile.write("cellBasisVector2    {0:>10.3f}{1:>10.3f}{2:>10.3f}\n".format(0., y, 0.))
    outFile.write("cellBasisVector3    {0:>10.3f}{1:>10.3f}{2:>10.3f}\n".format(0., 0., z))
    outFile.write("cellOrigin          {0:>10}{1:>10}{2:>10}\n\n".format(0, 0, 0))
    
    outFile.write("{0:<20}".format("wrapWater") + "off" + "\n")
    outFile.write("{0:<20}".format("wrapAll") + "off" + "\n")
    outFile.write("\n")

    outFile.write("# Set the input files\n\n")
    
    outFile.write("{0:<20}".format("structure") + pD['Config File Name'] + ".psf" + "\n")
    outFile.write("{0:<20}".format("coordinates") + pD['Config File Name'] + ".pdb" + "\n")
    outFile.write("\n")
    
    outFile.write("# Set the force field parameters\n\n")
    
    outFile.write("{0:<20}".format("paraTypeCharmm") + "on" + "\n")
    outFile.write("{0:<20}".format("parameters") + "par_all27_prot_lipid.prm" + "\n") 
    outFile.write("{0:<20}".format("exclude") + "scaled1-4" + "\n")
    outFile.write("{0:<20}".format("cutoff") + "12.0" + "\n")
    outFile.write("{0:<20}".format("pairlistdist") + "14.0" + "\n")
    outFile.write("{0:<20}".format("switching") + "on" + "\n")
    outFile.write("{0:<20}".format("switchdist") + "10.0" + "\n")
    if pD['PME'] == 'on':
        outFile.write("{0:<20}".format("PME") + "yes" + "\n")
        outFile.write("{0:<20}".format("PMEGridSpacing") + "1.0" + "\n")
        outFile.write("{0:<20}".format("FFTWWisdomFile") + wisdomFile + "\n")
    else:
        outFile.write("{0:<20}".format("PME") + "no" + "\n")
    outFile.write("\n")

    outFile.write("# Set the integration parameters\n\n")
    
    outFile.write("{0:<20}".format("timestep") + str(pD['dt (fs)']) + "\n")
    outFile.write("{0:<20}".format("nonbondedFreq") + "2" + "\n")
    outFile.write("{0:<20}".format("fullElectFrequency") + "4" + "\n")
    outFile.write("{0:<20}".format("stepspercycle") + "20" + "\n")
    outFile.write("{0:<20}".format("rigidBonds") + "water" + "\n")

    outFile.write("\n")

    outFile.write("{0:<20}".format("restartfreq") + str(pD['outputFreq']) + "\n")
    outFile.write("{0:<20}".format("dcdfreq") + str(pD['outputFreq']) + "\n")
    outFile.write("{0:<20}".format("veldcdfreq") + str(pD['outputFreq']) + "\n")
    outFile.write("{0:<20}".format("forcedcdfreq") + str(pD['outputFreq']) + "\n")
    
    outFile.write("\n")

    outFile.write("# Set the restraints on the carbon\n\n")
    
    if pD['Restraint'] == 0:
        outFile.write("{0:<20}".format("fixedAtoms") + "on" + "\n")
        outFile.write("{0:<20}".format("fixedAtomsFile") + pD['Config File Name'] + "-restraint.pdb" + "\n")
        outFile.write("{0:<20}".format("fixedAtomsCol") + "B" + "\n")
    else:
        outFile.write("{0:<20}".format("constraints") + "on" + "\n")
        outFile.write("{0:<20}".format("consref") + pD['Config File Name'] + "-restraint.pdb" + "\n")
        outFile.write("{0:<20}".format("conskfile") + pD['Config File Name'] + "-restraint.pdb" + "\n")
        outFile.write("{0:<20}".format("conskcol") + "O" + "\n")

    outFile.write("\n")
    
    if pD['Run Type'] == 'Minimize':
    
        outFile.write("{0:<20}".format("outputname") + pD['Data File Name'] + "\n")
        outFile.write("{0:<20}".format("temperature") + '0' + "\n")
        outFile.write("{0:<20}".format("minimize") + str(pD['Min Duration']) + "\n")
        
        outFile.close()
        return confFile
        
    elif preheat:
    
        pT = {}
        
        pT['5', '0'] = 5
        pT['5', '0.005'] = 8
        pT['5', '0.01'] = 20
        pT['5', '0.02'] = 80
        pT['5', '0.04'] = 150
        
        pT['20', '0'] = 20
        pT['20', '0.005'] = 25
        pT['20', '0.01'] = 40
        pT['20', '0.02'] = 80
        pT['20', '0.04'] = 160
        
        pT['80', '0'] = 80
        pT['80', '0.005'] = 90
        pT['80', '0.01'] = 95
        pT['80', '0.02'] = 120
        pT['80', '0.04'] = 180
        
        pT['300', '0'] = 300
        pT['300', '0.005'] = 320
        pT['300', '0.01'] = 340
        pT['300', '0.02'] = 340
        pT['300', '0.04'] = 360
        
        preheatTemp = str(pT[str(pD['Temperature (K)']), str(pD['Force (pN)'])])
        outFile.write("{0:<20}".format("langevin") + "on" + "\n")
        outFile.write("{0:<20}".format("langevinFile") + pD['Config File Name'] + "-langevin.pdb" + "\n")
        # outFile.write("{0:<20}".format("langevinTemp") + str(pD['Temperature (K)']) + "\n")
        outFile.write("{0:<20}".format("langevinTemp") + preheatTemp + "\n")
        if pD['Restraint'] == 0:
            outFile.write("{0:<20}".format("langevinCol") + "X" + "\n")
        else:
            outFile.write("{0:<20}".format("langevinCol") + "B" + "\n")
        
        outFile.write("\n")

        outFile.write("{0:<20}".format("outputname") + pD['Data File Name'] + "\n")
        outFile.write("{0:<20}".format("temperature") + preheatTemp + "\n")
        outFile.write("{0:<20}".format("bincoordinates") + pD['Config File Name'] + ".restart.coor" + "\n")
        outFile.write("{0:<20}".format("reinitvels") + preheatTemp + "\n")
        outFile.write("{0:<20}".format("run") + str(10000) + "\n")        

        outFile.close()
        return confFile

    outFile.write("# Set the external forcing\n\n")
    
    if pD['Force (pN)'] == 0:
        outFile.write("{0:<20}".format("constantForce") + "no" + "\n")
    else:
        outFile.write("{0:<20}".format("constantForce") + "yes" + "\n")
        outFile.write("{0:<20}".format("consForceScaling") + "0.014392621" + "\n")
        outFile.write("{0:<20}".format("consForceFile") + pD['Config File Name'] + "-forcing.pdb" + "\n")
    
    outFile.write("\n")

    outFile.write("# Set the Langevin thermostat\n\n")

    if pD['Restraint'] == 0 or pD['Thermostat'] == 'Off':
        outFile.write("{0:<20}".format("langevin") + "off" + "\n")
    else:
        outFile.write("{0:<20}".format("langevin") + "on" + "\n")
        outFile.write("{0:<20}".format("langevinFile") + pD['Config File Name'] + "-langevin.pdb" + "\n")
        outFile.write("{0:<20}".format("langevinCol") + "O" + "\n")
        outFile.write("{0:<20}".format("langevinTemp") + str(pD['Temperature (K)']) + "\n")
    
    outFile.write("\n")

    outFile.write("# Set the execution parameters\n\n")

    outFile.write("{0:<20}".format("outputname") + pD['Data File Name'] + "\n")
    outFile.write("{0:<20}".format("bincoordinates") + pD['Config File Name'] + ".restart.coor" + "\n")
    outFile.write("{0:<20}".format("binvelocities") + pD['Config File Name'] + ".restart.vel" + "\n")
    outFile.write("{0:<20}".format("run") + str(pD['Duration']) + "\n")

    outFile.close()
    return confFile
    
def linkFiles(paths, pD, wisdomFile, preheat = False):

    paramFile = "par_all27_prot_lipid.prm"
    
    print("################################################################\n" \
          + "Moving additional configuration files into " + paths['run'] + ".\n" \
          + "################################################################\n")
    
    # Make a link to the parameter file in the simulation run folder
    if not os.path.exists(paths['run'] + paramFile):
        os.link(paths['param'] + paramFile, paths['run'] + paramFile)
        
    # Make a link to the wisdom file in the simulation run folder
    if (not os.path.exists(paths['run'] + wisdomFile)) \
        and os.path.exists(paths['param'] + wisdomFile):
        os.link(paths['param'] + wisdomFile, paths['run'] + wisdomFile)
    
    # Make links to all the other configuration files in the simulation run folder
    if not os.path.exists(paths['run'] + pD['Config File Name'] + '.pdb'):
        os.link(my.configFile(paths['solvate'], pD['Config File Name'] + '.pdb'), paths['run'] + pD['Config File Name'] + '.pdb')
    if not os.path.exists(paths['run'] + pD['Config File Name'] + '.psf'):
        os.link(my.configFile(paths['solvate'], pD['Config File Name'] + '.psf'), paths['run'] + pD['Config File Name'] + '.psf')
    if not os.path.exists(paths['run'] + pD['Config File Name'] + '-restraint.pdb'):
        os.link(my.configFile(paths['restraint'], pD['Config File Name'] + '.pdb'), paths['run'] + pD['Config File Name'] + '-restraint.pdb')
    
    if pD['Run Type'] == 'New':
        if not os.path.exists(paths['run'] + pD['Config File Name'] + '-langevin.pdb') and pD['Restraint'] != 0:
            os.link(my.configFile(paths['temperature'], pD['Config File Name'] + '.pdb'), paths['run'] + pD['Config File Name'] + '-langevin.pdb')    
        if not preheat:
            if not os.path.exists(paths['run'] + pD['Config File Name'] + '-forcing.pdb'):
                os.link(my.configFile(paths['forcing'], pD['Config File Name'] + '.pdb'), paths['run'] + pD['Config File Name'] + '-forcing.pdb')
            if not os.path.exists(paths['run'] + pD['Config File Name'] + ".restart.coor"):
                os.link(paths['run'] + "Preheat/" + pD['Data File Name'] + ".restart.coor", paths['run'] + pD['Config File Name'] + ".restart.coor")
            if not os.path.exists(paths['run'] + pD['Config File Name'] + ".restart.vel"):
                os.link(paths['run'] + "Preheat/" + pD['Data File Name'] + ".restart.vel", paths['run'] + pD['Config File Name'] + ".restart.vel")
        else:
            if not os.path.exists(paths['run'] + pD['Config File Name'] + ".restart.coor"):
                os.link(paths['minimization'] + pD['Data File Name'] + ".restart.coor", paths['run'] + pD['Config File Name'] + ".restart.coor")
    
def runSim(executable, simPath, simFile, output, hostname):
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
    numThreads['nanotube'] = '+p16'
    numThreads['CGG-DPU-Laptop'] = '+p2'

    Namd2in=subprocess.Popen([executable, numThreads[hostname], '+isomalloc_sync', \
                simFile], stdin=subprocess.PIPE, stdout=logFile, stderr=logFile)
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

def getCNTBasis(inFile):
    """ getCNTBasis finds the basis of a nanotube with filename outFile. """

    # Opens the CNT prebond file and reads the header of the file

    with open(inFile) as basisFile:
        header = basisFile.readline()

    # Splits the first line of the CNT-prebond file, and finds the x, y, z basis vectors of the CNT 

    basis = re.split('\s+', header)
    xVec = eval(basis[1])
    yVec = eval(basis[2])
    zVec = eval(basis[3])
    
    return xVec, yVec, zVec
