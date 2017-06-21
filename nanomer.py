### Written by CGG, spring 2017
### Adapted from nanogen.py

# from __future__ import print_function
# from __future__ import division

from prody import *
import numpy as np
from scipy import *
import os
import socket
import datetime
import time
import shutil
import subprocess

import myutil as my

def run(pD):
    """ This is the main entry point for performing simulations with a carbon nanotube
    and a polymer. It takes a python dictionary of simulation parameters as its
    argument. Before calling run, you must call myutil.setPolymerParamDefaults()
    to create the dictionary."""
        
    # The paths dictionary holds the path name to each tier of the data hierarchy.

    paths = my.makePolymerPaths(pD)
    
    # Create a 'temporary' folder for the simulation. This will later be moved to
    # the correct place in the data hierarchy. So paths['run'] is the temporary
    # folder, while paths['data'] is the final destination.
    
    paths['run'] = os.getenv('HOME') + '/' + pD['Data Folder'] + '/'
    
    # If we're doing a production run (not a minimization) on nanotube, we want to
    # use the GPU, so we run the CUDA version of namd2. The CUDA version of namd seems
    # to give errors occasionally when minimizing, so we don't use it for that.
    
    hostname = socket.gethostname().partition('.')[0]
    
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

    # Make all the necessary pdb and psf files for the current simulation.

    makeConfigFiles(paths, pD)
    
    # Write out the NAMD configuration file.
    
    confFile = confWrite(paths, pD, hostname)
    
    # For a production run, set up the tcl forces to pull the polymer back and forth.
    
    if pD['Run Type'] == 'New':
    	tclForcesWrite(paths, pD)
    	
    # Link the various pdb, psf, parameter files, etc. to the temporary run folder.
    	
    linkFiles(paths, pD)
    
    # Actually run the simulation!
    
    result, runTime = runSim(executable, paths['run'], confFile, pD['Config File Name'], hostname)
    
    # Move the temporary folder into the correct place in the data hierarchy.
    
    shutil.move(paths['run'], paths['data'])
    
def makeConfigFiles(paths, pD):
    """ makeConfigFiles makes all the pdb and psf files necessary for the current
    simulation. It makes heavy use of the pdbWrite procedure in myutil.py. It
    does this by constructing a list of commands to be fed to VMD for each of 
    pdb and/or psf files that need to be created."""
    
    # This is the base name for all the various pdb and psf files we will generate.
    
    fileName = pD['Config File Name']

    # Create the desired nanotube. The files will be stored in the 'Config Files' folder
    # in the folder pointed to by paths['L'].
    
    package = "package require CNTtools 1.0\n"
    commands = [ "genNT " \
                + my.quoted(my.configFile(paths['L'], fileName)) + " " \
                + str(pD['L']/10) + " " + str(pD['n']) + " " + str(pD['m']) + "\n" ]
    commands.append("centerNT " + my.quoted(my.configFile(paths['L'], fileName)))
   
    my.pdbWrite(paths['L'], fileName, package, commands)
    
    # Join the nanotube to the polymer in single pdb and psf files. The files will be
    # stored in the 'Config Files' folder in the folder pointed to by paths['polymer'].
    
    package = "package require CNTPoly 1.0\n"
    commands = [ "joinMolecules " \
                + my.quoted(my.configFile(paths['L'], fileName)) + " " \
                + my.quoted(my.configFile(paths['L'], fileName)) + " " \
                + my.quoted(paths['pdb'] + pD['Polymer PSF']) + " " \
                + my.quoted(paths['pdb'] + pD['Polymer PDB']) + " " \
                + my.quoted(my.configFile(paths['polymer'], fileName)) + "\n" ]
    
    my.pdbWrite(paths['polymer'], fileName, package, commands)
    
    # Write out the pdb file with the restraint information for the nanotube. The
    # resulting pdb file will be written in the folder pointed to by paths['restraint'].
    
    package = "package require CNTtools 1.0\n"
    commands = [ "NTrestraint " \
                + my.quoted(my.configFile(paths['polymer'], fileName)) + " " \
                + my.quoted(my.configFile(paths['restraint'], fileName)) + " " \
                + str(pD['Restraint']) + "\n" ]
    
    my.pdbWrite(paths['restraint'], pD['Config File Name'], package, commands)

def confWrite(paths, pD, hostname):
    """ confWrite generates a .conf file to use as input to namd2. To organize simulations
    with different parameters, confWrite will create a directory for the simulation."""

    # paths['run'] contains the name of the temporary folder where the simulation will
    # be run. The first thing we need to do is set the configuration file name.
    
    confFile = paths['run'] + pD['Config File Name'] + '.conf'

    # Create the temporary run directory if it doesn't already exist.
    
    if not os.path.exists(paths['run']):
        os.makedirs(paths['run'])

    # Tell the world what we're doing.
    
    print("################################################################\n" \
          + "Writing NAMD configuration file " + confFile + ".\n" \
          + "################################################################\n")
    
    # Open the configuration file for writing.
    
    outFile = open(confFile, "w")
    
    # Mark the configuration file with the host, date, and time of creation.
    
    outFile.write("# Configuration file written on host " + hostname + "\n" \
        + "# on " + str(datetime.datetime.now().date()) \
        + " at " + str(datetime.datetime.now().time()).partition('.')[0] + ".\n\n")
        
    # Set the frequency at which energy is written to the log file.
        
    outFile.write("# Limit the length of the log file\n\n")
    
    outFile.write("{0:<20}".format("outputEnergies") + str(100*pD['outputFreq']) + "\n\n")
    
    # Set up the cylindrical boundary conditions
    
    radius = 4*2.46*np.sqrt(pD['n']**2 + pD['n']*pD['m'] + pD['m']**2)/3.1415926
    length = 2*pD['L'] + 25
    
    outFile.write("# Set up cylindrical boundary conditions\n\n")
    
    outFile.write("{0:<20}".format("cylindricalBC") + "on" + "\n")
    outFile.write("{0:<20}".format("cylindricalBCCenter") + "0 0 0" + "\n")
    outFile.write("{0:<20}".format("cylindricalBCAxis") + "z" + "\n")
    outFile.write("{0:<20}".format("cylindricalBCr1") + "{0:<8.3f}\n".format(radius))
    outFile.write("{0:<20}".format("cylindricalBCl1") + "{0:<10.3f}\n".format(length))
    outFile.write("{0:<20}".format("cylindricalBCk1") + "1.0" + "\n")
    
    outFile.write("\n")
    
    # Set up the input pdb and psf files.

    outFile.write("# Set the input files\n\n")
    
    outFile.write("{0:<20}".format("structure") + pD['Config File Name'] + ".psf" + "\n")
    outFile.write("{0:<20}".format("coordinates") + pD['Config File Name'] + ".pdb" + "\n")
    outFile.write("\n")
    
    # This section sets the force field parameters and the parameter files.
    
    outFile.write("# Set the force field parameters\n\n")
    
    outFile.write("{0:<20}".format("paraTypeCharmm") + "on" + "\n")
    outFile.write("{0:<20}".format("parameters") + "par_all27_prot_lipid.prm" + "\n") 
    outFile.write("{0:<20}".format("parameters") + "par_all35_ethers.inp" + "\n") 
    outFile.write("{0:<20}".format("exclude") + "scaled1-4" + "\n")
    outFile.write("{0:<20}".format("cutoff") + "12.0" + "\n")
    outFile.write("{0:<20}".format("pairlistdist") + "14.0" + "\n")
    outFile.write("{0:<20}".format("switching") + "on" + "\n")
    outFile.write("{0:<20}".format("switchdist") + "10.0" + "\n")
    outFile.write("{0:<20}".format("PME") + "no" + "\n")
    outFile.write("\n")
    
    # Some miscellaneous integration parameters.

    outFile.write("# Set the integration parameters\n\n")
    
    outFile.write("{0:<20}".format("timestep") + str(pD['dt (fs)']) + "\n")
    outFile.write("{0:<20}".format("nonbondedFreq") + "2" + "\n")
    outFile.write("{0:<20}".format("fullElectFrequency") + "4" + "\n")
    outFile.write("{0:<20}".format("stepspercycle") + "20" + "\n")
    outFile.write("{0:<20}".format("rigidBonds") + "water" + "\n")
    # outFile.write("{0:<20}".format("FullDirect") + "yes" + "\n")

    outFile.write("\n")
    
    # If this is a new run (not a minimization), we set up the tclForces file to pull
    # the polymer back and forth.

    if pD['Run Type'] == 'New':
    
        outFile.write("# Turn on tcl forces\n\n")
    
        outFile.write("{0:<20}".format("tclForces") + "on" + "\n")
        outFile.write("{0:<20}".format("tclForcesScript") + pD['tcl Script'] + "\n")

        outFile.write("\n")
        
    # If pD['Restraint'] is zero, we're using fixed atoms for the nanotube. Otherwise,
    # we're using harmonic constraints on the nanotube.

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
    
    # Now we set the frequency of the dcd files.

    outFile.write("# Set the frequency of the various dcd files\n\n")
    
    outFile.write("{0:<20}".format("restartfreq") + str(pD['outputFreq']) + "\n")
    outFile.write("{0:<20}".format("dcdfreq") + str(pD['outputFreq']) + "\n")
    outFile.write("{0:<20}".format("veldcdfreq") + str(pD['outputFreq']) + "\n")
    outFile.write("{0:<20}".format("forcedcdfreq") + str(pD['outputFreq']) + "\n")
    
    outFile.write("\n")
    
    # If this is a minimization run, we'll set the output file name, the temperature
    # (not used, but required by NAMD), and add a minimize command. That's all we need
    # to do, so we return.
    
    if pD['Run Type'] == 'Minimize':
    
        outFile.write("{0:<20}".format("outputname") + pD['Data File Name'] + "\n")
        outFile.write("{0:<20}".format("temperature") + '0' + "\n")
        outFile.write("{0:<20}".format("minimize") + str(pD['Min Duration']) + "\n")
        
        outFile.close()
        return confFile
        
    # If this is a production run, we might turn on a thermostat for the nanotube. Note
    # that right now there's no provision for using a thermostat in polymer runs. This
    # is currently left over from the water runs. We'll need to add that back in if we
    # really want to do polymer runs with a thermostat on.
        
    outFile.write("# Set the Langevin thermostat\n\n")

    if pD['Thermostat'] == 'Off':
        outFile.write("{0:<20}".format("langevin") + "off" + "\n")
    else:
        outFile.write("{0:<20}".format("langevin") + "on" + "\n")
        outFile.write("{0:<20}".format("langevinFile") + pD['Config File Name'] + "-langevin.pdb" + "\n")
        outFile.write("{0:<20}".format("langevinCol") + "O" + "\n")
        outFile.write("{0:<20}".format("langevinTemp") + str(pD['Temperature (K)']) + "\n")
    
    outFile.write("\n")
    
    # Finally, we set the output file name, the temperature, the coordinate file (from
    # a minimization run), and the duration of the run.

    outFile.write("# Set the execution parameters\n\n")

    outFile.write("{0:<20}".format("outputname") + pD['Data File Name'] + "\n")
    outFile.write("{0:<20}".format("temperature") + "{0:<8.2f}\n".format(pD['Temperature (K)']))
    outFile.write("{0:<20}".format("bincoordinates") + pD['Config File Name'] + ".restart.coor" + "\n")
    outFile.write("{0:<20}".format("run") + str(pD['Duration']) + "\n")

    outFile.close()
    return confFile
    
def tclForcesWrite(paths, pD):

    # We create the tcl forces script file in the temporary run directory.
    
    tclFile = paths['run'] + pD['tcl Script']
    if not os.path.exists(paths['run']):
        os.makedirs(paths['run'])

    print("################################################################\n" \
          + "Writing tcl Forces file " + tclFile + ".\n" \
          + "################################################################\n")
          
    # We make the assumption that we've already minimized the system. We need to
    # pull the coordinates for the first carbon atom in the polymer from the 
    # minimized dcd file using ProDy.
    
    x, y, z = getAtomCoords(paths['minimization'], pD, 'C1')
    
    # Open the tcl forces file.
    
    outFile = open(tclFile, "w")
    
    # Mark the file with the date and time.
    
    outFile.write("# tclForces file written" + "\n" \
        + "# on " + str(datetime.datetime.now().date()) \
        + " at " + str(datetime.datetime.now().time()).partition('.')[0] + ".\n\n")
        
    # The rest of this is just the tcl commands necessary to pull the end carbon
    # atom back and forth with a specified amplitude and frequency.
        
    outFile.write('set a [ atomid PO1 1 C1 ]\n')
    outFile.write("addatom $a\n")
    outFile.write("set x0 " + str(x) + "\n") 
    outFile.write("set y0 " + str(y) + "\n") 
    outFile.write("set z0 " + str(z) + "\n\n")
    outFile.write("set k 600\n")
    outFile.write("set f {0:<8.6f}\n".format(1.0/(pD['Driving Period (ps)']*1000)))
    outFile.write("set amp " + str(pD['Driving Amplitude (A)']) + "\n")
    outFile.write("set pi 3.1415926\n") 
    
    outFile.write("\n")
    
    outFile.write("proc calcforces {} {\n\n")
    
    outFile.write("    global a x0 y0 z0 k f amp pi\n\n")
    
    outFile.write("    loadcoords p\n\n")
    outFile.write("    set t [ getstep ]\n\n")
    outFile.write("    set z [ expr $z0 + $amp*sin(2*$pi*$f*$t) ]\n")
    outFile.write('    set r " $x0 $y0 $z "\n')
    outFile.write("    set diff [ vecsub $p($a) $r ]\n")
    outFile.write("    set force [ vecscale [ expr -$k ] $diff ]\n")
    
    outFile.write("    addforce $a $force\n")
    
    outFile.write("\n}\n")
        
    outFile.close()  

def linkFiles(paths, pD):

    paramFile1 = "par_all27_prot_lipid.prm"
    paramFile2 = "par_all35_ethers.inp"

    print("################################################################\n" \
          + "Moving additional configuration files into " + paths['run'] + ".\n" \
          + "################################################################\n")
    
    # Make a link to the parameter file(s) in the simulation run folder.
    
    if not os.path.exists(paths['run'] + paramFile1):
        os.link(paths['param'] + paramFile1, paths['run'] + paramFile1)

    if not os.path.exists(paths['run'] + paramFile2):
        os.link(paths['param'] + paramFile2, paths['run'] + paramFile2)

    # Make links to the pdb and psf files for the run.

    if not os.path.exists(paths['run'] + pD['Config File Name'] + '.pdb'):
        os.link(my.configFile(paths['polymer'], pD['Config File Name'] + '.pdb'), paths['run'] + pD['Config File Name'] + '.pdb')
    if not os.path.exists(paths['run'] + pD['Config File Name'] + '.psf'):
        os.link(my.configFile(paths['polymer'], pD['Config File Name'] + '.psf'), paths['run'] + pD['Config File Name'] + '.psf')
        
    # Make a link to the restraint file for the run.

    if not os.path.exists(paths['run'] + pD['Config File Name'] + '-restraint.pdb'):
        os.link(my.configFile(paths['restraint'], pD['Config File Name'] + '.pdb'), paths['run'] + pD['Config File Name'] + '-restraint.pdb')
        
    # Make a link to the input coordinate files for the run.

    if pD['Run Type'] == 'New':
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

def getAtomCoords(path, pD, atomName):

    pdbFile = path + pD['Config File Name'] + '.pdb'
    dcdFile = path + pD['Data File Name'] + '.dcd'

    # Now, open up the DCD file and get the ensembles for the carbons and oxygens

    structure = parsePDB(pdbFile)
    coordEnsemble = parseDCD(dcdFile)

    atom = structure.select("name " + atomName)
    
    coordEnsemble.setAtoms(atom)
    atomC = coordEnsemble.getCoordsets()

    coordX = atomC[-1, 0, 0]
    coordY = atomC[-1, 0, 1]
    coordZ = atomC[-1, 0, 2]
    
    return coordX, coordY, coordZ
