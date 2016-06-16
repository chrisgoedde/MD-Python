### Written by Matthew Kwiecien Jul, 2015
### Modified by CGG spring and summer 2016
package provide CNTtools 1.0
package require inorganicbuilder
package require nanotube
package require pbctools
package require psfgen

proc genNT {molnm path l n m} {
    # Create a nantube of length l with nxm structure
    nanotube -l $l -n $n -m $m
    set mymol [atomselect top all]
    
    # Move the nanotube so its CM is at the origin
    set cxyz [ measure center $mymol ] ; # find the coordinates of the center of the nanotube
    set ncxyz {} ; # make an empty list to hold the negative of the nanotube center
    foreach i $cxyz { lappend ncxyz [ expr { -$i } ] } ; # set ncxyz to -$cxyz
    $mymol moveby $ncxyz ; # move the nanotube so its CM is at the origin

    # Save the nanotube in the file molnm in the folder specified by path
    set molpath ${path}${molnm}
    $mymol writepsf $molpath.psf
    $mymol writepdb $molpath.pdb

}

proc fixNT {molnm} {
    # Open the files molnm.psf and molnm.pdb to load the molecule
    mol new [file normalize ${molnm}.psf] type psf autobonds off waitfor all
    mol addfile [file normalize ${molnm}.pdb] type pdb autobonds off waitfor all

    set all0 [atomselect top all]

    # Let nanotube be flexible by default; this will be changed later if necessary
    $all0 set beta 0
    
    # Set the resid of the nanotube to 1
    $all0 set resid 1

    # Write out the molecule back to file molnm
    set all [atomselect top all]
    $all writepsf $molnm.psf
    $all writepdb $molnm.pdb

}

proc removeLangevinWater {molnm} {
    #removes langevin thermostat for water molecules
    
    # Open the files molnm.psf and molnm.pdb to load the molecule
    mol new [file normalize ${molnm}.psf] type psf autobonds off waitfor all
    mol addfile [file normalize ${molnm}.pdb] type pdb autobonds off waitfor all

    # Set the occupancy column to zero
    set all0 [atomselect top all]
    $all0 set occupancy 0.00

    # Write the file back out to molnm
    set all [atomselect top all]
    $all writepsf $molnm.psf
    $all writepdb $molnm.pdb
}

proc pbcNT {molnm fileOut ntype} {
    # Based on procedure from Tom Sisan from Northwestern University
    # Modified by CGG Spring/Summer 2016
    
    # Read in cnt that does not have periodic bonds.
    mol new [file normalize ${molnm}.psf] type psf autobonds off waitfor all
    mol addfile [file normalize ${molnm}.pdb] type pdb autobonds off waitfor all
    set molid [molinfo top]

    # Get basis vectors (a,b,c) and origin (o).
    set basis [::inorganicBuilder::findBasisVectors $molid]
    foreach {o a b c} $basis {}

    # Set cutoff distance for finding bonded atoms.
    set CNTbondCutoff 1.6
    switch $ntype {
        bnt {set scalent 1.03}
        snt {set scalent 1.25}
        default {set scalent 1.0}
    }
    set CNTbondCutoff [expr $CNTbondCutoff*$scalent]

    # Define a system box.
    set mybox [::inorganicBuilder::defineMaterialBox $o [list $a $b $c] $CNTbondCutoff]

    # This was probably intended to be used for reading in molecules with the "filebonds off" setting.
    # It currently appears to be unused.
    # set fbonds off

    # Set the periodic box to the system box, and wrap atom coordinates into the box.
    # Any atoms previously outside the box will now be in the box
    ::inorganicBuilder::setVMDPeriodicBox $mybox
    ::inorganicBuilder::transformCoordsToBox $mybox $molid

    # Recalculate bonds for selected molecule (possibly not necessary because of later command?).
    mol bondsrecalc $molid

    # Set which dimensions are periodic (x,y,z).
    set periodicIn {false false true}
    
    # Build the spring bonds.
    set relabelBonds 0
    ::inorganicBuilder::buildBonds $mybox $relabelBonds $periodicIn $molid

    # Write pdb and psf files, with just linear bonds so far.
    set fname "${fileOut}"
    set fname0 ${fname}-prebond
    set mymol [atomselect top all]
    $mymol writepsf $fname0.psf
    $mymol writepdb $fname0.pdb

    # Build the angular bonds into output file fname.
    set dihedrals 1
    ::inorganicBuilder::buildAnglesDihedrals $fname0 $fname $dihedrals

    # Added this code to fix up the pbc box in the pdb file
    mol new [file normalize $fname.psf] type psf autobonds off waitfor all
    mol addfile [file normalize $fname.pdb] type pdb autobonds off waitfor all

    set mymol [ atomselect top all ]

    set cell [ lindex [ pbc get ] 0 ]
    set X [ lindex $a 0 ]
    set Y [ lindex $b 1 ]
    set Z [ lindex $c 2 ]
    
    lset cell 0 $X
    lset cell 1 $Y
    lset cell 2 $Z
    
    pbc set $cell
    
    set newCell [ lindex [ pbc get ] 0 ]
    
    $mymol writepdb $fname.pdb
    $mymol writepsf $fname.psf

}
