### Written by Matthew Kwiecien Jul, 2015
package provide CNTtools 1.0
package require inorganicbuilder
package require nanotube
package require pbctools
package require psfgen

proc genNT {molnm path l n m} {
  # Create a nantube of length l with nxm structure
  nanotube -l $l -n $n -m $m
  set mymol [atomselect top all]

  # saves nanotube in new directory
  set molpath ${path}${molnm}
  $mymol writepsf $molpath.psf
  $mymol writepdb $molpath.pdb

}

proc fixNT {molnm} {

  mol new [file normalize ${molnm}.psf] type psf autobonds off waitfor all
  mol addfile [file normalize ${molnm}.pdb] type pdb autobonds off waitfor all

  set all0 [atomselect top all]

  # Let nanotube be flexible by default; this will be changed later if necessary
  $all0 set beta 0
  
  # Set the resid of the nanotube to 1

  $all0 set resid 1

  set all [atomselect top all]
  $all writepsf $molnm.psf
  $all writepdb $molnm.pdb

}

proc removeLangevinWater {molnm} {
  #removes langevin thermostat for water molecules
  mol new [file normalize ${molnm}.psf] type psf autobonds off waitfor all
  mol addfile [file normalize ${molnm}.pdb] type pdb autobonds off waitfor all

  set all0 [atomselect top all]
  $all0 set occupancy 0.00

  set all [atomselect top all]
  $all writepsf $molnm.psf
  $all writepdb $molnm.pdb
}

proc pbcNT {molnm fileOut ntype} {
  ## Credit to Tom Sisan from Northwestern University
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

  # # Not sure why the below is here. Maybe I use it in the calling space.
  # ::inorganicBuilder::setVMDPeriodicBox $mybox
  # return $fname

}
