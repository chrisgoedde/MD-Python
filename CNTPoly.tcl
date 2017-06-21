package provide CNTPoly 1.0

package require psfgen
package require mergestructs

proc joinMolecules { firstPSF firstPDB secondPSF secondPDB output } {

    resetpsf
    
    readpsf $firstPSF.psf
    coordpdb $firstPDB.pdb
    readpsf $secondPSF.psf
    coordpdb $secondPDB.pdb
    writepsf $output.psf
    writepdb $output.pdb

}

proc mergeMolecules { first second } {

    set a [ list $first.psf $first.pdb ]
    set b [ list $second.psf $second.pdb ]
    ::MergeStructs::mergeMolecules $a $b $first-$second-merged
    
}

proc getAtomCoords { dcdFile psfFile frameNumber } {

    mol new [file normalize ${dcdFile}.dcd] type dcd autobonds off waitfor all
    mol addfile [file normalize ${psfFile}.psf] type psf autobonds off waitfor all
    
    set a [ atomselect top "segname PO1 and name C1" ]
    $a frame frameNumber
    
    set x $a get x
    set y $a get y
    set z $a get z
    
}