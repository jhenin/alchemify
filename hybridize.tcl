
proc hybridize { psfA pdbA psfB pdbB } {

    set mA [mol new $psfA]
    mol addfile $pdbA
    set mB [mol new $psfB]
    mol addfile $pdbB
    set merged [::TopoTools::mergemols [list $mA $mB]]
    if { $merged == -1 } {
        puts "Error merging molecules"
        return -1
    }

    set offset [molinfo $mA get numatoms]
    set ntotal [molinfo $merged get numatoms]

    set sA [atomselect $merged "index < $offset"]
    set sB [atomselect $merged "index >= $offset"]
    set all [atomselect $merged "all"]

    #- Detect common atoms, which have the same: segname, (resname), resid, atom name, atom type, atom charge
    set A_ids [$sA get { segname resid name type }]
    set qA [$sA get charge]
    set B_ids [$sB get { segname resid name type }]
    set qB [$sB get charge]

    #- atoms classified as A, B, C (common), BC (copy of common atom in B)
    # with beta factor -1, 1, 0 , 2
    # Mark all atoms as 1 (B) to start with
    set beta [lrepeat $ntotal 1]

    foreach i [$sA list] qA_i $qA {
        # Look for atom in molecule B with same signature
        set match [lsearch $B_ids [lindex $A_ids $i]]
        set qB_i   0
        set q_tol  1e-2  ;# tolerance for "equal" charges
        if { $match > -1 } {
            set qB_i [lindex $qB $match]
        }

        if { $match > -1 && [expr abs($qA_i - $qB_i) < $q_tol] } {
            lset beta $i 0
            # matching atom in B marked as redundant
            lset beta [expr {$match + $offset}] 2
            # Remember equivalent atom to reconnect bonds, angles etc.
            set common_map([expr {$match + $offset}]) $i
        } else {
            lset beta $i -1
        }
    }
    # Any atom left at 1 at this point is unique to B
    $all set beta $beta
    $all delete

    #- check that atoms in A, B, AC have unique segname/resid/name, rename as necessary

    #-  TODO TODO
    # for any parameter coupling BC and B: replace BC atoms with equivalent C atom id
    # Detect all bonds etc. touching B, and reassign labels for C
    $sB delete
    set sB [atomselect $merged "beta 1"]
    reassign_struct $sB [array get common_map]
    $sB delete

    # Delete all BC atoms by creating new molecule with only the atoms we want to keep
    set keep [atomselect $merged "beta < 1.5"]
    set hybrid [::TopoTools::selections2mol $keep]
    if { $hybrid == -1 } {
        puts "Error extracting relevant atoms"
        return -1
    }
    set all_hyb [atomselect $hybrid all]
    $all_hyb set beta [$keep get beta]
    $keep delete
    $all_hyb delete
    mol delete $merged

    animate write psf "hybrid.psf" $hybrid
    animate write pdb "hybrid.pdb" $hybrid
    # Get proper coloring in VMD
    mol reanalyze $hybrid
    mol modcolor 0 $hybrid Beta
    mol modstyle 0 $hybrid CPK 1.000000 0.300000 12.000000 12.000000

    return $hybrid
}


proc reassign_struct { sel lmap } {
    # reassign structure involving selected atoms, replacing common atoms to be deleted (group B)
    array set map $lmap
    set mol [$sel molid]
    set Bids [$sel list]

    foreach type { bond angle dihedral improper crossterm } {
        if { $type == "bond" } {
            set list [topo get${type}list none]
        } else {
            set list [topo get${type}list]
        }
        foreach l $list {
            if { $type == "angle" || $type == "improper" || $type == "dihedral" } {
                # remove type entry
                set l [lrange $l 1 end]
            }
            set has_CB 0
            set has_B 0
            set newl [list]
            foreach a $l {
                if [info exists map($a)] {
                    lappend newl $map($a)
                    set has_CB 1
                } else {
                    lappend newl $a
                }
                if {[lsearch -integer -sorted $Bids $a] >= 0 } {
                    set has_B 1
                }
            }     
            if { $has_B && $has_CB } {
                puts "$type  $l  ->  $newl"
                topo -molid $mol del$type {*}$l
                topo -molid $mol add$type {*}$newl
            }
        }
    }
}
