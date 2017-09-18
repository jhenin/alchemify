# alchemify

An X-PLOR PSF post-processor for alchemical free energy calculations in [NAMD](http://www.ks.uiuc.edu/Research/namd)  

Written by Jérôme Hénin  

* * *

## Updates:

Version 1.4 (June 2009) - tiny fix to handle PSF-CMAP files properly.  

As of NAMD 2.7b2, alchemify is not needed anymore, since all unwanted terms are discarded by NAMD at runtime. At this point, however, this does not apply to NAMD binaries compiled with the memory optimized option.  

* * *

## Motivation

NAMD features an [implementation](http://www.ks.uiuc.edu/Research/namd/2.12/ug/node61.html) of alchemical FEP that uses the dual topology approach. Dual topology requires that the initial and final groups do not interact. To prevent such interactions, the PSF file describing the system should:  

*   contain no bonded parameter coupling the initial and final groups;
*   contain exclusion lists ensuring that no nonbonded interactions between these groups exist.

Psfgen is not designed to handle nonbonded exclusion lists, and it is not straightforward to use it for removing unwanted bonded terms. Alchemify is a small program that takes care of these two tasks with minimal effort from the user.  
For simulating single amino-acid mutations in proteins, the setup can be handled automatically by the [VMD](http://www.ks.uiuc.edu/Research/vmd) plugin [Mutator](http://www.ks.uiuc.edu/Research/vmd/plugins/mutator).

## Syntax

`alchemify <input PSF> <output PSF> <FEPfile> [FEP column]`

## How to use it  

*   Generate an X-PLOR PSF file containing the dual topology, using either psfgen, or CHARMM, or any suitable tool. This PSF may contain auto-generated angles and dihedrals. The unavoidable part of the work is to define the dual topology more or less by hand, unless the alchemical transformation affects only whole molecules and not parts of them.  

*   Create an FEPfile, i.e. a PDB file for your system with a column (usually the B column) containing a '1' flag for atoms that appear and '-1' for atoms that disappear.

*   Run `alchemify in.psf out.psf FEPfile`.

* * *

## License

Alchemify is Free Software. It is distributed under the [GNU General Public License, version 2](http://www.gnu.org/copyleft/gpl.html).
