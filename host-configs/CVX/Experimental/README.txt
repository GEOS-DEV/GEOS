# --------------------------------------------------------------------------------------- #
# Chevron Experimental .cmake files
#
# michael.thomadakis@chevron.com HPC Innovation and R&D
#
# .../host-configs/CVX/Experimental
#
# Contains a number of .cmake files that are working but have been plaved here in order
# to avoid clutter in CVX's cmake files 
# Trilinos/ contains cmake files that build with Trilinos but not Hypre
# NVHPC/ contains cmake files that attempted to build using NVidia's HPC SDK stack 
#   currently cannot build som TPL
# IntelCompilers/ contains cmake files that use Intel classic or OneAPI 
#   currently cannot build som TPL

To use any of the cmake files, copy it to the directory above (Experimental/..) and
follow the standard process to build TPL/GEOS via Chevron's framework
