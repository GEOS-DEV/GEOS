from spack.package import *
from spack_repo.builtin.packages.hypre.package import Hypre as BuiltinHypre


class Hypre(BuiltinHypre):
    # HYPRE 3.1.0's MGR column-lumped restriction path copies the data out of
    # B_CF but never destroys that DenseBlockMatrix.  The sanitizer integrated
    # tests exercise this path through HypreDrive, so keep the fix local until
    # it is included in the upstream HYPRE package.
    patch("hypre_mgr_col_lumped_destroy.patch", when="@3.1.0")
