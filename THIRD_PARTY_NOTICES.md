# Third-Party Notices

KLS vendors SuiteSparse-derived KLU, BTF, AMD, CAMD, COLAMD, CCOLAMD, and
UFconfig sources from the Trilinos source tree under `third_party/suitesparse`.
These sources are built into KLS as the current in-tree sparse direct solver
engine and ordering toolkit.

The vendored SuiteSparse-derived sources are Copyright (c) 2004-2007,
University of Florida and their listed authors. See the per-package
`README.txt` and source headers under `third_party/suitesparse` for detailed
notices. SPDX-License-Identifier: LGPL-2.1-or-later for the LGPL-covered
packages; the license text is in `LICENSE` and `LICENSES/LGPL-2.1-or-later.txt`.
`UFconfig` states no licensing restrictions in its header.

KLS does not include CKTSO source or binaries. CKTSO is used only as an
optional external benchmark when the user provides a licensed local
installation.

KLS does not vendor HSL MC64 or solver-tree MC64 copies that retain HSL
redistribution restrictions. The bounded exact static-pivot assignment code in
`src/kls.c` is KLS in-tree code under the project LGPL-2.1-or-later license.
KLS builds BSD-licensed SPRAL Hungarian/auction matching and scaling support
from the pinned `third_party/spral` submodule by default, or can link to a
compatible system SPRAL installation after the builder explicitly acknowledges
that the selected system library is redistributable with LGPL-2.1-or-later KLS.
The bundled KLS build uses only SPRAL's scaling subset: `src/matrix_util.f90`,
`src/scaling.f90`, and
`interfaces/C/scaling.f90`, plus `include/spral_scaling.h`. Builds configured
with `KLS_ENABLE_SPRAL_SCALING=OFF` omit this component.

SPRAL is Copyright (c) 2014-2025, The Science and Technology Facilities
Council (STFC), and is licensed under BSD-3-Clause; see
`third_party/spral/LICENCE`.

KLS can build METIS nested-dissection ordering from pinned submodules under
`third_party/metis` and `third_party/gklib`. METIS is Copyright 1997, Regents
of the University of Minnesota, and is licensed under Apache-2.0; see
`third_party/metis/LICENSE`. GKlib is Copyright 1995-2018, Regents of the
University of Minnesota. Its primary license is Apache-2.0 and it also carries
LGPL-2.1-or-later and BSD-3-Clause files; see `third_party/gklib/LICENSE.txt`
and `third_party/gklib/LICENSES.md`.

KLS can build SCOTCH nested-dissection ordering from the pinned submodule under
`third_party/scotch`. SCOTCH is Copyright 2004-2021 by the listed SCOTCH
contributors and is licensed under CeCILL-C; see
`third_party/scotch/LICENSE_en.txt` and
`third_party/scotch/doc/CeCILL-C_V1-en.txt`.
