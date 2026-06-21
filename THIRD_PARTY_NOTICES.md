# Third-Party Notices

KLS vendors SuiteSparse-derived KLU, BTF, AMD, COLAMD, and UFconfig sources
from the Trilinos source tree under `third_party/suitesparse`. These sources
are built into KLS as the current in-tree serial sparse direct solver engine.

The vendored SuiteSparse-derived sources are Copyright (c) 2004-2007,
University of Florida and their listed authors. See the per-package
`README.txt` and source headers under `third_party/suitesparse` for detailed
notices. SPDX-License-Identifier: LGPL-2.1-or-later for the LGPL-covered
packages; the license text is in `LICENSE` and `LICENSES/LGPL-2.1-or-later.txt`.
`UFconfig` states no licensing restrictions in its header.

KLS does not include CKTSO source or binaries. CKTSO is used only as an
optional external benchmark when the user provides a licensed local
installation.
