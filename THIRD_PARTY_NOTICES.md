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

KLS can build METIS nested-dissection ordering from pinned submodules under
`third_party/metis` and `third_party/gklib`. METIS is Copyright 1997, Regents
of the University of Minnesota, and is licensed under Apache-2.0; see
`third_party/metis/LICENSE`. GKlib is Copyright 1995-2018, Regents of the
University of Minnesota. Its primary license is Apache-2.0 and it also carries
LGPL-2.1-or-later and BSD-3-Clause files; see `third_party/gklib/LICENSE.txt`
and `third_party/gklib/LICENSES.md`.
