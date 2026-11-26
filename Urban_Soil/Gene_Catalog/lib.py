#!/usr/bin/env python3
"""
Vendored copy of the small helper module shared with the Shanghai
Dogs gene-catalogue pipeline.

The Urban Soil scripts in this directory normally import these two
functions directly from:

    /work/microbiome/shanghai_dogs/resource_generation/lib.py

(via `sys.path.insert(0, ...)`) so that both catalogues stay in sync.
This file is kept here purely as a fallback/reference in case that
shared path is ever unavailable -- if you use it, just remove the
`sys.path.insert` + `from lib import ...` lines in the calling script
and drop this file next to it instead.

NOTE: `pad9()` here returns a plain zero-padded 9-digit ID
(e.g. US.95NT.000000001), which is what 05_postprocess_mmseqs95.py
and 07_US95_emapper.py actually use. 03_Run_US_100NT_catalog.py
defines its *own* local pad9() that groups digits with underscores
(e.g. US.100NT.000_000_001) instead of importing this one. Both ID
formats appear in the final US.clusters.tsv.xz cluster-mapping file,
which is why 06_make_US_cluster_table.py's numeric_id() parser
strips underscores before converting to int -- this is intentional,
not a bug, but worth knowing if you touch ID formatting.
"""

import lzma


def pad9(prefix, ix):
    return f"{prefix}.{ix:09d}"


def xz_out(path):
    """Return a writable, text-free (binary) lzma file handle.

    Use as: `with xz_out(path) as f: f.write(b"...")`.
    """
    return lzma.open(path, "wb")
