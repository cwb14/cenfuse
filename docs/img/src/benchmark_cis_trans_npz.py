"""scipy.sparse.load_npz without numpy's slow streamed read of zip members (same arrays, same matrix): each member is
read whole with zipfile and parsed with np.load from memory. Verified equal to load_npz on the tomato matrix (memo.md)."""
import io
import zipfile

import numpy as np
import scipy.sparse as sp


def load_npz(path):
    with zipfile.ZipFile(path) as z:
        a = {n[:-4]: np.load(io.BytesIO(z.read(n)), allow_pickle=False) for n in z.namelist()}
    fmt = a["format"].item()
    fmt = fmt.decode() if isinstance(fmt, bytes) else fmt
    shape = tuple(int(x) for x in a["shape"])
    if fmt == "csr":
        return sp.csr_matrix((a["data"], a["indices"], a["indptr"]), shape=shape)
    if fmt == "csc":
        return sp.csc_matrix((a["data"], a["indices"], a["indptr"]), shape=shape)
    if fmt == "coo":
        return sp.coo_matrix((a["data"], (a["row"], a["col"])), shape=shape)
    raise ValueError(f"unsupported sparse format {fmt}")
