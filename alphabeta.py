from moldenio import *
import numpy as np
import re


def load_matr(fname):
    # Get the total number of NAOs
    with open(fname) as f:
        print (next(f)) # skip the first line
        nNAOs = int(  next(f).strip().split("\t")[0] )
        NAO_labels = next(f).strip().split("\t")

    # Now load the S.D.S matrix in NAO basis created by JANPA:
    return np.loadtxt(fname, skiprows=3, usecols=range(nNAOs))
### 


fname = '1.m'
print '-'*80
Va, Vb, na, nb, Ea, Eb, atom_BF_ids, hdr_lines, geom =  load_molden(fname)
print '-'*80

IPython.embed()


save2molden(fname+'_Alpha.molden', Va, hdr_lines, Ea, 2. * na)

save2molden(fname+'_Beta.molden', Vb, hdr_lines, Eb, 2. * nb)
