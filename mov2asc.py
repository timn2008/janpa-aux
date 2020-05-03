# A Python replacement for NwChem's mov2asc program needed to convert NwChem's
#
# The structure of .movecs file seems to consist of records, each being:
#  record = [4-bytes-integer-size] [data bytes] [4-bytes-integer-size]
#
# (c) Tymofii Nikolaienko, 2017
#
# This file is a part of the JANPA project. 
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 
# 3. All advertising and/or published materials mentioning features or use 
#    of this software must display the following acknowledgement:
# 
#       This product includes components from JANPA package of programs
#       ( http://janpa.sourceforge.net/ ) developed by Tymofii Nikolaienko
# 
# 4. Neither the name of the developer, Tymofii Nikolaienko,  nor the
#    names of its contributors may be used to endorse or promote products
#    derived from this software without specific prior written permission.
# 
# 5. In case if the code of JANPA package code and/or its parts and/or any data 
#    produced with JANPA package of programs are published, the following citation 
#    for the JANPA package of programs should be given:
# 
#     T.Y.Nikolaienko, L.A.Bulavin, D.M.Hovorun; Comput.Theor.Chem. (2014)
#     V. 1050, P. 15-22, DOI: 10.1016/j.comptc.2014.10.002
# 
# THIS SOFTWARE IS PROVIDED BY ''AS IS'' AND ANY EXPRESS OR IMPLIED 
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
# IN NO EVENT SHALL Tymofii Nikolaienko BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 

import sys
import struct



def read_rec(doPrint = False, convert = None):
    """ reads one record from .movecs file (global variable f) and converts 
        the read bytes if necessary """        
    global f
    sz1 = struct.unpack("<L", f.read(4))[0]
    msg = f.read(sz1)
    sz2 = struct.unpack("<L", f.read(4))[0]
    if doPrint:
        print '%d -> "%s" <- %d' % (sz1, msg, sz2)
        print
    if convert is None:
        return msg
    else:
        return struct.unpack(convert, msg)[0]
    
def print_double_arr():
    """ Reads a series of doubles, determining the number of values to
        read as data_size/8 """
    double_arr = read_rec()
    i = 0
    while i < len(double_arr):
        print "%25.15E" % ( struct.unpack("<d", double_arr[i:i+8])[0] ),
        i += 8
        if (i % 24) == 0:
            print
    if (i % 24) != 0:
        print



if len(sys.argv) < 2:
    print "Usage: python mov2asc.py filename.movecs"
    quit()

        
fname = sys.argv[1]

with open(fname, 'rb') as f:    
    print "# This is an NWChem movecs file translated by mov2asc"
    
    
    id = read_rec()
    i = 0
    for _ in range(3):      # basissum, geomsum, bqsum
        print id[i:i+32]  
        i += 32        
    print id[i:i+20] # scftype20
    print id[i+20:]  # date
    
    
    method = read_rec()  # scftype20
    print method
    
    print "%10d" % read_rec(convert = "<Q")  # lentit
    
    title = read_rec()  # title 
    print title
    
    print "%10d" % read_rec(convert = "<Q")  # lenbas
        
    sect = read_rec()
    print sect   # basis_name
        
    nsets = read_rec(convert = "<Q")          # nsets
    print "%10d" % nsets   
    nBasis = read_rec(convert = "<Q")         # nbf
    print "%10d" % nBasis
    nMOs = []
    for jset in range(nsets):
        nMOs.append(read_rec(convert = "<Q")) # nmo(i) -- perhaps, the number if MOs in each set
        print "%10d" % nMOs[-1]   
    
    for jset in range(nsets):
        print_double_arr()   # Occupation numbers
        print_double_arr()   # Eigenvalues
        for _ in range(nBasis):
            print_double_arr()
    
    # additional two doubles
    print_double_arr()     # energy, enrep
    
    #print '-'*80
    #sz = struct.unpack("<L", f.read(4))[0]
    #print 'next: ',sz
    