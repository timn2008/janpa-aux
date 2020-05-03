# A simple program for parsing the results of CLPO analysis and preparing the JMOL script
# able to visualize covalent bonds
# Part of JANPA project: http://janpa.sf.net
# Cite as: [ Int. J. Quantum Chem. (2018), e25798, DOI: 10.1002/qua.25798 ]
#
# (c) Tymofii Nikolaienko, 2018
#

from __future__ import print_function
import re
import numpy as np
import sys
    
    
def parse_Janpa(fname):
    txt = open(fname).read()

    clpoRes = re.findall(r"\*\*\* Summary of CLPO results\s+(.*?)^$", txt, re.MULTILINE | re.DOTALL)
    if len(clpoRes) != 1:
        pritn ('ERROR2:', fname)
    bonds = re.findall(r"\(BD\)\s+([A-Za-z]+)(\d+)\-([A-Za-z]+)(\d+).*?Io = ([\d\.]+)\s+([\d.]+)\s+h(\d+)@([A-Za-z\d]+).*?\(([ \d+\-.]+)\)\s*\+\s*h(\d+)@([A-Za-z\d]+).*?\(([ \d+\-.]+)\)", clpoRes[0])
    # '15   (BD) BR1-H2, Io = 0.2133                1.96428         h15@BR1 * ( 0.7789) + h51@H2 * ( 0.6272)' ->
    #   -> ('BR', '1', 'H', '2', '0.2133', '1.96428', '15', 'BR1', ' 0.7789', '51', 'H2', ' 0.6272')
    int_bonds = []
    for (a1, I, a2, J, io, occup, h1, _, cA, h2, _, cB) in bonds:
        int_bonds.append([int(I), int(J)])

    max_atom_num = 0
    for ib in int_bonds:
        max_atom_num = max(max_atom_num, ib[0], ib[1])
        
    # convert that into adjacency matrix    
    adjMatr = np.zeros((max_atom_num+1,max_atom_num+1), dtype=int) # use 1-based atomic numbers for simplicity
    for ib in int_bonds:
        adjMatr[ib[0]][ib[1]] += 1
        adjMatr[ib[1]][ib[0]] += 1
    #print (adjMatr)

    print('// JMol script produced by %s from: %s' % (sys.argv[0], fname))
    print('connect none')
    for i, row in enumerate(adjMatr):
        for k in range(i+1, len(row)):
            if row[k] > 0:
                assert row[k]<=3, 'Bonds of multiplicities higher than 3 are not supported!'
                print('connect (atomno=%d) (atomno=%d) %s ' % (i,k, ['single', 'double', 'triple'][row[k]-1] ))#, end='')
    print() 
        

if len(sys.argv) < 2:
    exit(1)

parse_Janpa(sys.argv[1]) 
                
