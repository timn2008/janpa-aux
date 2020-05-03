################################################################
# Loads basis set assignment and MO data from MOLDEN file
# (c) Tymofii Nikolaienko, 2020
################################################################


import numpy as np
import re
import IPython


def load_molden(fname):
    current_section = None
    current_atom_Id = None
    num_L_funcs = {'s': 1, 'p': 3, 'd': 6, 'f': 10, 'g': 15} # default is to use Cartesian functions
    atomicBFs = {}
    MOs = []
    current_MO = None
    save_hdr = True
    hdr_lines = ''
    geom = []
    with open(fname) as f:
        for line in f:
            line0 = line
            line = line.strip()
            if not line:
                if save_hdr:
                    hdr_lines += line0
                continue # skip empty lines
            
            hdr = re.match('\\[([\\w\\s\\d]+)\\]', line) # '^' is not necessary for .match
            if hdr is not None:
                current_section = hdr.group(1)
                
                # Some special treatment for 1-line sections:                
                # Process the Cartesian/Pure keywords
                if current_section in ['5D', '5D7F']:
                    num_L_funcs['d'] = 5
                    num_L_funcs['f'] = 7
                elif current_section == '5D10F':
                    num_L_funcs['d'] = 5
                elif current_section == '7F':
                    num_L_funcs['f'] = 7
                elif current_section == '9G':
                    num_L_funcs['g'] = 9
            else:
                # process geometry section
                if current_section == 'Atoms':
                    
                    atm = re.match('([\\w]+).+\\s+([\\d+-.ED]+)\\s+([\\d+-.ED]+)\\s+([\\d+-.ED]+)$', line)
                    if atm is None:
                        print 'ERROR in line: ', line
                    else:
                        geom.append( [atm.group(1)] + [float(atm.group(i)) for i in [2,3,4] ] )
                    
                # process the basis set section
                elif current_section == 'GTO':
                    line = line.lower()
                    atNum = re.match('([0-9]+)\\s*[0-9]+$', line)
                    if atNum is not None:
                        current_atom_Id = int(atNum.group(1))
                    bf_type = re.match('([spdfghi]+)\\s+[\\d]+\\s*[\\d.]*', line)
                    if bf_type is not None:
                        # add this basis function L value to the current atom's list
                        if current_atom_Id not in atomicBFs:
                            atomicBFs[current_atom_Id] = []
                        atomicBFs[current_atom_Id].append( bf_type.group(1) )
                
                
                # process MO coefficient section
                elif current_section == 'MO':
                    save_hdr = False
                    mo_property = re.match('([\\w]+)\\s*=\\s*([\\d.\\w+-]+)', line)
                    
                    # is this line a kind of ' Sym=     1a' ?
                    if mo_property is not None:
                        prop_name = mo_property.group(1)
                        if prop_name not in ['Sym', 'Ene', 'Spin', 'Occup']:
                            print 'ERROR: unknown MO property: ', prop_name
                        else:
                            if current_MO is not None and (current_MO['coefs']):
                                #print 'NEW MO', line
                                # if we encountered MO property for a previously 'started'
                                # MO which has non-empty current_MO['coefs'] dict
                                
                                MOs.append( current_MO ) # Ok, a prev. MO has been completed
                                current_MO = None
                            
                            if current_MO is None:
                                current_MO = {'coefs': {}}
                                
                            current_MO[prop_name] = mo_property.group(2) # save the value of the property
                            
                    # is this line a kind of '  1       0.708282036084' ?
                    else:
                        cf_line = re.match('([\\d]+)\\s*([\\d.+-EDed]+)', line)
                        if cf_line is not None:
                            current_MO['coefs'][int(cf_line.group(1))] = float(cf_line.group(2))
                        else:
                            print 'ERROR: bad MO coefficient line: ', line
                        
                else:
                    print 'Skipping [%s]: ' % current_section, line
            if save_hdr:
                hdr_lines += line0
                    

    # get the total number of BFs
    print num_L_funcs
    tot_BFs = 0
    atom_BF_ids = {}
    for atom in sorted(atomicBFs.keys()):
        atom_BF_ids[atom] = []        
        bf_at_atom = 0
        for bf in atomicBFs[atom]:
            bf_at_atom += num_L_funcs[bf]
            for i in range(num_L_funcs[bf]):
                atom_BF_ids[atom] . append(tot_BFs + i)                
            tot_BFs += num_L_funcs[bf]
        
        print atom, 'has', bf_at_atom, 'basis functions'
        
    print 'All atoms have', tot_BFs, 'basis functions'
    print 'Distribution of BFs over atoms: ',atom_BF_ids
    
    # save the last MO
    if current_MO is not None:
        MOs.append( current_MO ) # Ok, a prev. MO has been completed

    # sort basis function numbers in MOs and convert coefs. to lists
    MO_cfs = { 'Alpha': [], 'Beta': [] }
    MO_occ = { 'Alpha': [], 'Beta': [] }
    MO_ene = { 'Alpha': [], 'Beta': [] }
    # TODO: order MOs by occupancy
    for i, m in enumerate(MOs):
        m['coefs'] = [cf[1] for cf in sorted(m['coefs'].items(), key = lambda x: x[0])]        
        cf_len = len(m['coefs'])
        if cf_len != tot_BFs:
            print ('ERROR: number of basis functions (%d) is different '+\
                   'from the number of MO %d coefs (%d)') % (tot_BFs, i, cf_len)
        MO_cfs[ m['Spin'] ].append(m['coefs'])
        MO_occ[ m['Spin'] ].append(float(m['Occup']))
        MO_ene[ m['Spin'] ].append(float(m['Ene']))
        #MO_n.append()

    Va = np.array(MO_cfs['Alpha']).T
    Vb = np.array(MO_cfs['Beta']).T
    na = np.array(MO_occ['Alpha'])
    nb = np.array(MO_occ['Beta'])
    Ea = np.array(MO_ene['Alpha'])
    Eb = np.array(MO_ene['Beta'])
    #IPython.embed()
    return Va, Vb, na, nb, Ea, Eb, atom_BF_ids, hdr_lines, geom
#-----------------



def save2molden(fname, v_ao, pre_MO_header, ene = None, occ = None):
        
    with open(fname, 'w') as f:
        f.write(pre_MO_header)
        for i, cf in enumerate(v_ao.T):
            f.write((' Sym=     1a\n'+
             ' Ene= %.12f\n'+
             ' Spin= Alpha\n'+
             ' Occup= %.5f\n') % (i if ene is None else ene[i], 0 if occ is None else occ[i]))
            for j,c in enumerate(cf):
                f.write('%3d   %18.12f\n' % (j+1, c))
