
sel = []
sel.append(['aligned_CAs','protein and name CA'])
sel.append(['full_protein','protein and not name H*'])
sel.append(['full_backbone','backbone'])
sel.append(['full_protein_no_sym_issues','protein and not name H* and not (resname ARG and name NH1 NH2) and not (resname ASP and name OD1 OD2) and not (resname GLU and name OE1 OE2) and not (resname LEU and name CD1 CD2) and not (resname PHE TYR and name CD1 CD2 CE1 CE2) and not (resname VAL and name CG1 CG2)'])

