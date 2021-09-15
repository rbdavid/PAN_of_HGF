
# ----------------------------------------
# USAGE:
# python3 rmsd.ref.py config_file

# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import MDAnalysis
from MDAnalysis.analysis.align import *
from distance_functions import *
from selection_list import *

flush = sys.stdout.flush

config_file = sys.argv[1]

# ----------------------------------------
# RESIDUE LISTS:

nucleic = ['A5','A3','A','G5','G3','G','C5','C3','C','T5','T3','T','U5','U3','U']
triphosphate = ['atp','adp','PHX']
other = ['MG']

# ----------------------------------------
# HOMEMADE ATOM SELECTION STRINGS:

sugar = "name C5' C4' O4' C1' C3' C2' O2' " + " C5* C4* O4* C1* C3* O3* C2* O2* "		# NO HYDROGENS; DOES NOT INCLUDE THE O5' atom (which I will include in the phosphate atom selection string...; the atoms with * are found in triphosphates;
sugar_5= sugar + " O5'"		# NO HYDROGENS
sugar_3= sugar + " O3' "	# NO HYDROGENS
base = 'name N9 C8 N7 C5 C6 N6 N1 C2 N3 C4 O6 N4 C2 O2 O4'	# NO HYDROGENS; selection string that will select all appropriate atoms for any of the nucleic residues...
a_phos = 'name O5* O2A O1A PA O3A'
b_phos = 'name PB O1B O2B O3B'
g_phos = 'name PG O1G O2G O3G'
inorg_phos = 'name P O1 O2 O3 O4'	# NO HYDROGENS

# ----------------------------------------
# FUNCTIONS: (NOT INCLUDING MAKE_SELECTIONS FUNCTION WHICH IS FOUND AT THE BOTTOM OF THIS SCRIPT)

def ffprint(string):
    print('%s'%(string))
    flush()

necessary_parameters = ['ref_pdb','pdb','traj_loc','start','end','Wrapped','outputname','selection_file']
all_parameters = ['ref_pdb','pdb','traj_loc','start','end','Wrapped','outputname','selection_file','alignment','wrapping_selection','substrates','homemade_selections','write_summary','summary_filename','selection_output']
def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
    for i in range(len(necessary_parameters)):
        parameters[necessary_parameters[i]] = ''
        
    # SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
    parameters['alignment'] = 'protein'
    parameters['wrapping_selection'] = 'not (resname WAT or resname Na+ or resname Cl- or protein)'
    parameters['substrates'] = None
    parameters['homemade_selections'] = None
    parameters['write_summary'] = False 
    parameters['summary_filename'] = 'rmsd.summary'
    parameters['selection_output'] = 'selections.txt'
    
    # GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
    with open(config_file) as f:
        exec(compile(f.read(),config_file,'exec'),parameters)
    
    for key, value in list(parameters.items()):
        if value == '':
            print('%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key))
            sys.exit()

def summary():
    with open('%s' %(parameters['summary_filename']),'w') as f:
        f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
        f.write('To recreate this analysis, run this line:\n')
        for i in range(len(sys.argv)):
            f.write('%s ' %(sys.argv[i]))
        f.write('\n\n')
        f.write('Parameters used:\n')
        for i in all_parameters:
            f.write('%s = %s \n' %(i,parameters[i]))
        f.write('\n\n')
        f.write('output is written to:\n')
        f.write('	%s\n' %(parameters['outputname']))
        f.write('\nNumber of steps analyzed: %d\n' %(nSteps))
        f.write('\nAtom selections analyzed have been written out to %s\n' %(parameters['selection_output']))

def make_selections(analysis_universe,ref_universe,resname,resid,output_file,selection_list,nAtoms,ref_pos):
    """A function that takes in a residue name and creates a non-standard MDAnalysis atom selection
    
    Usage: new_sel = make_selection(........)
    
    Arguments:
    	analysis_universe: MDAnalysis Universe object to be used as the analysis universe.
    	reference_universe: MDAnalysis Universe object to be used as the reference universe.
    	resname: string of the residue name;
    	resid: int of the residue ID number;
    	output_file: file object that is to be written to;
    """
    
    global count
    
    # ----------------------------------------
    # CREATING THE NUCLEIC SELECTIONS
    if resname in nucleic:
        # CREATING THE SLECTION FOR THE BASE OF NUCLEIC RESIDUES
        sel_string = 'resname %s and resid %d and %s' %(resname,resid,base)
        u_temp = analysis_universe.select_atoms(sel_string)
        selection_list.append(u_temp)
        nAtoms.append(u_temp.n_atoms)
        ref_temp = ref_universe.select_atoms(sel_string)
        ref_pos.append(ref_temp.positions)
        if u_temp.n_atoms != ref_temp.n_atoms:
            ffprint('Number of atoms do not match for selection %02d, %s, %s' %(count,resname,sel_string))
            sys.exit()
        output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
        count +=1
        
        # CREATING THE SLECTION FOR THE SUGAR OF NUCLEIC RESIDUES
        if resname in ['A5','G5','C5','T5','U5']:
            sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar_5)
            u_temp = analysis_universe.select_atoms(sel_string)
            selection_list.append(u_temp)
            nAtoms.append(u_temp.n_atoms)
            ref_temp = ref_universe.select_atoms(sel_string)
            ref_pos.append(ref_temp.positions)
            if u_temp.n_atoms != ref_temp.n_atoms:
                ffprint('Number of atoms do not match for selection %02d, %s, %s' %(count,resname,sel_string))
                sys.exit()
            output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
            count +=1
            return
            
        elif resname in ['A3','U3','C3','G3','T3']:
            sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar_3)
            u_temp = analysis_universe.select_atoms(sel_string)
            selection_list.append(u_temp)
            nAtoms.append(u_temp.n_atoms)
            ref_temp = ref_universe.select_atoms(sel_string)
            ref_pos.append(ref_temp.positions)
            if u_temp.n_atoms != ref_temp.n_atoms:
                ffprint('Number of atoms do not match for selection %02d, %s, %s' %(count,resname,sel_string))
                sys.exit()
            output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
            count +=1
            
        else:
            sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar)
            u_temp = analysis_universe.select_atoms(sel_string)
            selection_list.append(u_temp)
            nAtoms.append(u_temp.n_atoms)
            ref_temp = ref_universe.select_atoms(sel_string)
            ref_pos.append(ref_temp.positions)
            if u_temp.n_atoms != ref_temp.n_atoms:
                ffprint('Number of atoms do not match for selection %02d, %s, %s' %(count,resname,sel_string))
                sys.exit()
            output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
            count +=1
            
        # CREATING THE SLECTION FOR THE PHOSPHATE OF NUCLEIC RESIDUES
        sel_string = "(resname %s and resid %s and name P OP1 OP2 O5') or (resid %s and name O3')" %(resname,resid,analysis_universe.residues[resid-1].resid) 
        u_temp = analysis_universe.select_atoms(sel_string)
        selection_list.append(u_temp)
        nAtoms.append(u_temp.n_atoms)
        ref_temp = ref_universe.select_atoms(sel_string)
        ref_pos.append(ref_temp.positions)
        output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
        count += 1
        return
        
    # ----------------------------------------
    # CREATING THE TRIPHOSPHATE ATOM SELECTIONS
    elif resname in triphosphate:
        if resname in ['atp','adp']:
            sel_string = 'resname %s and resid %d and %s' %(resname,resid,base)
            u_temp = analysis_universe.select_atoms(sel_string)
            selection_list.append(u_temp)
            nAtoms.append(u_temp.n_atoms)
            ref_temp = ref_universe.select_atoms(sel_string)
            ref_pos.append(ref_temp.positions)
            if u_temp.n_atoms != ref_temp.n_atoms:
                ffprint('Number of atoms do not match for selection %02d, %s, %s' %(count,resname,sel_string))
                sys.exit()
            output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
            count +=1
            
            sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar)
            u_temp = analysis_universe.select_atoms(sel_string)
            selection_list.append(u_temp)
            nAtoms.append(u_temp.n_atoms)
            ref_temp = ref_universe.select_atoms(sel_string)
            ref_pos.append(ref_temp.positions)
            if u_temp.n_atoms != ref_temp.n_atoms:
                ffprint('Number of atoms do not match for selection %02d, %s, %s' %(count,resname,sel_string))
                sys.exit()
            output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
            count +=1
            
        if resname == 'atp':
            sel_string = 'resname %s and resid %d and %s' %(resname,resid,a_phos)
            u_temp = analysis_universe.select_atoms(sel_string)
            selection_list.append(u_temp)
            nAtoms.append(u_temp.n_atoms)
            ref_temp = ref_universe.select_atoms(sel_string)
            ref_pos.append(ref_temp.positions)
            output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
            count +=1
            
            sel_string = 'resname %s and resid %d and %s' %(resname,resid,b_phos)
            u_temp = analysis_universe.select_atoms(sel_string)
            selection_list.append(u_temp)
            nAtoms.append(u_temp.n_atoms)
            ref_temp = ref_universe.select_atoms(sel_string)
            ref_pos.append(ref_temp.positions)
            output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
            count +=1
            
            sel_string = 'resname %s and resid %d and %s' %(resname,resid,g_phos)
            u_temp = analysis_universe.select_atoms(sel_string)
            selection_list.append(u_temp)
            nAtoms.append(u_temp.n_atoms)
            ref_temp = ref_universe.select_atoms(sel_string)
            ref_pos.append(ref_temp.positions)
            output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
            count +=1
            return
            
        elif resname == 'adp':
            sel_string = 'resname %s and resid %d and %s' %(resname,resid,a_phos)
            u_temp = analysis_universe.select_atoms(sel_string)
            selection_list.append(u_temp)
            nAtoms.append(u_temp.n_atoms)
            ref_temp = ref_universe.select_atoms(sel_string)
            ref_pos.append(ref_temp.positions)
            output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
            count +=1
            
            sel_string = 'resname %s and resid %d and %s' %(resname,resid,b_phos)
            u_temp = analysis_universe.select_atoms(sel_string)
            selection_list.append(u_temp)
            nAtoms.append(u_temp.n_atoms)
            ref_temp = ref_universe.select_atoms(sel_string)
            ref_pos.append(ref_temp.positions)
            output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
            count +=1
            return
            
        elif resname == 'PHX':
            sel_string = 'resname %s and resid %d and %s' %(resname,resid,inorg_phos)
            u_temp = analysis_universe.select_atoms(sel_string)
            selection_list.append(u_temp)
            nAtoms.append(u_temp.n_atoms)
            ref_temp = ref_universe.select_atoms(sel_string)
            ref_pos.append(ref_temp.positions)
            output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
            count +=1
            return
            
    # ----------------------------------------
    # CREATING ANY REMAINING ATOM SELECTIONS...
    elif resname in other:
        sel_string = 'resname %s and resid %d' %(resname,resid)
        u_temp = analysis_universe.select_atoms(sel_string)
        selection_list.append(u_temp)
        nAtoms.append(u_temp.n_atoms)
        ref_temp = ref_universe.select_atoms(sel_string)
        ref_pos.append(ref_temp.positions)
        output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
        count +=1
        return
    
# ----------------------------------------
# MAIN:
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

# ----------------------------------------
# LOAD IN THE REFERENCE AND ANALYSIS UNIVERSES
ref = MDAnalysis.Universe(parameters['ref_pdb'])
ref_align = ref.select_atoms(parameters['alignment'])
ref_all = ref.select_atoms('all')
ref_all.translate(-ref_align.center_of_mass())
pos0 = ref_align.positions

u = MDAnalysis.Universe(parameters['pdb'])
u_align = u.select_atoms(parameters['alignment'])
u_all = u.select_atoms('all')
if not parameters['Wrapped']:
    rest = u.select_atoms(parameters['wrapping_selection'])
    rest_nRes = rest.n_residues

if u_align.n_atoms != ref_align.n_atoms:
    ffprint('Alignment atom selections do not have the same number of atoms.')
    sys.exit()

# ----------------------------------------
# CREATION OF ATOM SELECTIONS FOR RMSD ANALYSIS

nSel = len(sel)

selection_list = []
nAtoms = []
ref_pos = []
with open('%s' %(parameters['selection_output']),'w') as f:
    for i in range(nSel):
        u_temp = u.select_atoms(sel[i][1])
        selection_list.append(u_temp)
        nAtoms.append(u_temp.n_atoms)
        ref_temp = ref.select_atoms(sel[i][1])
        ref_pos.append(ref_temp.positions)
        if u_temp.n_atoms != ref_temp.n_atoms:
            ffprint('Number of atoms do not match for selection %02d %s' %(i,sel[i][0]))
            sys.exit()
            
        f.write("%02d   %s   '%s'\n" %(i,sel[i][0],sel[i][1]))
        
    count = nSel
    
    if parameters['substrates'] != None:
        u_subs = u.select_atoms(parameters['substrates'])
        for i in range(u_subs.n_residues):
            temp_resname = u_subs.residues[i].resname
            temp_id = u_subs.residues[i].resid
            if temp_resname in parameters['homemade_selections']:
                make_selections(u,ref,temp_resname,temp_id,f,selection_list,nAtoms,ref_pos)
                
            else:
                u_temp = u.select_atoms('resname %s and resid %d' %(temp_resname,temp_id))
                selection_list.append(u_temp)
                nAtoms.append(u_temp.n_atoms)
                ref_temp = ref.select_atoms('resname %s and resid %d' %(temp_resname, temp_id))
                ref_pos.append(ref_temp.positions)
                if u_temp.n_atoms != ref_temp.n_atoms:
                    ffprint('Number of atoms do not match for selection %02d %s' %(i,sel[i][0]))
                    sys.exit()
                    
                f.write('%02d   %s   resid %d\n' %(count,temp_resname,temp_id))
                count += 1
                
nSteps = 0
start = int(parameters['start'])
end = int(parameters['end'])
with open(parameters['outputname'],'w') as f:
    ffprint('Beginning trajectory analysis')
    while start <= end:
        ffprint('Loading trajectory %s' %(start))
        u.load_new(parameters['traj_loc']%(start))
        nSteps += len(u.trajectory)
        # Loop through trajectory
        for ts in u.trajectory:
            # Align to reference (moves COM of alignment to origin)
            u_all.translate(-u_align.center_of_mass())
            
            # CALCULATIONS that are unnecessary if the trajectory is wrapped.
            if not parameters['Wrapped']:		# Test to see if the 'Wrapped' key is equal to False
                dims = u.dimensions[:3]		
                dims2 = dims/2.0
                
                for i in range(rest_nRes):
                    COM = rest.residues[i].center_of_mass()
                    t = wrapping(COM,dims,dims2)
                    rest.residues[i].atoms.translate(t)
                    
            # Calculate the rotational matrix to align u to the ref
            R, rmsd = rotation_matrix(u_align.positions, pos0)
            # Apply rotation matrix to atoms within u
            u_all.rotate(R)
            
            # loop through selections and compute RMSD
            for i in range(count):
                temp_coords = selection_list[i].positions
                
                rmsd = RMSD(temp_coords,ref_pos[i],nAtoms[i])
                f.write('%10.6f   ' %(rmsd))
                
            f.write('\n')
        start += 1
        
if parameters['write_summary']:
    summary()

