
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

import sys
import numpy as np
from numpy.linalg import *
import MDAnalysis
from MDAnalysis.analysis.align import *
from distance_functions import *

# ----------------------------------------
# VARIABLE DECLARATION

pdb_file = sys.argv[1]
traj_loc = sys.argv[2]
model = int(sys.argv[3])
start = int(sys.argv[4])
end = int(sys.argv[5])
out_file = sys.argv[6]

alignment = 'protein and name CA'
important = 'protein and not name H*'

zeros = np.zeros
dot_prod = np.dot
sqrt = np.sqrt
flush = sys.stdout.flush

thresh = 1E-5
maxIter = 100

# ----------------------------------------
# FUNCTIONS:

def ffprint(string):
    print('%s' %(string))
    flush()

# ----------------------------------------
# MAIN PROGRAM:

# LOAD IN THE UNIVERSE TO BE ANALYZED AND INITIATE THE NECESSARY ATOM SELECTIONS
u = MDAnalysis.Universe(pdb_file)
u_all = u.select_atoms('all')
u_align = u.select_atoms(alignment)
u_important = u.select_atoms('protein and not name H*')

# INITIATE THE PDB FILE TO BE USED TO SAVE THE AVERAGE STRUCTURE
avg_pdb = MDAnalysis.Universe(pdb_file)
avg_important = avg_pdb.select_atoms('protein and not name H*')

# GRABBING IMPORTANT NUMBERS FROM THE UNIVERSE
u_important_atoms = len(u_important.atoms)
u_align_atoms = len(u_align.atoms)

# DETERMINING THE NUMBER OF STEPS TO BE AVERAGED OVER
temp = start
nSteps = 0
while temp <= end:
    u.load_new('%shgf_wt_model_%s.prod_%02d.trunc.dcd' %(traj_loc,model,start))
    nSteps += len(u.trajectory)
    temp += 1

ffprint('Number of steps to be averaged over : %d' %(nSteps))

# ARRAY DECLARATION
all_coord = zeros((nSteps,u_important_atoms,3),dtype=np.float32)
avgCoord = zeros((u_important_atoms,3),dtype=np.float32)
all_align = zeros((nSteps,u_align_atoms,3),dtype=np.float32)
avgAlign = zeros((u_align_atoms,3),dtype=np.float32)

# Trajectory Analysis: 
ffprint('Beginning trajectory analysis')
temp = 0 
while start <= end:
    ffprint('Loading trajectory %s (%shgf_wt_model_%s.prod_%02d.trunc.dcd)' %(start,traj_loc,model,start))
    u.load_new('%shgf_wt_model_%s.prod_%02d.trunc.dcd' %(traj_loc,model,start))
    
    for ts in u.trajectory:
        u_all.translate(-u_align.center_of_mass())
        
        avgCoord += u_important.positions
        avgAlign += u_align.positions
        all_coord[temp] = u_important.positions
        all_align[temp] = u_align.positions
        temp += 1
    start += 1

ffprint(nSteps)
if temp != nSteps:
    ffprint('Failed to analyze all timesteps; fucked shit up')

avgCoord /= float(nSteps)
avgAlign /= float(nSteps)
ffprint(len(all_coord))
ffprint(len(all_coord[0]))
ffprint('Finished with the trajectory analysis')

# Calculating and Aligning to the average positions
iteration = 0
residual = thresh + 10.0 					# arbitrary assignment greater than thresh
ffprint('Beginning iterative process of calculating average positions and aligning to the average')
while residual > thresh and iteration < maxIter:
    tempAvgCoord = zeros((u_important_atoms,3),dtype=np.float32)		# zeroing out the tempAvgCoord array every iteration
    tempAvgAlign = zeros((u_align_atoms,3),dtype=np.float32)
    for i in range(nSteps):
        R, d = rotation_matrix(all_align[i,:,:],avgAlign)
        all_align[i,:,:] = dot_prod(all_align[i,:,:],R.T)
        all_coord[i,:,:] = dot_prod(all_coord[i,:,:],R.T)
        tempAvgAlign += all_align[i,:,:]
        tempAvgCoord += all_coord[i,:,:]			# recalculate the average coordinates to optimize the average position
    tempAvgCoord /= float(nSteps)				# finishing the average
    tempAvgAlign /= float(nSteps)
    residual = RMSD(avgAlign, tempAvgAlign, u_align_atoms)			# calculating the rmsd between avg and tempAvg to quantify our iterative optimization of the average positions	
    rmsd_all = RMSD(avgCoord, tempAvgCoord, u_important_atoms)
    iteration += 1
    avgCoord = tempAvgCoord
    avgAlign = tempAvgAlign
    ffprint('Steps: %d, RMSD btw alignment atoms: %e, RMSD btw all atoms: %e' %(iteration,residual,rmsd_all))
ffprint('Average structure has converged')				# Now have the iteratively aligned avgCoord array, as well as the iteratively aligned (COG-corrected and rotated) allCoord array

# Print out pdb of average structure
ffprint('Writing a pdb of the average structure.')
avg_important.positions = avgCoord
avg_important.write('%03d.%03d.avg_structure.pdb' %(int(sys.argv[4]),end))
ffprint('Finished writing pdb of the average structure')

# APPENDING INFORMATION TO A SUMMARY FILE
with open('%s' %(out_file),'a') as out:
    out.write('%d   %d   %d\n' %(int(sys.argv[4]), end, nSteps))

