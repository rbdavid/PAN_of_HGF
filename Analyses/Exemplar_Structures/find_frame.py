#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

import numpy as np
import MDAnalysis
import sys

pdb_file = sys.argv[1]
rmsd_file = sys.argv[2]

rmsd_data = np.loadtxt(rmsd_file)

frame = np.argmin(rmsd_data[:,4])

print frame

traj_num = int(frame/5000)	# 5000 frames per trajectory...

frame_num = frame%5000

print traj_num, frame_num

u = MDAnalysis.Universe(pdb_file, '../Trajectories/production.%s/production.%s.dcd' %(traj_num, traj_num))

#system = u.select_atoms('protein or nucleic or resname A5 or resname A3 or resname U5 or resname atp or resname adp or resname PHX or resname MG')
system = u.select_atoms('all')
backbone = u.select_atoms('backbone')

u.trajectory[frame_num]
system.translate(-backbone.center_of_mass())

system.write('frame_%s.pdb' %(frame))

