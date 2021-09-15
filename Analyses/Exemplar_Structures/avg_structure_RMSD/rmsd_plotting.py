
# USAGE:
# python3 rmsd_plotting.py rmsd_data_file selection_list descriptor

# PREAMBLE:

from plotting_functions import *

data_file = sys.argv[1]  
selection_list = sys.argv[2]
system = sys.argv[3]

# ----------------------------------------
# MAIN:

datalist = np.loadtxt(data_file)
nSels = len(datalist[0])
nSteps = len(datalist)
print('Number of selections: %d, number of steps: %d' %(nSels,nSteps))

time = np.zeros(nSteps)
for i in range(nSteps):
    time[i] = i*0.002		# units of time in ns; each frame is separated by 0.002 ns 

selection_titles = []
with open(selection_list,'r') as f:
    for line in f:
        temp = line.split('   ')
        selection_titles.append('%02d.%s'%(int(temp[0]),temp[1]))

for i in range(nSels):
    scat_hist(time[:],datalist[:,i],'k','Time (ns)','RMSD',selection_titles[i],'%s' %(system),yunits='$\AA$',x_lim=(0,200),y_lim=(0,4.0))

