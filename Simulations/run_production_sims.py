
import sys
import time
import simtk
import parmed

prmtop = sys.argv[1]
chk_file = sys.argv[2]
output_dir = sys.argv[3]
prod_id = int(sys.argv[4])

time_step   = 2*simtk.unit.femtosecond	# simulation timestep
ts_float    = 2*10**-6			# femtoseconds to nanosecond
temperature = 300*simtk.unit.kelvin     # simulation temperature
friction    = 1/simtk.unit.picosecond   # collision rate
pressure    = 1.01325*simtk.unit.bar	# simulation pressure 
mcbarint    = 100			# number of steps between volume change attempts
num_steps   = 5000000                   # number of integration steps to run; 10 ns per run
trj_freq    = 1000                      # number of steps per written trajectory frame
data_freq   = 1000                      # number of steps per written simulation statistics

# SETTING up the simulation system
print('Loading the AMBER files')
parmed_system = parmed.load_file(prmtop)
print('Creating OpenMM system')
system = parmed_system.createSystem(nonbondedMethod=simtk.openmm.app.PME, nonbondedCutoff=1.2*simtk.unit.nanometer,constraints=simtk.openmm.app.HBonds)

# SETTING up the Langevin dynamics thermostat.
integrator = simtk.openmm.LangevinIntegrator(temperature, friction, time_step)

# SETTING up the Monte Carlo barostat.
if sum(parmed_system.get_box()[0][:3]) != 0.:
    barostat = simtk.openmm.MonteCarloBarostat(pressure,temperature,mcbarint)
    system.addForce(barostat)
else:
    print('Pressure cannot be implemented because no PBC are defined')

# SETTING the simulation platform .
platform = simtk.openmm.Platform.getPlatformByName('CUDA')
prop     = dict(CudaPrecision='mixed') # Use mixed single/double precision

# SETTING up an OpenMM simulation.
simulation = simtk.openmm.app.Simulation(parmed_system.topology, system, integrator, platform, prop)

# LOADING a checkpoint file
with open(chk_file,'rb') as f:
    simulation.context.loadCheckpoint(f.read())

# SETTING up output files.
file_naming = output_dir + prmtop.split('/')[-1].split('.')[0] 
simulation.reporters.append(simtk.openmm.app.dcdreporter.DCDReporter(file_naming+'.prod_%02d.dcd'%(prod_id),trj_freq))
simulation.reporters.append(simtk.openmm.app.statedatareporter.StateDataReporter(file_naming+'.prod_%02d.out'%(prod_id),data_freq,step=True,potentialEnergy=True,kineticEnergy=True,temperature=True,volume=True,density=True,speed=True))

# RUNNING the simulation
print("Starting simulation")
start = time.process_time()
simulation.step(num_steps)
end = time.process_time()
print("Elapsed time %.2f seconds\nAverage speed: %.3f ns day^{-1}" % (end-start,(num_steps*ts_float*86400)/(end-start)))
print("Done!")

simulation.saveCheckpoint(file_naming+'.prod_%02d.chkpt'%(prod_id))
simulation.saveState(file_naming+'.prod_%02d.stt'%(prod_id))

