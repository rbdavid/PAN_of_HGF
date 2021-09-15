
import sys
import time
import simtk
import parmed

prmtop = sys.argv[1]
rst7   = sys.argv[2]
output_dir = sys.argv[3]

time_step   = 2*simtk.unit.femtosecond	# simulation timestep
ts_float    = 2*10**-6			# femtoseconds to nanosecond
temperature = 300*simtk.unit.kelvin     # simulation temperature
friction    = 1/simtk.unit.picosecond   # collision rate
pressure    = 1.01325*simtk.unit.bar	# simulation pressure 
mcbarint    = 100			# number of steps between volume change attempts
num_steps   = 10000                     # number of integration steps to run
trj_freq    = 1000                      # number of steps per written trajectory frame
data_freq   = 1000                      # number of steps per written simulation statistics

# SETTING up the simulation system
print('Loading the AMBER files')
parmed_system = parmed.load_file(prmtop,rst7)
print('Creating OpenMM system')
system = parmed_system.createSystem(nonbondedMethod=simtk.openmm.app.PME, nonbondedCutoff=1.2*simtk.unit.nanometer,constraints=simtk.openmm.app.HBonds)

# SETTING up the Langevin dynamics thermostat.
integrator = simtk.openmm.LangevinIntegrator(temperature, friction, time_step)

# SETTING the simulation platform .
platform = simtk.openmm.Platform.getPlatformByName('CUDA')
prop     = dict(CudaPrecision='mixed') # Use mixed single/double precision

# SETTING up an OpenMM simulation.
simulation = simtk.openmm.app.Simulation(parmed_system.topology, system, integrator, platform, prop)

# SETTING the initial positions.
simulation.context.setPositions(parmed_system.positions)

# SETTING the velocities from a Boltzmann distribution at a given temperature.
simulation.context.setVelocitiesToTemperature(temperature)

# SETTING up output files.
file_naming = output_dir + prmtop.split('/')[-1].split('.')[0] 
simulation.reporters.append(simtk.openmm.app.dcdreporter.DCDReporter(file_naming+'.nvt_equilib.dcd',trj_freq))
simulation.reporters.append(simtk.openmm.app.statedatareporter.StateDataReporter(file_naming+'.nvt_equilib.out',data_freq,step=True,potentialEnergy=True,kineticEnergy=True,temperature=True,volume=True,density=True,speed=True))

# RUNNING the simulation
print("Starting simulation")
start = time.process_time()
simulation.step(num_steps)
end = time.process_time()
print("Elapsed time %.2f seconds\nAverage speed: %.3f ns day^{-1}" % (end-start,(num_steps*ts_float*86400)/(end-start)))
print("Done!")

simulation.saveCheckpoint(file_naming+'.nvt_equilib.chkpt')
simulation.saveState(file_naming+'.nvt_equilib.stt')

