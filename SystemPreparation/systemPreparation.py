#!/usr/bin/env python
# coding: utf-8

import BioSimSpace as BSS

# Using files in AMBER or GROMACS format, BioSimSpcae can handel both.

system = BSS.IO.readMolecules(['your/path/to/the/files', 'your/path/to/files'])

## Information about the system 

# make sure you read the correct file type

system.fileFormat()

# check the number of molecules in the system

system.nMolecules()


# how many of these molecules are water molecules

system.nWaterMolecules()

# the total charge of the system

system.charge()

# see info about all the molecules in the system eg nAtoms, nResidues, charge 

print(f'{"ID":5}{"Number of Atoms":20}{"Number of Residues"}{"The charge of the molecules":>30}')

for molecule in range(system.nMolecules()):
    print(f'{molecule+1:<5}{system[molecule].nAtoms():< 25}{system[molecule].nResidues():<20}{system[molecule].charge()}')


## Applying a force field to the system

# When applying a force filed, you cant pass the whole system as it will give u an error 
# Because the force field used for protien not the same as for the molcules
# here we will not apply a force field as the it is already applied 

#gaff_system = BSS.Parameters.gaff(system).getMolecule()


## Solvation of the system and neutralization of the charge 
##### NB: Solvation can occur to the whole system at once not as applying force field


# #### First calculate the box size in which the system will be sovated

# Get the minimium and maximum coordinates of the bounding box that
# minimally encloses the protein.
box_min, box_max = system.getAxisAlignedBoundingBox()

# Work out the box size from the difference in the coordinates.
box_size = [y - x for x, y in zip(box_min, box_max)]

# How much to pad each side of the protein? (Nonbonded cutoff = 10 A)
padding = 15 * BSS.Units.Length.angstrom

# Work out an appropriate box. This will used in each dimension to ensure
# that the cutoff constraints are satisfied if the molecule rotates.
box_length = max(box_size) + 2*padding


# #### After calculating the size of the box, we place the system, water, counter ions inside it.

# tip3 water model and ~150 nM of Nacl is used to neutralize and solvate the system 

solavted_system = BSS.Solvent.tip3p(system, box = 3 * [box_length], ion_conc=0.15)


# checking the charge of the system after adding counterions to nutralize the charge of the system
# Ideally it will be zero, but very small charge that is ~ zero is allowed

solavted_system.charge()

## Viewing of the solvated system

view = BSS.Notebook.View(solavted_system)
view.system()


# Before running the Simulation process, we have to do Energy Minimization and Equiliburation for the system

## Firstly, Energy Minimization

#### The goal of Energy minimization is to remove any steric clashes or unusual geometry which would artificially raise the energy of the system, we must relax the structure by running an energy minimization (EM) algorithm.

#Energy Minimisation for the whole system

minimisation_protocol = BSS.Protocol.Minimisation(steps=1000)

# Choosing the engine to do the minimisation process

minimisation_process = BSS.Process.Amber(solavted_system, minimisation_protocol)

# Checking the source of the engine used to carry out the process

minimisation_process.exe()


# Start the minimization process

minimisation_process.start()

# Checking if it still running or not

minimisation_process.isRunning()


# Checking the time taken to finish the process

minimisation_process.runTime()

# Getting the system after finishing the minimization process

minimisation_system = minimisation_process.getSystem(block=True)

# checking whether the process done successfully or not

minimisation_process.isError()


# Plotting energy vs time

plot = BSS.Notebook.plot(minimisation_process.getTime(time_series=True), minimisation_process.getTotalEnergy(time_series=True))


# Saving the minimized system

BSS.IO.saveMolecules('path/to/the/folder/to/save/the/minimized/system', minimisation_system, ["prm7", "rst7","pdb"])

# ## Secondaly, Equilibiration

# #### At this point equilibration of the solvent around the solute (i.e. the protein) is necessary. This is performed in two stages: equilibration under an NVT ensemble, followed by an NPT ensemble. Use of the NVT ensemble entails maintaining constant number of particles, volume and temperature, while the NPT ensemble maintains constant number of particles, pressure and temperature. (The NVT ensemble is also known as the isothermal-isochoric ensemble, while the NPT ensemble is also known as the isothermal-isobaric ensemble).

# ### Firstly, NVT Equilibiration

# #### NVT equilibration for 5 ps while restraining all non-solvent atoms


equil_protocol_1 = BSS.Protocol.Equilibration(

                                runtime=5*BSS.Units.Time.picosecond, 
                                temperature_start=0*BSS.Units.Temperature.kelvin, 
                                temperature_end=300*BSS.Units.Temperature.kelvin,
                                restraint="all"
                                
)

# pass the system and the protocol to the engine

equil_process_1 = BSS.Process.Amber(minimisation_system, equil_protocol_1, exe='/your/path/to/pmemd.cuda')


# Checking which engine will carry out the process eg sandr or pmemd

equil_process_1.exe()

equil_process_1.start()

equil_process_1.isRunning()

equil_system_1 = equil_process_1.getSystem(block=True)

equil_process_1.isError()

type(equil_system_1)

# Generate a plot of time vs temperature.
plot1 = BSS.Notebook.plot(equil_process_1.getTime(time_series=True), equil_process_1.getTemperature(time_series=True))

# Generate a plot of time vs energy.
plot2 = BSS.Notebook.plot(equil_process_1.getTime(time_series=True), equil_process_1.getTotalEnergy(time_series=True))


# #### Sander NVT equilibration for 50 ps without restrains


equil_protocol_2 = BSS.Protocol.Equilibration(
                                runtime=50*BSS.Units.Time.picosecond, 
                                temperature=300*BSS.Units.Temperature.kelvin,
                                restraint="backbone"
                                )

equil_process_2 = BSS.Process.Amber(equil_system_1, equil_protocol_2, exe='/your/path/to/pmemd.cuda')

equil_process_2.exe()

equil_process_2.start()

equil_process_2.isRunning()

equil_system_2 = equil_process_2.getSystem(block=True)

equil_process_2.isError()

# Generate a plot of time vs temperature.
plot1 = BSS.Notebook.plot(equil_process_2.getTime(time_series=True), equil_process_2.getTemperature(time_series=True))

# Generate a plot of time vs energy.
plot2 = BSS.Notebook.plot(equil_process_2.getTime(time_series=True), equil_process_2.getTotalEnergy(time_series=True))


# #### NVT equilibiration for 50 ps without restraints



equil_protocol_3 = BSS.Protocol.Equilibration(

                                runtime=50*BSS.Units.Time.picosecond, 
                                temperature_end=300*BSS.Units.Temperature.kelvin
)


# #### Sander NPT equilibration for 200 ps while restraining non-solvent heavy atoms

equil_process_3 = BSS.Process.Amber(equil_system_2, equil_protocol_3, exe='/your/path/to/pmemd.cuda')

equil_process_3.exe()

equil_process_3.start()

equil_process_3.isRunning()

equil_system_3 = equil_process_3.getSystem(block=True)

equil_process_3.isError()

# Generate a plot of time vs temperature.
plot1 = BSS.Notebook.plot(equil_process_3.getTime(time_series=True), equil_process_3.getTemperature(time_series=True))

# Generate a plot of time vs energy.
plot2 = BSS.Notebook.plot(equil_process_3.getTime(time_series=True), equil_process_3.getTotalEnergy(time_series=True))


#### pmemd.cuda NPT equilibration for 200 ps while restraining non-solvent heavy atoms
equil_protocol_4 = BSS.Protocol.Equilibration(
                                runtime=200*BSS.Units.Time.picosecond, 
                                pressure=1*BSS.Units.Pressure.atm,
                                temperature_end=300*BSS.Units.Temperature.kelvin,
                                restraint="heavy"
                                )

equil_process_4 = BSS.Process.Amber(equil_system_3, equil_protocol_4, exe='/your/path/to/pmemd.cuda')

equil_process_4.exe()

equil_process_4.start()

equil_process_4.isRunning()

equil_system_4 = equil_process_4.getSystem(block=True)

equil_process_4.isError()

# Generate a plot of time vs temperature.
plot1 = BSS.Notebook.plot(equil_process_4.getTime(time_series=True), equil_process_4.getTemperature(time_series=True))

# Generate a plot of time vs energy.
plot2 = BSS.Notebook.plot(equil_process_4.getTime(time_series=True), equil_process_4.getTotalEnergy(time_series=True))

equil_protocol_5 =  BSS.Protocol.Equilibration(
                                runtime=200*BSS.Units.Time.picosecond, 
                                pressure=1*BSS.Units.Pressure.atm,
                                temperature=300*BSS.Units.Temperature.kelvin
                                )


equil_process_5 = BSS.Process.Amber(equil_system_4, equil_protocol_5, exe='/your/path/to/pmemd.cuda')

equil_process_5.exe()

equil_process_5.start()

equil_process_5.isRunning()

equil_system_5 = equil_process_5.getSystem(block=True)

equil_process_5.isError()

# Generate a plot of time vs temperature.
plot1 = BSS.Notebook.plot(equil_process_5.getTime(time_series=True), equil_process_5.getTemperature(time_series=True))

# Generate a plot of time vs energy.
plot2 = BSS.Notebook.plot(equil_process_5.getTime(time_series=True), equil_process_5.getTotalEnergy(time_series=True))

# Saving the NVT system

BSS.IO.saveMolecules('path/to/the/folder/to/save/final_system', equil_system_5, ["prm7", "rst7","pdb"])
