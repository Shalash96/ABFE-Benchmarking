import BioSimSpace as BSS
import sys


def simulation (file1 , file2):

    # read the system which you want to do simulation for i.e the system produced from the system preparation step.
    system = BSS.IO.readMolecules([file1, file2])

    # simulatoin id done by using production protocol
    # specify the length of the si,ulation process
    protocol = BSS.Protocol.Production(runtime=50 * BSS.Units.Time.nanosecond)

    # choose the engine you want to carry on the simulation process eg AMBER-pmemd.cuda
    # work_dir is the folder where the products of the simulation (eg trajactory file -nc-) process are present after the finishing the simulation. 
    simProcess = BSS.Process.Amber(system, protocol, exe='/your/path/to/pmemd.cuda', work_dir='/path/to/the/folder/to/the/productOfSimulationProcess')
    simProcess.start()

    # check if the process is running or not
    # in case if not running, output files will be created containing the error 
    if simProcess.isRunning() == False:
        simProcess.getOutput('/path/to/the/folder/to/thefailedSimulationRunning')
        sys.exit()

    simSystem = simProcess.getSystem(block=True)

    # check if the process completed correctly or with errors
    # in case of error, output files will be created containing the error 

    if simProcess.isError() == True:
        simProcess.getOutput(f'/path/to/the/folder/to/the/failedSimulation')
        sys.exit()


    
    # save the system as pdb, rst7 and prm7 formats
    BSS.IO.saveMolecules('/path/to/the/folder/to/the/simulatedSystem',simSystem, ["prm7", "rst7","pdb"])


if __name__ == "__main__":
    simulation(sys.argv[1], sys.argv[2])
