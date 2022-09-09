from numpy import number
from fileTransfer import *
from handelMORPHs import *
from numberOfRuns import *
from editCFG import *
from boreshRestraints import *
import shutil


# take arguments by input from the command line
print('Hello to run00X generator')
prm7 = input('Enter the path to the prm7 file: ')
rst7 = input('Enter the path to the rst7 file: ')
state = input('Enter the state/mode of the calculations (free, bound or all): ')
number = int(input('Enter the number of runs: '))

if state == 'bound' or state == 'all':
    question1 = input('Do you have a boresch restraints file? (y/n): ')
    if question1.startswith('y'):
        answer1 = input('Enter the path to the boresch restraints txt file: ')
        # open the txt file and read the first line
        with open(answer1, 'r') as f:
            boresch_restrains = f.readline()
            editCFG(boresch_restrains, state)
            print('''
************************************************************************
************************************************************************
********************* Run00X generator is running***********************
************************************************************************
************************************************************************
''')

    else:
        shutil.copy(prm7, './SYSTEM.top')
        nc_file = input('Enter the path to the nc file to calculate the boresch restrains: ')

        # calculate the boresch restraints from the nc file
        getboresch_restraints(nc_file)

        # the output of this function is a txt file, so we need to open it and read the first line
        with open('the_boresch_parm_used.txt', 'r') as f:
            boresch_restrains = f.readline()
            editCFG(boresch_restrains, state)
            print('''
************************************************************************
************************************************************************
********************* Run00X generator is running***********************
************************************************************************
************************************************************************
''')



fileTransfer([prm7, rst7], state)
handelMORPHs(prm7, rst7, state)
numberOfRuns(number, state)
shutil.copy('run_all.sh', './free')
shutil.copy('run_all.sh', './bound')
