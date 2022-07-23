from configparser import ConfigParser

# # read lines from a file
# with open('the_boresch_parm_used.txt', 'r') as f:
#     # just read the first line, as it only contains one line
#     boresch_restrains = f.readline()


def editCFG(boresch_restrains, state='all'):
    if state != 'free':
        parser = ConfigParser()
        parser.read('bound_template/discharge/input/sim.cfg')
        parser.set('HEAD', 'boresch restraints dictionary', boresch_restrains)
        with open('bound_template/discharge/input/sim.cfg', 'w') as configfile:
            parser.write(configfile)
        
        
        parser.read('bound_template/restrain/input/sim.cfg')
        parser.set('HEAD', 'boresch restraints dictionary', boresch_restrains)
        with open('bound_template/restrain/input/sim.cfg', 'w') as configfile:
            parser.write(configfile)

        parser.read('bound_template/vanish/input/sim.cfg')
        parser.set('HEAD', 'boresch restraints dictionary', boresch_restrains)
        with open('bound_template/vanish/input/sim.cfg', 'w') as configfile:
            parser.write(configfile)


