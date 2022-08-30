from multiprocessing.connection import wait
import shutil, os
from morph_discharge import discharge
from morph_vanish import vanish
from morph_dummy import dummy

def free(prm7, rst7):
    discharge(prm7, rst7)
    vanish(prm7, rst7)
    shutil.copy("MORPH.discharge.pert", "./free_template/discharge/input/")
    shutil.copy("MORPH.vanish.pert", "./free_template/vanish/input/")
    
def bound(prm7, rst7):
    discharge(prm7, rst7)
    vanish(prm7, rst7)
    dummy(prm7, rst7)
    shutil.copy("MORPH.discharge.pert", "./bound_template/discharge/input/")
    shutil.copy("MORPH.vanish.pert", "./bound_template/vanish/input/")
    shutil.copy("MORPH.dummy.pert", "./bound_template/restrain/input/")

def all(prm7, rst7):
    free(prm7, rst7)
    bound(prm7, rst7)


def handelMORPHs(prm7, rst7, state='all'):

    if state == 'free':
        free(prm7, rst7)
    elif state == 'bound':
        bound(prm7, rst7)
    elif state == 'all':
        free(prm7, rst7)
        bound(prm7, rst7)

