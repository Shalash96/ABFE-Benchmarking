import os, shutil, sys, re, glob, threading, time
import subprocess as subp
#Check this package:  (should be parmed: import parmed.amber /from parmed.amber import *   )
from parmed.amber import *

def discharge(prm7, rst7):

    base = AmberParm(prm7, rst7)
    atoms = base.residues[0].atoms
    natoms = len(atoms)
    fout = open("MORPH.discharge.pert", "w")
    fout.write("version 1\n molecule %s\n" % base.residues[0].name)

    for atom in atoms:

         name   = atom.name
         attype = atom.type
         charge = atom.charge
         sigma  = atom.sigma
         eps    = atom.epsilon


         fout.write("    atom\n")
         fout.write("            name %s\n" % name)
         fout.write("            initial_type    %s\n" % attype)
         fout.write("            final_type      %s\n" % attype)
         if charge < 0:
             fout.write("                initial_charge %.5f\n" % charge)
         else:
             fout.write("                initial_charge  %.5f\n" % charge )

         fout.write("                 final_charge    %.5f\n" % 0.00000)

         fout.write("            initial_LJ      %.5f  %.5f\n" %(sigma,eps))
         fout.write("            final_LJ        %.5f  %.5f\n" %(sigma,eps))
         fout.write("    endatom\n")



    fout.write("end molecule")
