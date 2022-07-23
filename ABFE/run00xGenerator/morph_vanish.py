import os, shutil, sys, re, glob, threading, time
import subprocess as subp
#Check this package:  (should be parmed: import parmed.amber /from parmed.amber import *   )
from parmed.amber import *

def vanish(prm7, rst7):

     base = AmberParm(prm7, rst7)
     #Strip other molecules. E.g. if you have water: base.strip("WAT") 
     atoms = base.residues[0].atoms
     #Save morph
     fout = open("MORPH.vanish.pert", "w")

     fout.write("version 1\n molecule %s\n" %  base.residues[0].name)

     for atom in atoms:

         name = atom.name
         attype = atom.type
         sigma  = atom.sigma
         eps    = atom.epsilon


         fout.write("    atom\n")
         fout.write("            name %s\n" % name)
         fout.write("            initial_type    %s\n" % attype)
         fout.write("            final_type      %s\n" % "du")
         fout.write("            initial_charge  0.00000\n")
         fout.write("            final_charge    %s\n" % "0.00000")
         fout.write("            initial_LJ      %.5f  %.5f\n" %(sigma,eps))
         fout.write("            final_LJ        %s  %s\n" %("0.00000","0.00000"))
         fout.write("    endatom\n")




     fout.write("end molecule")
