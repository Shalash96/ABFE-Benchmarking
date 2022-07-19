# Hydrogen Bond Analysis
I used cpptraj, which is the main program in Amber for processing coordinate trajectories and data files, to carry on the h bond analysis. 

The ***hbond.traj*** file contain the commands that are used for the analysis process

To run the analysis in your terminal type `cpptraj hbond.traj` then hit enter to start the analysis. Make sure you have amber tools installed on your machine, and also edit the ***hbond.traj*** to read your own system that you want to analyse.

After the analysis is done, a lot of files will be produced that contain information about hydrogen bonds formed between solute-solute and solute-solvent. I am interested in the information related to the LIG only, so I wrote a bash script to extract the info related to LIG from the files and put them in new txt files. Just run `./findLIG.sh` from your terminal


Click [here](https://amberhub.chpc.utah.edu/hydrogen-bond-analysis-within-a-protein/) for the full tutorial on how to analyse the h bond
