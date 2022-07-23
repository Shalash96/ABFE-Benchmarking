# Run00XGenerator 
In our group they create the files and folder required to run the ABFE calculation manually, and they are good in gooing this, but me on the other side did some stupid annoying mistakes such as forgetting to change the boresch restraints values in the .cfg files or forget to change the name of the systems to match the path in .cfg files, etc... . All of these costed me a lot of time to spot the error also much time of running wrong calculation. So i decided to automate the whole process. I write some python scripts that just need you to pass the path of prm7 and rst7 through an interactive command line the
 - It first ask u for the prm7 file
 - Then for the rst7 file
 - Then the mode of the calculations, bound, free, or all
 - Then the number of runs u want
 - Then will ask you if u have a txt file that contains ur boresch restraints if not u will pass the nc file to calaculate it for you


After you submit these files the whole process will be carried out for you, where
 - The prm7 and rst7 files will be copied to the input folder of vanish, discharge and restrain with the name of SYSTEM
 - the morphs files will be created and transferred to the crossponding correct folder
 - Boresch restraints values willbe calculated and save in txt file
 - Also the graphs and charts from the .ipynb file will be saved as png, as i changed it to run as a normal .py file
 - Also there was a small bug in this script where u have everytime to change the value of iteration loop for plotting the graphs by hand to match the same number of the frames in the nc file, so i made it to change automatically.
 - In case of bound or all the cfg file will be edited and updated by the new boresch restraints values that are saved to the txt file or from the txt file that u passed throught the terminal 
 - Then finally the templates folder will be copied to the numbers of runs you specified
 - Then the run_all.sh file will be copied to the free and bound folders that contain the run00x folders. 

And it's simple as running `python run00xGenerator.py` and follow the instructions. 
