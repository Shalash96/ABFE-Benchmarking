import shutil


# copy file from one folder to another and rename it
def free(file):
     # check the file extension
            if file.endswith(".prm7"):
                shutil.copy(file, "./free_template/discharge/input/SYSTEM.top")    
                shutil.copy(file, "./free_template/vanish/input/SYSTEM.top")
            
            elif file.endswith(".rst7"):
                shutil.copy(file, "./free_template/discharge/input/SYSTEM.crd")    
                shutil.copy(file, "./free_template/vanish/input/SYSTEM.crd")
            

def bound (file):
    if file.endswith(".prm7"):
        shutil.copy(file, "./bound_template/discharge/input/SYSTEM.top")    
        shutil.copy(file, "./bound_template/vanish/input/SYSTEM.top")
        shutil.copy(file, "./bound_template/restrain/input/SYSTEM.top")
    
    elif file.endswith(".rst7"):
        shutil.copy(file, "./bound_template/discharge/input/SYSTEM.crd")    
        shutil.copy(file, "./bound_template/vanish/input/SYSTEM.crd")
        shutil.copy(file, "./bound_template/restrain/input/SYSTEM.crd")


# fuction to check the state and call the appropriate function
def fileTransfer(src=[], state='all'):
    for file in src:
        # check if state is free or bound or all
        if state == 'free':
            free(file)
        elif state == 'bound':
            bound(file)

        elif state == 'all':
            free(file)
            bound(file)
