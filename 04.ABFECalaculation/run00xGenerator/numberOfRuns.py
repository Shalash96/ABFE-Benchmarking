import shutil

# check the the state and then call the appropriate function
def numberOfRuns(number=5, state='all'):
    if state == 'free':
        for i in range(number):
            # copy the whole folder
            shutil.copytree("./free_template", f'./free/run00{i+1}')

    elif state == 'bound':
        for i in range(number):
            # copy the whole folder
            shutil.copytree("./bound_template", f'./bound/run00{i+1}')
            
    elif state == 'all':
        for i in range(number):
            # copy the whole folder
            shutil.copytree("./free_template", f'./free/run00{i+1}')
            shutil.copytree("./bound_template", f'./bound/run00{i+1}')

