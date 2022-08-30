import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist
from MDAnalysis.lib.distances import calc_dihedrals
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import nglview as nv



def getboresch_restraints(nc_file, top='SYSTEM.top'):
    print('''
************************************************************************
************************************************************************
********************* CALCULATING BORESCH RESTRAINTS *******************
************************************************************************
************************************************************************
''')
    u = mda.Universe(top, nc_file)

    lig_heavy = u.select_atoms("resname LIG and not name H*")

    # anchors dict of dict. For each ligand heavy atom there is a dictionary of protein heavy atoms,
    # for each of which there is a dictionary of average distance and standard deviation

    anchors_dict = {}
    for lig_atom in lig_heavy:
        for prot_atom in u.select_atoms(f"(protein or resname PRT) and (around 10 index {lig_atom.index}) and (not name H*)"): # protein does not recognise PRT
            anchors_dict[(lig_atom.index,prot_atom.index)]={}
            anchors_dict[(lig_atom.index, prot_atom.index)]["dists"]=[]
        


    # Compute average distance and standard deviation for each pair of anchors
    for frame in u.trajectory:
        for lig_atom_index, prot_atom_index in anchors_dict.keys():
            distance = dist(mda.AtomGroup([u.atoms[lig_atom_index]]), mda.AtomGroup([u.atoms[prot_atom_index]]), box=frame.dimensions)[2][0]
            anchors_dict[(lig_atom_index,prot_atom_index)]["dists"].append(distance)




    # change lists to numpy arrays
    for pair in anchors_dict.keys():
        anchors_dict[pair]["dists"] = np.array(anchors_dict[pair]["dists"])




    # calculate average and SD
    for pair in anchors_dict.keys():
        anchors_dict[pair]["avg_dist"] = anchors_dict[pair]["dists"].mean()
        anchors_dict[pair]["sd_dist"] = anchors_dict[pair]["dists"].std()



    # get n pairs with lowest SD and write the values to a txt file
    pairs_ordered_sd=[]
    for item in sorted(anchors_dict.items(), key=lambda item: item[1]["sd_dist"]):
        pairs_ordered_sd.append(item[0])
        # write pair, avg distance and SD to file
        with open("anchors_ordered_sd.txt", "a") as f:
            f.write(f"Pair: {item[0]}, avg_dist: {item[1]['avg_dist']:.2f}, SD: {item[1]['sd_dist']:.2f}\n")
        f.close()




    # For Pairs with Lowest Pairwise RMSDs, find Adjacent Heavy Atoms
    def get_anchor_ats(a1_idx,u):
        """Takes in index of anchor atom 1 and universe and returns
        list of all three anchor atoms, which are chosen to be bonded
        and not H"

        Args:
            a1_idx (int): Index of the first anchor atom
            u (mda universe): The mda universe

        Returns:
            ints: The indices of all three anchor points
        """

        a1_at = u.atoms[a1_idx]
        bonded_heavy_at = a1_at.bonded_atoms.select_atoms("not name H*")
        a2_idx = bonded_heavy_at[0].index

        if len(bonded_heavy_at)>1:
            # not at end of chain
            a3_idx = bonded_heavy_at[1].index
            # Might be better to return all possible combinations
        else:
            # at end of chain, get next heavy atom along
            a3_idx = bonded_heavy_at[0].bonded_atoms.select_atoms("not name H*")[0].index

        return a1_idx, a2_idx, a3_idx



    # Use These As Anchors and Plot Variance of Associated Degrees of Freedom
    def get_distance(idx1, idx2, u):
        """ Distance in Angstrom"""
        distance = dist(mda.AtomGroup([u.atoms[idx1]]), mda.AtomGroup([u.atoms[idx2]]), box=u.dimensions)[2][0]
        return distance


    def get_angle(idx1, idx2, idx3, u):
        """Angle in rad"""
        C = u.atoms[idx1].position 
        B = u.atoms[idx2].position 
        A = u.atoms[idx3].position 
        BA = A - B
        BC = C - B
        angle = np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC)))
        return angle


    def get_dihedral(idx1, idx2, idx3, idx4, u):
        """Dihedral in rad"""
        positions =[u.atoms[idx].position for idx in [idx1,idx2,idx3,idx4]]
        dihedral = calc_dihedrals(positions[0], positions[1], positions[2], positions[3], box = u.dimensions)
        return dihedral


    def get_boresch_dof(l1,l2,l3,r1,r2,r3,u):
        """Calculate Boresch degrees of freedom from indices of anchor atoms"""
        # Ordering of connection of anchors is r3,r2,r1,l1,l2,l3
        r = get_distance(r1,l1,u)
        thetaA = get_angle(r2,r1,l1,u)
        thetaB = get_angle(r1,l1,l2,u)
        phiA = get_dihedral(r3,r2,r1,l1,u)
        phiB = get_dihedral(r2,r1,l1,l2,u)
        phiC = get_dihedral(r1,l1,l2,l3,u)
        # Not restrained but distance from coolinearity must be checked
        thetaR = get_angle(r3,r2,r1,u) # Receptor internal angle
        thetaL = get_angle(l1,l2,l3,u) # Ligand internal angle
        return r, thetaA, thetaB, phiA, phiB, phiC, thetaR, thetaL


    lig_anchors = get_anchor_ats(pairs_ordered_sd[0][0],u)
    prot_anchors = get_anchor_ats(pairs_ordered_sd[0][1],u)
    get_boresch_dof(lig_anchors[0],lig_anchors[1],lig_anchors[2],prot_anchors[0],prot_anchors[1],prot_anchors[2],u)


    # get values of degrees of freedom for lowest SD pairs across whole trajectory

    boresch_dof_dict = {}
    for pair in pairs_ordered_sd[:200]:
        boresch_dof_dict[pair]={}
        l1_idx, r1_idx = pair
        _, l2_idx, l3_idx = get_anchor_ats(l1_idx,u)
        _, r2_idx, r3_idx = get_anchor_ats(r1_idx,u)
        boresch_dof_dict[pair]["anchor_ats"]=[l1_idx,l2_idx,l3_idx,r1_idx,r2_idx,r3_idx]

        boresch_dof_list = ["r","thetaA","thetaB","phiA","phiB","phiC","thetaR","thetaL"]

        # Add sub dictionaries for each Boresch degree of freedom
        for dof in boresch_dof_list:
            boresch_dof_dict[pair][dof]={}
            boresch_dof_dict[pair][dof]["values"]=[]

        # Populate these dictionaries with values from trajectory
        n_frames = len(u.trajectory)

        for i, frame in enumerate(u.trajectory):
            r, thetaA, thetaB, phiA, phiB, phiC, thetaR, thetaL = get_boresch_dof(l1_idx,l2_idx,l3_idx,r1_idx,r2_idx,r3_idx,u)
            boresch_dof_dict[pair]["r"]["values"].append(r)
            boresch_dof_dict[pair]["thetaA"]["values"].append(thetaA)
            boresch_dof_dict[pair]["thetaB"]["values"].append(thetaB)
            boresch_dof_dict[pair]["phiA"]["values"].append(phiA)
            boresch_dof_dict[pair]["phiB"]["values"].append(phiB)
            boresch_dof_dict[pair]["phiC"]["values"].append(phiC)
            boresch_dof_dict[pair]["thetaR"]["values"].append(thetaR)
            boresch_dof_dict[pair]["thetaL"]["values"].append(thetaL)

            if i == n_frames-1:
                boresch_dof_dict[pair]["tot_var"]=0
                for dof in boresch_dof_list:
                    boresch_dof_dict[pair][dof]["values"]=np.array(boresch_dof_dict[pair][dof]["values"])
                    boresch_dof_dict[pair][dof]["avg"]=boresch_dof_dict[pair][dof]["values"].mean()
                    # For dihedrals, compute variance and mean based on list of values corrected for periodic boundary at 
                    # pi radians, because there is no problem with dihedrals in this region
                    if dof[:3] == "phi":
                        avg = boresch_dof_dict[pair][dof]["avg"]

                        # correct variance - fully rigorous
                        corrected_values_sd = []
                        for val in boresch_dof_dict[pair][dof]["values"]:
                            dtheta = abs(val - avg)
                            corrected_values_sd.append(min(dtheta, 2*np.pi-dtheta))
                        corrected_values_sd = np.array(corrected_values_sd) 
                        boresch_dof_dict[pair][dof]["sd"]=corrected_values_sd.std()

                        # Correct mean (will fail if very well split above and below 2pi)
                        # get middle of interval based on current mean
                        corrected_values_avg=[]
                        periodic_bound = avg - np.pi
                        if periodic_bound < -np.pi:
                            periodic_bound+=2*np.pi
                        # shift vals from below periodic bound to above
                        for val in boresch_dof_dict[pair][dof]["values"]:
                            if val < periodic_bound:
                                corrected_values_avg.append(val+2*np.pi)
                            else:
                                corrected_values_avg.append(val)
                        corrected_values_avg = np.array(corrected_values_avg)
                        mean_corrected = corrected_values_avg.mean()
                        #shift mean back to normal range
                        if mean_corrected > np.pi:
                            boresch_dof_dict[pair][dof]["avg"]=mean_corrected-2*np.pi
                        else:
                            boresch_dof_dict[pair][dof]["avg"]=mean_corrected
                            
                    else:
                        boresch_dof_dict[pair][dof]["sd"]=boresch_dof_dict[pair][dof]["values"].std()
                    # Exclude variance of internal angles as these are not restrained
                    if (dof != "thetaR" and dof != "thetaL"):
                        boresch_dof_dict[pair]["tot_var"]+=boresch_dof_dict[pair][dof]["sd"]**2
                    # Assume Gaussian distributions and calculate force constants for harmonic potentials
                    # so as to reproduce these distributions
                    boresch_dof_dict[pair][dof]["k"]=0.593/(boresch_dof_dict[pair][dof]["sd"]**2) # RT at 289 K is 0.593 kcal mol-1


    # Save dictionary of Boresch DOFs to file
    with open("boresch_dof_dict.txt", "w") as f:
        for key, value in boresch_dof_dict.items():
            f.write(str(key) + " " + str(value) + "\n \n")
        f.close()



    # Plot avg and sd of dof for Boresch dof from top 10 lowest SD pairs
    boresch_dof_list = ["r","thetaA","thetaB","phiA","phiB","phiC","thetaR","thetaL"]
    num_pairs = 5
    n_dof = len(boresch_dof_list)

    fig, axs = plt.subplots(1,n_dof, figsize=(2.6*n_dof,6))
    for i, dof in enumerate(boresch_dof_list):
        for j, pair in enumerate(boresch_dof_dict.keys()):
            axs[i].plot([x for x in range(len(u.trajectory))], boresch_dof_dict[pair][dof]["values"],label=f"Pair {pair}")
            if dof == "r":
                axs[i].set_ylabel("r ($\AA$)")
            else:
                axs[i].set_ylabel(f"{dof} (rad)")
            if j == num_pairs-1:
                break
        axs[i].set_xlabel("Frame No")
        axs[i].legend()
    plt.savefig("boresch_dof_avg_sd.png")


    # Filter, Pick Optimum Degrees of Freedom, and Select Force Constants Based on Variance, Select Equilibrium Values
    # Order pairs according to variance 
    pairs_ordered_boresch_var=[]
    for item in sorted(boresch_dof_dict.items(), key=lambda item: item[1]["tot_var"]):
        pairs_ordered_boresch_var.append(item[0])

    # Filter out r <1, theta >150 or < 30 
    selected_pairs_boresch = []
    for pair in pairs_ordered_boresch_var:
        cond_dist = boresch_dof_dict[pair]["r"]["avg"] > 1
        avg_angles =[]
        #angles = ["thetaA", "thetaB", "thetaR","thetaL"] # also check internal angles
        angles = ["thetaA", "thetaB"] # May also be good to check internal angles
        for angle in angles:
            avg_angles.append(boresch_dof_dict[pair][angle]["avg"])
        cond_angles = list(map(lambda x: (x<2.62 and x >0.52),avg_angles))
        if cond_dist and all(cond_angles):
            selected_pairs_boresch.append(pair)


    # save selected pairs boresch to file
    with open("selected_pairs_boresch.txt", "w") as f:
        for pair in selected_pairs_boresch:
            f.write(str(pair) + "\n")
        f.close()


    # Plot histograms

    fig, axs = plt.subplots(len(selected_pairs_boresch[:3]),6, figsize=(16,4*len(selected_pairs_boresch[:3])))
    for pair_idx, pair in enumerate(selected_pairs_boresch[:3]):
        #pair_idx=0
        for i, dof in enumerate(["r","thetaA","thetaB","phiA","phiB","phiC"]):
            axs[pair_idx][i].hist(boresch_dof_dict[pair][dof]["values"],bins=10,label=f"Pair {pair}")
            axs[pair_idx][i].axvline(x=boresch_dof_dict[pair][dof]["avg"], color='r', linestyle='dashed', linewidth=2,label="average")
            if dof == "r":
                axs[pair_idx][i].set_xlabel("r ($\AA$)")
            else:
                axs[pair_idx][i].set_xlabel(f"{dof} (rad)")
            axs[pair_idx][i].set_ylabel("Num Vals")
            axs[pair_idx][i].legend()
    plt.savefig("boresch_dof_hist.png")


    def print_boresch_params(pair):
        l1 = boresch_dof_dict[pair]["anchor_ats"][0]
        l2 = boresch_dof_dict[pair]["anchor_ats"][1]
        l3 = boresch_dof_dict[pair]["anchor_ats"][2]
        r1 = boresch_dof_dict[pair]["anchor_ats"][3]
        r2 = boresch_dof_dict[pair]["anchor_ats"][4]
        r3 = boresch_dof_dict[pair]["anchor_ats"][5]
        r0 = boresch_dof_dict[pair]["r"]["avg"]
        thetaA0 = boresch_dof_dict[pair]["thetaA"]["avg"]
        thetaB0 = boresch_dof_dict[pair]["thetaB"]["avg"]
        phiA0 = boresch_dof_dict[pair]["phiA"]["avg"]
        phiB0 = boresch_dof_dict[pair]["phiB"]["avg"]
        phiC0 = boresch_dof_dict[pair]["phiC"]["avg"]
        kr = boresch_dof_dict[pair]["r"]["k"]
        kthetaA = boresch_dof_dict[pair]["thetaA"]["k"]
        kthetaB = boresch_dof_dict[pair]["thetaB"]["k"]
        kphiA = boresch_dof_dict[pair]["phiA"]["k"]
        kphiB = boresch_dof_dict[pair]["phiB"]["k"]
        kphiC = boresch_dof_dict[pair]["phiC"]["k"]

        return (f'{{"anchor_points":{{"r1":{r1}, "r2":{r2}, "r3":{r3}, "l1":{l1}, "l2":{l2}, "l3":{l3}}},\
        "equilibrium_values":{{"r0":{r0:.2f}, "thetaA0":{thetaA0:.2f}, "thetaB0":{thetaB0:.2f},"phiA0":{phiA0:.2f}, "phiB0":{phiB0:.2f}, "phiC0":{phiC0:.2f}}},\
        "force_constants":{{"kr":{kr:.2f}, "kthetaA":{kthetaA:.2f}, "kthetaB":{kthetaB:.2f}, "kphiA":{kphiA:.2f}, "kphiB":{kphiB:.2f}, "kphiC":{kphiC:.2f}}}}}')




    the_boresch_parm_used = print_boresch_params(selected_pairs_boresch[0])
    # save the_boresch_parm_used to file
    with open("the_boresch_parm_used.txt", "w") as f:
        f.write(str(the_boresch_parm_used))
        f.close()