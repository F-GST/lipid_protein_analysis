#calculate the DE index
# make trajectory that contains only the the GL1 atoms of lipids and the protein 
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import pandas as pd
import sys
import os

trajectory=sys.argv[2]
topology_file=sys.argv[3]
output=sys.argv[1]
traj=md.load(trajectory,top=topology_file)

def get_distance(p, q):
    s_sq_difference = 0
    for p_i,q_i in zip(p,q):
        s_sq_difference += (p_i - q_i)**2
        distance = s_sq_difference**0.5
    return distance

def compute_lipid_interaction_matrix(traj,topology_file,output,lipid_name):
    top=traj.topology
    all_lipids=traj.xyz[:,top.select("name GL1")]
    prot_xyz_all=traj.xyz[:,top.select("name SC1")]
    frames=len(all_lipids)
    lipid_interaction_matrix=np.
    
    
  
    DE_file=open(output,'w')
    DE_file.write("lipid_name DE_0.7 DE_1.4 DE_2.1 ")
    DE_file.write("\n")
    cutoff=[0.7,1.4,2.1]
    all_lipids=traj.xyz[:,top.select("name GL1")]
    prot_xyz_all=traj.xyz[:,top.select("protein")]
    frames=len(all_lipids)
    all_lipid_distances_0_7=0
    all_lipid_distances_1_4=0
    all_lipid_distances_2_1=0  
    
    for i in range(frames):
        print("frame %i"%(i))             
        lipid=all_lipids[i]
        prot=prot_xyz_all[i]
        for liptm in lipid:
            for prottm in prot: 
                distance=get_distance(liptm,prottm)
                if distance<0.7:
                    all_lipid_distances_0_7=all_lipid_distances_0_7+1
                elif distance < 1.4:
                    all_lipid_distances_1_4=all_lipid_distances_1_4+1 
                elif distance < 2.1:
                    all_lipid_distances_2_1=all_lipid_distances_2_1+1


    avg_all_lipid_distances_0_7=all_lipid_distances_0_7/frames
    avg_all_lipid_distances_1_4=all_lipid_distances_1_4/frames
    avg_all_lipid_distances_2_1=all_lipid_distances_2_1/frames  


    for l in range(len(lipids_list)):
        print(lipid_head_name[l])
        lipid_distances_0_7=0
        lipid_distances_2_1=0
        lipid_distances_1_4=0
        splipid=traj.xyz[:,top.select(lipids_list[l])]
        #print(splipid)
        for i in range(frames):
            lipid=splipid[i]
            prot=prot_xyz_all[i]
            for liptm in lipid:
                for prottm in prot:
                    distance=get_distance(liptm,prottm)
                    if distance<0.7:
                        lipid_distances_0_7=lipid_distances_0_7+1
                    elif distance < 1.4:
                        lipid_distances_1_4=lipid_distances_1_4+1
                    elif distance < 2.1:
                        lipid_distances_2_1=lipid_distances_2_1+1
        avg_lipid_distances_0_7=lipid_distances_0_7/frames
        avg_lipid_distances_1_4=lipid_distances_1_4/frames
        avg_lipid_distances_2_1=lipid_distances_2_1/frames
        ratio_l_1_7=avg_lipid_distances_0_7/avg_all_lipid_distances_0_7
        ratio_l_1_4=avg_lipid_distances_1_4/avg_all_lipid_distances_1_4
        ratio_l_2_1=avg_lipid_distances_2_1/avg_all_lipid_distances_2_1
        ratio_l_bulk=len(splipid[0])/len(all_lipids[0])
        DE1_7=ratio_l_1_7/ratio_l_bulk
        DE1_4=ratio_l_1_4/ratio_l_bulk
        DE2_1=ratio_l_2_1/ratio_l_bulk
        print(DE1_7)
        print(DE1_4)
        print(DE2_1)
        
        
        
        DE_file.write("%s %f %f %f \n"%(lipid_head_name[l],DE1_7,DE1_4,DE2_1))


lipids_list=[
"resname OPPG and name GL1",
"resname OPGG and name GL1",
"resname OPMG and name GL1",
"resname OPSG and name GL1",
"resname DPSG and name GL1"
]

lipid_head_name=["OPPG","OPGG","OPMG","OPSG","DPSG"]


#plot 2D lateral density maps


            
computeDE(traj, topology_file,output,lipids_list,lipid_head_name)
