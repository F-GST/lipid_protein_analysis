import numpy as np
import mdtraj as md
import itertools
from itertools import *
import matplotlib.pyplot as plt
from scipy import sparse


###### selection parameters, input and output files nd data paths #######
lipid='OPMG'
sel_1 = "resname %s name GL1" %(lipid)
sel_2 = "name SC1"

threshold = 1 # distance threshold in Angstrom


trajectory = '1jb0_mon_native_SC1_GL1_up.xtc' # path to the trajectory file
topology_file = '1jb0_mon_native_SC1_GL1_up.gro' # path to the topology file

#output = 'C:/user/folder/' # path where results will be saved


######## loading input files ########################

traj=md.load(trajectory,top=topology_file)
traj.unitcell_angles = None
traj.unitcell_lengths = None
traj.unitcell_vectors = None	

###################### compute distances matrix ###########################

idx_1 = traj.topology.select(sel_1)    
idx_2 = traj.topology.select(sel_2)
all_lipids=[]
for i in idx_1:
    idx_pairs=[]
    for j in idx_2:
        idx_pairs.append((i,j))
    arr = np.array(idx_pairs)
    distances=md.compute_distances(traj, arr)
    binary_interactions=(distances<threshold).astype(int) # convert distances to binary information/contact map
    trans=binary_interactions.transpose()
    sA = sparse.csr_matrix(trans) 
    print(sA)
    correlation_matrix = np.corrcoef(trans)
    print(correlation_matrix)
    
    all_lipids.append(binary_interactions.transpose())
    
    correlation_matrix = np.corrcoef(binary_interactions.transpose())
    
all_lipids_arr=np.array(all_lipids)
contact_info = coo_matrix((data, (row, col)), shape=(self._nresi_per_protein, ncol_start))
self.interaction_corrcoef = sparse_corrcoef(contact_info)





def compute_lipid_interaction_matrix(traj,topology_file,output,lipid_name):
    top=traj.topology
    all_lipids=traj.xyz[:,top.select("name GL1")]
    prot_xyz_all=traj.xyz[:,top.select("name SC1")]
    frames=len(all_lipids)
    lipid_interaction_matrix=np.zeros((len(#residues),frames))
    for i in range(bframes):
        for j in range(residues): 
            lipid_interaction_matrix[i,j]=get_distance(i,j)
    return lipid_interaction_matrix

def compute_corr_matrix(lipid_interaction_matrix):
    








def compute_correlation_matrix(A,B):
    A = sparse.vstack((A, B), format='csr')
    A = A.astype(np.float64)
    n = A.shape[1]
    # Compute the covariance matrix
    rowsum = A.sum(1)
    centering = rowsum.dot(rowsum.T.conjugate()) / n
    C = (A.dot(A.T.conjugate()) - centering) / (n - 1)
    # The correlation coefficients are given by
    # C_{i,j} / sqrt(C_{i} * C_{j})
    d = np.diag(C)
    corrcoefs = C / np.sqrt(np.outer(d, d))
    return corrcoefs
            
    
    
    
  
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
