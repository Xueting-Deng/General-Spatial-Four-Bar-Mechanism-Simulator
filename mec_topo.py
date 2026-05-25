import numpy as np
import copy

#mec_order_list = ["rscr","rsur","rups","rrus","rrsc","rsup","rusr","rrcs","rsrc","rrsu","rurs","rsru","rusp","rsup","ruuu","rcsr","rcrs","rpus","rpsu","rscp","rspc","rcsp","rcps","rpcs","rpsc","rccc"]
###########################################################################################################################
JJ_4 =  np.array([
    [2, 1, 0, 1],
    [1, 0, 1, 0],
    [0, 1, 0, 1],
    [1, 0, 1, 2],
], dtype=int)

JC_4 = np.array([
    [0, 0, 0], #j0
    [0, 0, 0], #j1
    [0, 0, 0], #j2
    [0, 0, 0], #j3
], dtype=np.float64)

JA_4 = np.array([
    [[0, 0, 1], [0, 0, 0]], # R
    [[0, 0, 0], [0, 0, 0]],
    [[0, 0, 0], [0, 0, 0]],
    [[0, 0, 0], [0, 0, 0]],
], dtype=np.float64)
###########################################################################################################################
JJ_5 = np.array([
    [2, 1, 0, 0, 1],
    [1, 0, 1, 0, 0],
    [0, 1, 0, 1, 0],
    [0, 0, 1, 0, 1],
    [1, 0, 0, 1, 2],
], dtype=int)

JC_5 = np.array([
    [0, 0, 0], #j0
    [0, 0, 0], #j1
    [0, 0, 0], #j2
    [0, 0, 0], #j3
    [0, 0, 0], #j4    
], dtype=np.float64)

JA_5 = np.array([
    [[0, 0, 1], [0, 0, 0]], # R
    [[0, 0, 0], [0, 0, 0]],
    [[0, 0, 0], [0, 0, 0]],
    [[0, 0, 0], [0, 0, 0]],
    [[0, 0, 0], [0, 0, 0]],
], dtype=np.float64)
###########################################################################################################################
JJ_6 = np.array([
    [2, 1, 0, 0, 0, 1],
    [1, 0, 1, 0, 0, 0],
    [0, 1, 0, 1, 0, 0],
    [0, 0, 1, 0, 1, 0],
    [0, 0, 0, 1, 0, 1],
    [1, 0, 0, 0, 1, 2],
], dtype=int)

JC_6 = np.array([
    [0, 0, 0], #j0
    [0, 0, 0], #j1
    [0, 0, 0], #j2
    [0, 0, 0], #j3
    [0, 0, 0], #j4
    [0, 0, 0], #j5    
], dtype=np.float64)

JA_6 = np.array([
    [[0, 0, 1], [0, 0, 0]], # R
    [[0, 0, 0], [0, 0, 0]],
    [[0, 0, 0], [0, 0, 0]],
    [[0, 0, 0], [0, 0, 0]],
    [[0, 0, 0], [0, 0, 0]],
    [[0, 0, 0], [0, 0, 0]],
], dtype=np.float64)

JJ_7 = np.array([
    [2, 1, 0, 0, 0, 0, 1],
    [1, 0, 1, 0, 0, 0, 0],
    [0, 1, 0, 1, 0, 0, 0],
    [0, 0, 1, 0, 1, 0, 0],
    [0, 0, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 1, 0, 1],
    [1, 0, 0, 0, 0, 1, 2],
], dtype=int)

JC_7 = np.array([
    [0, 0, 0], #j0
    [0, 0, 0], #j1
    [0, 0, 0], #j2
    [0, 0, 0], #j3
    [0, 0, 0], #j4
    [0, 0, 0], #j5  
    [0, 0, 0], #j6  
], dtype=np.float64)

JA_7 = np.array([
    [[0, 0, 1], [0, 0, 0]], # R
    [[0, 0, 0], [0, 0, 0]],
    [[0, 0, 0], [0, 0, 0]],
    [[0, 0, 0], [0, 0, 0]],
    [[0, 0, 0], [0, 0, 0]],
    [[0, 0, 0], [0, 0, 0]],
    [[0, 0, 0], [0, 0, 0]],
], dtype=np.float64)


# R-0, P-1, U-2, C-3, S-4, E (plane)-5, W (welded) -6
Mechanism = {
    "RSCR": {#0
        "JJ" : copy.deepcopy(JJ_5),
        "JT" : np.array([0,4,6,3,0], dtype=int),
        "JC" : copy.deepcopy(JC_5),
        "JA" : copy.deepcopy(JA_5),
    },

    "RSUR": {#1
        "JJ" : copy.deepcopy(JJ_4),
        "JT" : np.array([0,4,2,0], dtype=int),
        "JC" : copy.deepcopy(JC_4),
        "JA" : copy.deepcopy(JA_4),       
    },

    "RUPS": {#2
        "JJ" : copy.deepcopy(JJ_5),
        "JT" : np.array([0,2,6,1,4], dtype=int),
        "JC" : copy.deepcopy(JC_5),
        "JA" : copy.deepcopy(JA_5),       
    },

    "RRUS": {#3
        "JJ" : copy.deepcopy(JJ_4),
        "JT" :np.array([0,0,2,4], dtype=int),
        "JC" : copy.deepcopy(JC_4),
        "JA" : copy.deepcopy(JA_4),       
    },

    "RRSC": {#4
        "JJ" : copy.deepcopy(JJ_5),
        "JT" : np.array([0,0,4,6,3], dtype=int),
        "JC" : copy.deepcopy(JC_5),
        "JA" : copy.deepcopy(JA_5),       
    },

    "RSUP": {#5
        "JJ" : copy.deepcopy(JJ_5),
        "JT" : np.array([0,4,2,6,1], dtype=int),
        "JC" : copy.deepcopy(JC_5),
        "JA" : copy.deepcopy(JA_5),       
    },

    "RUSR": {#6
        "JJ" : copy.deepcopy(JJ_4),
        "JT" : np.array([0,2,4,0], dtype=int),
        "JC" : copy.deepcopy(JC_4),
        "JA" : copy.deepcopy(JA_4),       
    },

    "RRCS": {#7
        "JJ" : copy.deepcopy(JJ_5),
        "JT" : np.array([0,0,6,3,4], dtype=int),
        "JC" : copy.deepcopy(JC_5),
        "JA" : copy.deepcopy(JA_5),       
    },    

    "RSRC": {#8
        "JJ" : copy.deepcopy(JJ_5),
        "JT" : np.array([0,4,0,6,3], dtype=int),
        "JC" : copy.deepcopy(JC_5),
        "JA" : copy.deepcopy(JA_5),       
    },    

    "RRSU": {#9
        "JJ" : copy.deepcopy(JJ_4),
        "JT" : np.array([0,0,4,2], dtype=int),
        "JC" : copy.deepcopy(JC_4),
        "JA" : copy.deepcopy(JA_4),       
    },    

    "RURS": {#10
        "JJ" : copy.deepcopy(JJ_4),
        "JT" : np.array([0,2,0,4], dtype=int),
        "JC" : copy.deepcopy(JC_4),
        "JA" : copy.deepcopy(JA_4),       
    },
   
    "RSRU": {#11
        "JJ" : copy.deepcopy(JJ_4),
        "JT" : np.array([0,4,0,2], dtype=int),
        "JC" : copy.deepcopy(JC_4),
        "JA" : copy.deepcopy(JA_4),       
    },   

    "RUSP": {#12
        "JJ" : copy.deepcopy(JJ_5),
        "JT" : np.array([0,2,4,6,1], dtype=int),
        "JC" : copy.deepcopy(JC_5),
        "JA" : copy.deepcopy(JA_5),       
    },    

    "RSPU": {#13
        "JJ" : copy.deepcopy(JJ_5),
        "JT" : np.array([0,4,6,1,2], dtype=int),
        "JC" : copy.deepcopy(JC_5),
        "JA" : copy.deepcopy(JA_5),       
    },    

    "RUUU": {#14
        "JJ" : copy.deepcopy(JJ_4),
        "JT" : np.array([0,2,2,2], dtype=int),
        "JC" : copy.deepcopy(JC_4),
        "JA" : copy.deepcopy(JA_4),       
    },    

    "RCSR": {#15
        "JJ" : copy.deepcopy(JJ_5),
        "JT" : np.array([0,6,3,4,0], dtype=int),
        "JC" : copy.deepcopy(JC_5),
        "JA" : copy.deepcopy(JA_5),       
    },    

    "RCRS": {#16
        "JJ" : copy.deepcopy(JJ_5),
        "JT" : np.array([0,6,3,0,4], dtype=int),
        "JC" : copy.deepcopy(JC_5),
        "JA" : copy.deepcopy(JA_5),       
    },    

    "RPUS": {#17
        "JJ" : copy.deepcopy(JJ_5),
        "JT" : np.array([0, 6, 1, 2, 4], dtype=int),
        "JC" : copy.deepcopy(JC_5),
        "JA" : copy.deepcopy(JA_5),       
    },    

    "RPSU": {#18
        "JJ" : copy.deepcopy(JJ_5),
        "JT" : np.array([0, 6, 1, 4, 2], dtype=int),
        "JC" : copy.deepcopy(JC_5),
        "JA" : copy.deepcopy(JA_5),       
    },    

    "RSCP": {#19
        "JJ" : copy.deepcopy(JJ_6),
        "JT" : np.array([0,4,6,3,6,1], dtype=int),
        "JC" : copy.deepcopy(JC_6),
        "JA" : copy.deepcopy(JA_6),       
    },    

    "RSPC": {#20
        "JJ" : copy.deepcopy(JJ_6),
        "JT" : np.array([0,4,6,1,6,3], dtype=int),
        "JC" : copy.deepcopy(JC_6),
        "JA" : copy.deepcopy(JA_6),       
    },    

    "RCSP": {#21
        "JJ" : copy.deepcopy(JJ_6),
        "JT" : np.array([0,6,3,4,6,1], dtype=int),
        "JC" : copy.deepcopy(JC_6),
        "JA" : copy.deepcopy(JA_6),       
    },    

    "RCPS": {#22
        "JJ" : copy.deepcopy(JJ_6),
        "JT" : np.array([0,6,3,6,1,4], dtype=int),
        "JC" : copy.deepcopy(JC_6),
        "JA" : copy.deepcopy(JA_6),       
    },  

    "RPCS": {#23
        "JJ" : copy.deepcopy(JJ_6),
        "JT" : np.array([0,6,1,6,3,4], dtype=int),
        "JC" : copy.deepcopy(JC_6),
        "JA" : copy.deepcopy(JA_6),       
    },    

    "RPSC": {#24
        "JJ" : copy.deepcopy(JJ_6),
        "JT" : np.array([0,6,1,4,6,3], dtype=int),
        "JC" : copy.deepcopy(JC_6),
        "JA" : copy.deepcopy(JA_6),       
    },    

    "RCCC": {#25
        "JJ" : copy.deepcopy(JJ_7),
        "JT" : np.array([0,6,3,6,3,6,3], dtype=int),
        "JC" : copy.deepcopy(JC_7),
        "JA" : copy.deepcopy(JA_7),       
    },    

}

def update_jc_ja(mec_name, JC, JA, Joint1, Joint2, Joint4, rotational_axis, rotational_axis_vet, m, k, n):
    if mec_name == "RSCR": #0
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint2[m] #j3
        JC[4,:] = Joint4[m] #j4
        JA[3,0,:] = rotational_axis[k] # u3
        JA[4,0,:] = rotational_axis[n] # u4
        save_inx = np.array([1, 2, 4, 3, 4])

    elif mec_name == "RSUR": #1
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint4[m] #j3
        JA[2,0,:] = rotational_axis[k] #u2
        JA[2,1,:] = rotational_axis_vet[k] #u2_2
        JA[3,0,:] = rotational_axis[n] #u3
        save_inx = np.array([1, 2, 3, 2, 3]) #" Why not 23 here???"

    elif mec_name == "RUPS": #2
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint2[m] #j3=j2
        JC[4,:] = Joint4[m] #j4
        JA[1,0,:] = rotational_axis[k] #u1_1
        JA[1,1,:] = rotational_axis_vet[k] #u1_2
        JA[3,0,:] = rotational_axis[n] #u3_1
        save_inx = np.array([1, 2, 4, 1, 3])

    elif mec_name == "RRUS": #3
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint4[m] #j3
        JA[1,0,:] = rotational_axis[k] #u1
        JA[2,0,:] = rotational_axis[n] #u2
        JA[2,1,:] = rotational_axis_vet[n] #u2_2
        save_inx = np.array([1, 2, 3, 1, 2])

    elif mec_name == "RRSC": #4
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint4[m] #j3=j4
        JC[4,:] = Joint4[m] #j3
        JA[1,0,:] = rotational_axis[k] #u1
        JA[4,0,:] = rotational_axis[n] #u4
        save_inx = np.array([1, 2, 3, 1, 4])

    elif mec_name == "RSUP": #5
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint4[m] #j3=j4
        JC[4,:] = Joint4[m] #j4
        JA[2,0,:] = rotational_axis[k] #u2_1
        JA[2,1,:] = rotational_axis_vet[k] #u2_2
        JA[4,0,:] = rotational_axis[n] #u4
        save_inx = np.array([1, 2, 3, 2, 4])

    elif mec_name == "RUSR": #6
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint4[m] #j3=j4
        JA[1,0,:] = rotational_axis[k] #u1_1
        JA[1,1,:] = rotational_axis_vet[k] #u1_2
        JA[3,0,:] = rotational_axis[n] #u3
        save_inx = np.array([1, 2, 3, 1, 3])

    elif mec_name == "RRCS": #7
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint2[m] #j3=j2
        JC[4,:] = Joint4[m] #j4
        JA[1,0,:] = rotational_axis[k] #u1_1
        JA[3,0,:] = rotational_axis[n] #u3
        save_inx = np.array([1, 2, 4, 1, 3])

    elif mec_name == "RSRC": #8
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint4[m] #j3
        JC[4,:] = Joint4[m] #j4=j3
        JA[2,0,:] = rotational_axis[k] #u2
        JA[4,0,:] = rotational_axis[n] #u4
        save_inx = np.array([1, 2, 3, 2, 4])

    elif mec_name == "RRSU": #9
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint4[m] #j3
        JA[1,0,:] = rotational_axis[k] # u1
        JA[3,0,:] = rotational_axis[n] # u3
        JA[3,1,:] = rotational_axis_vet[n] #u3_2
        save_inx = np.array([1, 2, 3, 1, 3])

    elif mec_name == "RURS": #10
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint4[m] #j3
        JA[1,0,:] = rotational_axis[k] # u1
        JA[1,1,:] = rotational_axis_vet[k] #u1_2
        JA[2,0,:] = rotational_axis[n] # u2
        save_inx = np.array([1, 2, 3, 1, 2])

    elif mec_name == "RSRU": #11
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint4[m] #j3
        JA[2,0,:] = rotational_axis[k] # u2
        JA[3,0,:] = rotational_axis[n] # u3
        JA[3,1,:] = rotational_axis_vet[n] #u3_2
        save_inx = np.array([1, 2, 3, 2, 3])

    elif mec_name == "RUSP": #12
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint4[m] #j3
        JC[4,:] = Joint4[m] #j4
        JA[1,0,:] = rotational_axis[k] # u1
        JA[1,1,:] = rotational_axis_vet[k] #u1_2
        JA[4,0,:] = rotational_axis[n] # u4
        save_inx = np.array([1, 2, 3, 1, 4])

    elif mec_name == "RSPU": #13
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint2[m] #j3=j2
        JC[4,:] = Joint4[m] #j4
        JA[3,0,:] = rotational_axis[k] # u3
        JA[4,0,:] = rotational_axis[n] # u4
        JA[4,1,:] = rotational_axis_vet[n] # u4
        save_inx = np.array([1, 2, 4, 3, 4])

    elif mec_name == "RUUU": #14
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint4[m] #j3
        JA[1,0,:] = rotational_axis[k] # u1
        JA[1,1,:] = rotational_axis_vet[k] # u1_2
        JA[2,0,:] = rotational_axis[n] # u2
        JA[2,1,:] = rotational_axis_vet[n] # u2_2
        JA[3,0,:] = np.array([0, 0, 1]) # fix u3
        JA[3,1,:] = np.array([0, 1, 0])
        save_inx = np.array([1, 2, 3, 1, 2])

    elif mec_name == "RCSR": #15
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint1[m] #j2=j1
        JC[3,:] = Joint2[m] #j3
        JC[4,:] = Joint4[m] #j4
        JA[2,0,:] = rotational_axis[k] # u2
        JA[4,0,:] = rotational_axis[n] # u4
        save_inx = np.array([1, 3, 4, 2, 4])

    elif mec_name == "RCRS": #16
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint1[m] #j2=j1
        JC[3,:] = Joint2[m] #j3
        JC[4,:] = Joint4[m] #j4
        JA[2,0,:] = rotational_axis[k] # u2
        JA[3,0,:] = rotational_axis[n] # u3
        save_inx = np.array([1, 3, 4, 2, 3])

    elif mec_name == "RPUS": #17
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint1[m] #j2=j1
        JC[3,:] = Joint2[m] #j3
        JC[4,:] = Joint4[m] #j4
        JA[2,0,:] = rotational_axis[k] # u2
        JA[3,0,:] = rotational_axis[n] # u3
        JA[3,1,:] = rotational_axis_vet[n] # u3_2
        save_inx = np.array([1, 3, 4, 2, 3])

    elif mec_name == "RPSU": #18
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint1[m] #j2=j1
        JC[3,:] = Joint2[m] #j3
        JC[4,:] = Joint4[m] #j4
        JA[2,0,:] = rotational_axis[k] # u2
        JA[4,0,:] = rotational_axis[n] # u4
        JA[4,1,:] = rotational_axis_vet[n] # u4_2
        save_inx = np.array([1, 3, 4, 2, 4])

    elif mec_name == "RSCP": #19
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint2[m] #j3=j2
        JC[4,:] = Joint4[m] #j4
        JC[5,:] = Joint4[m] #j5=j4
        JA[3,0,:] = rotational_axis[k] # u3
        JA[5,0,:] = rotational_axis[n] # u5
        save_inx = np.array([1, 2, 4, 3, 5])

    elif mec_name == "RSPC": #20
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint2[m] #j2
        JC[3,:] = Joint2[m] #j3=j2
        JC[4,:] = Joint4[m] #j4
        JC[5,:] = Joint4[m] #j4=j4
        JA[3,0,:] = rotational_axis[k] # u3
        JA[5,0,:] = rotational_axis[n] # u5
        save_inx = np.array([1, 2, 4, 3, 5])

    elif mec_name == "RCSP": #21
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint1[m] #j2=j1
        JC[3,:] = Joint2[m] #j3
        JC[4,:] = Joint4[m] #j4
        JC[5,:] = Joint4[m] #j5=j4       
        JA[2,0,:] = rotational_axis[k] # u2
        JA[5,0,:] = rotational_axis[n] # u5
        save_inx = np.array([1, 3, 4, 2, 5])

    elif mec_name == "RCPS": #22
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint1[m] #j2=j1
        JC[3,:] = Joint2[m] #j3
        JC[4,:] = Joint2[m] #j4=j3
        JC[5,:] = Joint4[m] #j5    
        JA[2,0,:] = rotational_axis[k] # u2
        JA[4,0,:] = rotational_axis[n] # u4
        save_inx = np.array([1, 3, 5, 2, 4])

    elif mec_name == "RPCS": #23
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint1[m] #j2=j1
        JC[3,:] = Joint2[m] #j3
        JC[4,:] = Joint2[m] #j4=j3
        JC[5,:] = Joint4[m] #j5    
        JA[2,0,:] = rotational_axis[k] # u2
        JA[4,0,:] = rotational_axis[n] # u4
        save_inx = np.array([1, 3, 5, 2, 4])
        
    elif mec_name == "RPSC": #24
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint1[m] #j2=j1
        JC[3,:] = Joint2[m] #j3
        JC[4,:] = Joint4[m] #j4
        JC[5,:] = Joint4[m] #j5=j4    
        JA[2,0,:] = rotational_axis[k] # u2
        JA[5,0,:] = rotational_axis[n] # u5
        save_inx = np.array([1, 3, 4, 2, 5])
        
    elif mec_name == "RCCC": #25
        JC[1,:] = Joint1[m] #j1
        JC[2,:] = Joint1[m] #j2=j1
        JC[3,:] = Joint2[m] #j3
        JC[4,:] = Joint2[m] #j4=j3
        JC[5,:] = Joint4[m] #j5    
        JC[6,:] = Joint4[m] #j6
        JA[2,0,:] = rotational_axis[k] # u2
        JA[4,0,:] = rotational_axis[n] # u4

        JA[6,0,:] = np.array([0, 0, 1]) # fix u3
        save_inx = np.array([1, 3, 5, 2, 4])
        
    else:
        print("The input mec is not in the dictionary")
    return save_inx

def recover_jc_ja(mec_name, JC, JA, mec_t):
    Joint1, Joint2, Joint4 = mec_t[0], mec_t[1], mec_t[2]
    rotational_axis_1, rotational_axis_2 = mec_t[3], mec_t[4]
    rotational_axis_1 = rotational_axis_1 / np.linalg.norm(rotational_axis_1)
    rotational_axis_2 = rotational_axis_2 / np.linalg.norm(rotational_axis_2)

    rotational_axis_1_vet =  np.cross(rotational_axis_1,np.array([1,0,0]))
    rotational_axis_1_vet = rotational_axis_1_vet / np.linalg.norm(rotational_axis_1_vet)

    rotational_axis_2_vet =  np.cross(rotational_axis_2,np.array([1,0,0]))
    rotational_axis_2_vet = rotational_axis_2_vet / np.linalg.norm(rotational_axis_2_vet)

    # rotational_axis_1_vet =  np.cross(rotational_axis_1,np.array([rotational_axis_1[0],0,rotational_axis_1[2]]))
    # rotational_axis_1_vet = rotational_axis_1_vet / np.linalg.norm(rotational_axis_1_vet)

    # rotational_axis_2_vet =  np.cross(rotational_axis_2,np.array([rotational_axis_2[0],0,rotational_axis_2[2]]))
    # rotational_axis_2_vet = rotational_axis_2_vet / np.linalg.norm(rotational_axis_2_vet)

    if mec_name == "RSCR": #0
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint2 #j3
        JC[4,:] = Joint4 #j4
        JA[3,0,:] = rotational_axis_1 # u3
        JA[4,0,:] = rotational_axis_2  # u4

    elif mec_name == "RSUR": #1
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint4 #j3
        JA[2,0,:] = rotational_axis_1 #u2
        JA[2,1,:] = rotational_axis_1_vet #u2_2
        JA[3,0,:] = rotational_axis_2 #u3

    elif mec_name == "RUPS": #2
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint2 #j3=j2
        JC[4,:] = Joint4 #j4
        JA[1,0,:] = rotational_axis_1 #u1_1
        JA[1,1,:] = rotational_axis_1_vet #u1_2
        JA[3,0,:] = rotational_axis_2 #u3_1

    elif mec_name == "RRUS": #3
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint4 #j3
        JA[1,0,:] = rotational_axis_1 #u1
        JA[2,0,:] = rotational_axis_2 #u2
        JA[2,1,:] = rotational_axis_2_vet #u2_2

    elif mec_name == "RRSC": #4
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint4 #j3=j4
        JC[4,:] = Joint4 #j3
        JA[1,0,:] = rotational_axis_1 #u1
        JA[4,0,:] = rotational_axis_2 #u4

    elif mec_name == "RSUP": #5
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint4 #j3=j4
        JC[4,:] = Joint4 #j4
        JA[2,0,:] = rotational_axis_1 #u2_1
        JA[2,1,:] = rotational_axis_1_vet #u2_2
        JA[4,0,:] = rotational_axis_2 #u4

    elif mec_name == "RUSR": #6
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint4 #j3=j4
        JA[1,0,:] = rotational_axis_1 #u1_1
        JA[1,1,:] = rotational_axis_1_vet #u1_2
        JA[3,0,:] = rotational_axis_2 #u3
        
    elif mec_name == "RRCS": #7
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint2 #j3=j2
        JC[4,:] = Joint4 #j4
        JA[1,0,:] = rotational_axis_1 #u1_1
        JA[3,0,:] = rotational_axis_2 #u3

    elif mec_name == "RSRC": #8
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint4 #j3
        JC[4,:] = Joint4 #j4=j3
        JA[2,0,:] = rotational_axis_1 #u2
        JA[4,0,:] = rotational_axis_2 #u4

    elif mec_name == "RRSU": #9
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint4 #j3
        JA[1,0,:] = rotational_axis_1 # u1
        JA[3,0,:] = rotational_axis_2 # u3
        JA[3,1,:] = rotational_axis_2_vet #u3_2

    elif mec_name == "RURS": #10
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint4 #j3
        JA[1,0,:] = rotational_axis_1 # u1
        JA[1,1,:] = rotational_axis_1_vet #u1_2
        JA[2,0,:] = rotational_axis_2 # u2

    elif mec_name == "RSRU": #11
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint4 #j3
        JA[2,0,:] = rotational_axis_1 # u2
        JA[3,0,:] = rotational_axis_2 # u3
        JA[3,1,:] = rotational_axis_2_vet #u3_2

    elif mec_name == "RUSP": #12
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint4 #j3
        JC[4,:] = Joint4 #j4
        JA[1,0,:] = rotational_axis_1 # u1
        JA[1,1,:] = rotational_axis_1_vet #u1_2
        JA[4,0,:] = rotational_axis_2 # u4

    elif mec_name == "RSPU": #13
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint2 #j3=j2
        JC[4,:] = Joint4 #j4
        JA[3,0,:] = rotational_axis_1 # u3
        JA[4,0,:] = rotational_axis_2 # u4
        JA[4,1,:] = rotational_axis_2_vet # u4

    elif mec_name == "RUUU": #14
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint4 #j3
        JA[1,0,:] = rotational_axis_1 # u1
        JA[1,1,:] = rotational_axis_1_vet # u1_2
        JA[2,0,:] = rotational_axis_2 # u2
        JA[2,1,:] = rotational_axis_2_vet # u2_2
        JA[3,0,:] = np.array([0, 0, 1]) # fix u3
        JA[3,1,:] = np.array([0, 1, 0])

    elif mec_name == "RCSR": #15
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint1 #j2=j1
        JC[3,:] = Joint2 #j3
        JC[4,:] = Joint4 #j4
        JA[2,0,:] = rotational_axis_1 # u2
        JA[4,0,:] = rotational_axis_2 # u4

    elif mec_name == "RCRS": #16
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint1 #j2=j1
        JC[3,:] = Joint2 #j3
        JC[4,:] = Joint4 #j4
        JA[2,0,:] = rotational_axis_1 # u2
        JA[3,0,:] = rotational_axis_2 # u3

    elif mec_name == "RPUS": #17
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint1 #j2=j1
        JC[3,:] = Joint2 #j3
        JC[4,:] = Joint4 #j4
        JA[2,0,:] = rotational_axis_1 # u2
        JA[3,0,:] = rotational_axis_2 # u3
        JA[3,1,:] = rotational_axis_2_vet # u3

    elif mec_name == "RPSU": #18
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint1 #j2=j1
        JC[3,:] = Joint2 #j3
        JC[4,:] = Joint4 #j4
        JA[2,0,:] = rotational_axis_1 # u2
        JA[4,0,:] = rotational_axis_2 # u4
        JA[4,1,:] = rotational_axis_2_vet # u4_2

    elif mec_name == "RSCP": #19
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint2 #j3=j2
        JC[4,:] = Joint4 #j4
        JC[5,:] = Joint4 #j5=j4
        JA[3,0,:] = rotational_axis_1 # u3
        JA[5,0,:] = rotational_axis_2 # u5

    elif mec_name == "RSPC": #20
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint2 #j2
        JC[3,:] = Joint2 #j3=j2
        JC[4,:] = Joint4 #j4
        JC[5,:] = Joint4 #j4=j4
        JA[3,0,:] = rotational_axis_1 # u3
        JA[5,0,:] = rotational_axis_2 # u5

    elif mec_name == "RCSP": #21
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint1 #j2=j1
        JC[3,:] = Joint2 #j3
        JC[4,:] = Joint4 #j4
        JC[5,:] = Joint4 #j5=j4       
        JA[2,0,:] = rotational_axis_1 # u2
        JA[5,0,:] = rotational_axis_2 # u5

    elif mec_name == "RCPS": #22
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint1 #j2=j1
        JC[3,:] = Joint2 #j3
        JC[4,:] = Joint2 #j4=j3
        JC[5,:] = Joint4 #j5    
        JA[2,0,:] = rotational_axis_1 # u2
        JA[4,0,:] = rotational_axis_2 # u4

    elif mec_name == "RPCS": #23
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint1 #j2=j1
        JC[3,:] = Joint2 #j3
        JC[4,:] = Joint2 #j4=j3
        JC[5,:] = Joint4 #j5    
        JA[2,0,:] = rotational_axis_1 # u2
        JA[4,0,:] = rotational_axis_2 # u4
        
    elif mec_name == "RPSC": #24
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint1 #j2=j1
        JC[3,:] = Joint2 #j3
        JC[4,:] = Joint4 #j4
        JC[5,:] = Joint4 #j5=j4    
        JA[2,0,:] = rotational_axis_1 # u2
        JA[5,0,:] = rotational_axis_2 # u5
        
    elif mec_name == "RCCC": #25
        JC[1,:] = Joint1 #j1
        JC[2,:] = Joint1 #j2=j1
        JC[3,:] = Joint2 #j3
        JC[4,:] = Joint2 #j4=j3
        JC[5,:] = Joint4 #j5    
        JC[6,:] = Joint4 #j6
        JA[2,0,:] = rotational_axis_1 # u2
        JA[4,0,:] = rotational_axis_2 # u4
        JA[6,0,:] = np.array([0, 0, 1]) # fix u3

    else:
        print("The input mec is not in the dictionary")
    return JC, JA

def Rigid_ref(mec_name, step_para, step_actuated_link):
    if mec_name == "RSCR": #0
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_para[6:9] 
    elif mec_name == "RSUR": #1
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_para[6:9] 
    elif mec_name == "RUPS": #2
        ref_pos_1 = step_para[3:6]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_para[9:12] 
    elif mec_name == "RRUS": #3
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_para[3:6] 
    elif mec_name == "RRSC": #4
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_actuated_link[1,:]
    elif mec_name == "RSUP": #5
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_para[3:6] 
    elif mec_name == "RUSR": #6
        ref_pos_1 = step_para[3:6]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_para[0:3] 
    elif mec_name == "RRCS": #7
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_actuated_link[1,:]
    elif mec_name == "RSRC": #8
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_para[3:6] 
    elif mec_name == "RRSU": #9
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_actuated_link[1,:]
    elif mec_name == "RURS": #10
        ref_pos_1 = step_para[3:6]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_para[6:9] 
    elif mec_name == "RSRU": #11
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_para[3:6] 
    elif mec_name == "RUSP": #12
        ref_pos_1 = step_para[3:6]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_para[0:3] 
    elif mec_name == "RSPU": #13
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_para[6:9] 
    elif mec_name == "RUUU": #14
        ref_pos_1 = step_para[3:6]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_para[0:3] 
    elif mec_name == "RCSR": #15
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_para[3:6]
        ref_axis = step_actuated_link[2,:]
    elif mec_name == "RCRS": #16
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_para[3:6]
        ref_axis = step_actuated_link[2,:]
    elif mec_name == "RPUS": #17
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_para[3:6]
        ref_axis = step_actuated_link[2,:]
    elif mec_name == "RPSU": #18
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_para[3:6]
        ref_axis = step_actuated_link[2,:]
    elif mec_name == "RSCP": #19
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_para[6:9] 
    elif mec_name == "RSPC": #20
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_actuated_link[0,:]
        ref_axis = step_para[6:9] 
    elif mec_name == "RCSP": #21
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_para[3:6]
        ref_axis = step_actuated_link[2,:]
    elif mec_name == "RCPS": #22
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_para[3:6]
        ref_axis = step_actuated_link[2,:]
    elif mec_name == "RPCS": #23
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_para[3:6]
        ref_axis = step_actuated_link[2,:]
    elif mec_name == "RPSC": #24
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_para[3:6]
        ref_axis = step_actuated_link[2,:]
    elif mec_name == "RCCC": #25
        ref_pos_1 = step_para[0:3]
        ref_pos_2 = step_para[3:6]
        ref_axis = step_actuated_link[2,:]
    else:
        print("Label doesn't exist")
        
    return ref_pos_1, ref_pos_2, ref_axis





