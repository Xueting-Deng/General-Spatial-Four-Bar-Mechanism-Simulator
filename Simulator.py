import numpy as np
import copy
import math
import time
import sympy as sp

from Simulator_Functions import *
from mec_topo import *

# NOTE: If you want to generate the mechanism collected in mec_topo, please enter the mec_label or mec_name, and change JC and JA matrix to initialize the joints' coordinates
# NOTE: If you want to generate mechanisms that is not in mec_topo, please define your own JJ, JT, JC, JA. Other topology hasen't been fully tested yet, so you might occur error

# please input the mec_label or mec_name, take reference from mec_typo for 1-dof mec

Root_accuracy = 1e-5 # set up accuracy for the calculated result
mec_label = 1  # take RSUR as an example

name_to_label = {name: i for i, name in enumerate(Mechanism.keys())}
label_to_name = {i: name for name, i in name_to_label.items()}

mec_name = label_to_name[mec_label]
print("__________________________________________")
print("Start:",mec_label, mec_name)
m = Mechanism[mec_name]
Actuator = 0 # Put the actuator at J0
JJ = m["JJ"]
JT = m["JT"]
JC = m["JC"]
JA = m["JA"]

if (JJ != (np.transpose(JJ)) ).all():
    print ("Joint-to-Joint matrix is wrong.")

# The JC and JA below are from an RSUR example
JC = np.array([[ 0, 0, 0],
       [ 1, 2, 0],
       [ 4.2266283,  1.53373432, -1.99758577],
       [ 3,  0,  0]])

JA = np.array([[[ 0.        ,  0.        ,  1.        ],
        [ 0.        ,  0.        ,  0.        ]],

       [[ 0.        ,  0.        ,  0.        ],
        [ 0.        ,  0.        ,  0.        ]],

       [[ 0.4759185 , -0.0929362 ,  0.87456529],
        [-0.        ,  0.99440119,  0.10567063]],

       [[ 0.4759185 , -0.0929362 ,  0.87456529],
        [ 0.        ,  0.        ,  0.        ]]])


# Calculate the configuration
Bi_links, Tri_links, Ground_link, Ground_inx, Actuate_link, Actuate_inx = find_link(JJ, JT, Actuator)
LJ = find_link_joint_table(JJ, Bi_links, Tri_links)
UN = generate_unknown_table(JJ)
Bi_link_len, Bi_link_ang= compute_bi_link_length_angle (LJ, JT, Bi_links, JC, JA)
Tri_link_len, Tri_link_ang = compute_tri_link_length_angle (LJ, JT, Tri_links, JC, JA)
Spacial_p_ang = compute_p_angle(JT, Bi_links, Tri_links, JC, JA)
Type_angle = compute_type_angle(JT, JA)


step = 360
path_data = np.zeros((27, step, 3))

Step_para, Step_actuated_link = Open_curve_simulation (JJ, JT, JC, JA, Actuator, step, Root_accuracy)
for i in range(len(Step_para)):
    Ref_pos_1, Ref_pos_2, Ref_axis = Rigid_ref(mec_name, Step_para[i], Step_actuated_link[i])
    Step_path = RigidbodyMesh(Ref_pos_1, Ref_pos_2, Ref_axis)
    for j in range(27):
        path_data[j, i, :] = Step_path[j,:] 
        
        
##### plot
" Might need to copy the plot functions from Simulator_Functions.py to this script to run correctly "
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
inxt = 0
fig = plt.figure(figsize=(8, 8))
ax = plt.subplot(1, 1, 1, projection='3d')
plotMec(JT, JC, JA, Bi_links, Tri_links, Ground_link, JC[0] ,scale=0.5)
#plotPath(Step_para,color = 'red',linestyle = 'line1')
plotPath(path_data[inxt],color = 'red',linestyle = 'point1')
print(inxt)

##### Animation
" Might need to copy the animation functions from Simulator_Functions.py to this script to run correctly "
from vpython import *
scene = canvas(width=800,height=500,center=vector(0,0,0),background=color.white);

an_Joint, an_Link = animation_setting (JT, JC, JA, Bi_links, Tri_links)
an_Update_joint, an_Update_axis = find_unknowns_inx(Ground_link, Ground_inx, Actuate_link, JT, Bi_links)

while True:
    for i in range(0,len(Step_para)):
        #print(i)
        rate(20)
        an_Update_jc, an_Update_ja = animation_update_jc_ja (Actuate_link, an_Update_joint, an_Update_axis, JC, JA, Step_para[i],Step_actuated_link[i])
        #print(an_Update_jc)
        animation_update(an_Joint, an_Link, JT, an_Update_jc, an_Update_ja, Bi_links, Tri_links)
        