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
mec_label = 22

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


Bi_links, Tri_links, Ground_link, Ground_inx, Actuate_link, Actuate_inx = find_link(JJ, JT, Actuator)
LJ = find_link_joint_table(JJ, Bi_links, Tri_links)
UN = generate_unknown_table(JJ)
Bi_link_len, Bi_link_ang= compute_bi_link_length_angle (LJ, JT, Bi_links, JC, JA)
Tri_link_len, Tri_link_ang = compute_tri_link_length_angle (LJ, JT, Tri_links, JC, JA)
Spacial_p_ang = compute_p_angle(JT, Bi_links, Tri_links, JC, JA)
Type_angle = compute_type_angle(JT, JA)

# Find Constraint Eqs
Updated_UN = assign_value_for_knowns(Ground_link, Ground_inx, Actuate_link, JT, UN, JA, JC, Bi_links) #update UN
UN_actuator_inx = find_actuator_inx (UN, Actuate_link)

Bi_constraint = build_bi_link_constraints (Updated_UN, JT, Bi_links, Bi_link_len, Bi_link_ang, Ground_link)
Tri_constraint= build_tri_link_constraints (Updated_UN, JT, Tri_links, Tri_link_len, Tri_link_ang)
Type_constraint = build_type_angle_constraints (Updated_UN, Type_angle)
Unit_constraint = build_unit_angle_constraints (Updated_UN, JT)
Special_p_constraint = build_p_constraints(JT, Bi_links, Bi_link_len, Tri_links, Tri_link_len, JC, JA, Spacial_p_ang, Updated_UN)
Final_constraint = remove_unvalid_constraint(np.array(Bi_constraint + Tri_constraint+ Type_constraint + Unit_constraint + Special_p_constraint))

# Calculate the initial state
phi = 0
Updated_actuated_link = update_actuated_link(Actuate_link, Actuator, Actuate_inx, JA, JC, JT, Bi_link_len, Tri_link_len,phi)

#Find the initial unknowns and initial guess
Unknown_paras, Initial_guess =find_real_unknowns_and_initial_guess(Ground_link, Ground_inx, Actuate_link, JT, UN, JA, JC, Bi_links, Updated_actuated_link)
Updated_constraint = assign_value_for_actuator(Final_constraint, UN_actuator_inx, Updated_actuated_link)

Initial_para, condition = solve_constraint_equations(Updated_constraint, Unknown_paras, Initial_guess)
if condition == False:
    print("Initial_para doesn't have result, something must be wrong")

phi = 0
step = 360
Step_para = np.zeros((step,len(Unknown_paras))) # store every position of point_4
Step_actuated_link = np.zeros((step,len(Updated_actuated_link),len(Updated_actuated_link[0])))
path_data = np.zeros((27, step, 3))


Step_para[0,:] = Initial_para
Step_actuated_link[0,:] = Updated_actuated_link


for i in range(step-1):
    phi = phi + 360/step
    Updated_actuated_link = update_actuated_link(Actuate_link, Actuator, Actuate_inx, JA, JC, JT, Bi_link_len, Tri_link_len,phi)
    Updated_constraint = assign_value_for_actuator(Final_constraint, UN_actuator_inx, Updated_actuated_link)
    temp_para, condition = solve_constraint_equations(Updated_constraint, Unknown_paras, Step_para[i])
    if condition == True:
        #print('-------------------------')
        Step_para[i+1,:] = temp_para
        Step_actuated_link[i+1,:] = Updated_actuated_link
    else:
        print("No full rotation")
        print("step:",i)
        break
if condition == True:
    print("Mec fully rotated")
for i in range(step):
    Ref_pos_1, Ref_pos_2, Ref_axis = Rigid_ref(mec_name, Step_para[i], Step_actuated_link[i])
    Step_path = RigidbodyMesh(Ref_pos_1, Ref_pos_2, Ref_axis)
    for j in range(27):
        path_data[j, i, :] = Step_path[j,:]

##### plot
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
from vpython import *
scene = canvas(width=800,height=500,center=vector(0,0,0),background=color.white);

an_Joint, an_Link = animation_setting (JT, JC, JA, Bi_links, Tri_links)
an_Update_joint, an_Update_axis = find_unknowns_inx(Ground_link, Ground_inx, Actuate_link, JT, Bi_links)

while True:
    for i in range(0,200):
        #print(i)
        rate(20)
        an_Update_jc, an_Update_ja = animation_update_jc_ja (Actuate_link, an_Update_joint, an_Update_axis, JC, JA, Step_para[i],Step_actuated_link[i])
        #print(an_Update_jc)
        animation_update(an_Joint, an_Link, JT, an_Update_jc, an_Update_ja, Bi_links, Tri_links)