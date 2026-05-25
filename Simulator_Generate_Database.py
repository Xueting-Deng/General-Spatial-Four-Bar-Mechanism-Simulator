import numpy as np
import copy
import math
import time
import sympy as sp


from Simulator_Functions import *
from mec_topo import *
import sys


#mec_label = 0
#start = 0
mec_label = int(sys.argv[1])
start = int(sys.argv[2])
Root_accuracy = 1e-5 # set up accuracy for the calculated result
#num_simulation = 1000
num_simulation = 8000 
mec_expect_num = 1800

name_to_label = {name: i for i, name in enumerate(Mechanism.keys())}
label_to_name = {i: name for name, i in name_to_label.items()}

mec_name = label_to_name[mec_label]
m = Mechanism[mec_name]

print("Start:",mec_label, mec_name, start)

Actuator = 0 # Put the actuator at J0
JJ = m["JJ"]
JT = m["JT"]
JC = m["JC"]
JA = m["JA"]
    
filename_mec = f"/gpfs/scratch/xudeng/{mec_name}_mec_{start}_closed.npy"
filename_path = f"/gpfs/scratch/xudeng/{mec_name}_path_{start}_closed.npy"

##################### Random set Joint ##############################
J2_x = [-1,1]
J2_y = [-1,1]
J3_x = [-1,1]
J3_y = [-1,1]
J4_x = [-1,1]
J4_y = [-1,1]
partition = []
for j2x in J2_x:
    for j2y in J2_y:
        for j3x in J3_x:
            for j3y in J3_y:
                for j4x in J4_x:
                    for j4y in J4_y:
                        partition.append([j2x, j2y, j3x, j3y, j4x, j4y])

Joint1_x = np.random.rand(num_simulation)
Joint1_x = Joint1_x*partition[start][0]

Joint1_y = np.random.rand(num_simulation)
Joint1_y = Joint1_y*partition[start][1]

Joint2_x = np.random.rand(num_simulation)
Joint2_x = Joint2_x*partition[start][2]

Joint2_y = np.random.rand(num_simulation)
Joint2_y = Joint2_y*partition[start][3]

Joint4_x = np.random.rand(num_simulation)
Joint4_x = Joint4_x*partition[start][4]

Joint4_y = np.random.rand(num_simulation)
Joint4_y = Joint4_y*partition[start][5]

Joint4_z = np.hstack((np.random.rand(int(num_simulation/2)),-1*np.random.rand(int(num_simulation/2))))
np.random.shuffle(Joint4_z)

Joint1 = np.zeros((num_simulation,3), dtype=np.float32)
Joint2 = np.zeros((num_simulation,3), dtype=np.float32)
Joint4 = np.zeros((num_simulation,3), dtype=np.float32)
for i in range(num_simulation):
    Joint1[i] = np.array([Joint1_x[i],Joint1_y[i],0], dtype=np.float32)
    Joint2[i] = np.array([Joint2_x[i],Joint2_y[i],0], dtype=np.float32)
    Joint4[i] = np.array([Joint4_x[i],Joint4_y[i],Joint4_z[i]], dtype=np.float32)


def fun(u, p):
    # u: 1D array for unknowns (x27..x32), p: 1D array for params (x9,x10,x11)
    vals = (*u, *p)
    return np.asarray(F_num(*vals), dtype=float).ravel()

# # Calculate the configuration

mec_num = mec_expect_num
temp_count = 0
step = 360 # this step decides how to devide 2pi
mec_data = np.zeros((mec_num*27,6,3), dtype=np.float32) # only store j1, j2, j3, u0, u1, coupler
path_data = np.zeros((mec_num*27,step,3), dtype=np.float32)
print("here")

#path_mm = open_memmap(filename_path, mode='w+', dtype=np.float16, shape=(mec_num*27, step, 3))
#mec_mm  = open_memmap(filename_mec,  mode='w+', dtype=np.float16, shape=(mec_num*27, 6, 3))

############################## Start Generating ##################################
start_time = time.time()
fullyrotated_num = 0

for m in range(num_simulation):
    if fullyrotated_num == mec_expect_num:
        break  

    angle1 = np.random.rand(3)*np.pi
    angle2 = np.random.rand(3)*np.pi
    rotational_axis = []
    rotational_axis_vet = [] # To find a value for u5 that is vertical to u6 and u6's projection on xz plane
    for ang1 in angle1:
        uy = np.cos(ang1)
        for ang2 in angle2:
            ux = np.sin(ang1)*np.cos(ang2)
            uz = np.sin(ang1)*np.sin(ang2)
            rotational_axis.append(np.array([ux,uy,uz]))
            vet = np.cross(np.array([ux,uy,uz]),np.array([1,0,0]))
            vet = vet/np.linalg.norm(vet)
            rotational_axis_vet.append(vet)
    rotational_axis = np.array(rotational_axis, dtype=np.float32)
    rotational_axis_vet = np.array(rotational_axis_vet, dtype=np.float32)

    for k in range(len(rotational_axis)):
        if fullyrotated_num == mec_expect_num:
            break
        for n in range(len(rotational_axis)): 
            if fullyrotated_num == mec_expect_num:
                break
            try:
                temp_count += 1
                # print(temp_count)
               # print("\r",temp_count,end='')
                Save_inx = update_jc_ja(mec_name, JC, JA, Joint1, Joint2, Joint4, rotational_axis, rotational_axis_vet, m, k, n)
                if JT[2] == 4 and JC[2,0] == 0 and JC[2,1] == 0:
                    # and J2 is a spherical joint
                    # print("JC_2",JC[2])
                    # print("JA_2",JA[2])
                    continue                
                elif JC[2,0] == 0 and JC[2,1] == 0 and (JT[2] == 0 or JT[2] == 2 or JT[2] == 3):
                    if (JA[2,0,:]==np.array([0,0,1])).all() or (JA[2,1,:]==np.array([0,0,1])).all():
                    # and J2 is a r/u/c joint, and J2's rotational axis is along z axis
                        # print("JC_2",JC[2])
                        # print("JA_2",JA[2])
                        continue 
                # Calculate the configuration
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

                #try:
                    #Find the initial unknowns and initial guess
                Unknown_paras, Initial_guess =find_real_unknowns_and_initial_guess(Ground_link, Ground_inx, Actuate_link, JT, UN, JA, JC, Bi_links, Updated_actuated_link)
                arglist = [*Unknown_paras, *UN_actuator_inx]
                F_filtered = sp.Matrix([f for f in Final_constraint if f.free_symbols & set(Unknown_paras)])
                F_num = sp.lambdify(arglist, F_filtered, modules = 'numpy') #change the constraint equations from symbolic to numerical
                F_wrapped = lambda u, p: fun(u, p)
                result_temp = root(F_wrapped, Initial_guess, args=(Updated_actuated_link.flatten()), tol = Root_accuracy)
                Initial_para, condition = result_temp.x, result_temp.success
                if condition == False:
                    continue
                phi = 0

                Step_para = np.zeros((step,len(Unknown_paras))) # store every position of point_4
                Step_actuated_link = np.zeros((step,len(Updated_actuated_link),len(Updated_actuated_link[0])))
                Step_para[0,:] = Initial_para
                Step_actuated_link[0,:] = Updated_actuated_link


                for i in range(step-1):
                    phi = phi + 360/step
                    Updated_actuated_link = update_actuated_link(Actuate_link, Actuator, Actuate_inx, JA, JC, JT, Bi_link_len, Tri_link_len,phi)
                    result_temp = root(F_wrapped, Step_para[i], args=(Updated_actuated_link.flatten()), tol = Root_accuracy)
                    temp_para, condition = result_temp.x, result_temp.success    
                    if condition == True:
                        #print('-------------------------')
                        Step_para[i+1,:] = temp_para
                        Step_actuated_link[i+1,:] = Updated_actuated_link
                    else:
    #                     print("step:",i)
                        break
                if condition == True:
                    for i in range(step):
                        Ref_pos_1, Ref_pos_2, Ref_axis = Rigid_ref(mec_name, Step_para[i], Step_actuated_link[i])
                        Step_path = RigidbodyMesh(Ref_pos_1, Ref_pos_2, Ref_axis)
                        for j in range(27):
                            path_data[fullyrotated_num*27+j, i, :] = Step_path[j,:]

                    for j in range(27):
                        mec_data[fullyrotated_num*27+j,0,:] = JC[Save_inx[0],:] # joint1
                        mec_data[fullyrotated_num*27+j,1,:] = JC[Save_inx[1],:] # joint2
                        mec_data[fullyrotated_num*27+j,2,:] = JC[Save_inx[2],:] # joint4
                        mec_data[fullyrotated_num*27+j,3,:] = JA[Save_inx[3],0,:] # u3
                        mec_data[fullyrotated_num*27+j,4,:] = JA[Save_inx[4],0,:] # u4
                        mec_data[fullyrotated_num*27+j,5,:] = path_data[fullyrotated_num*27+j, 0, :] #coupler

                    fullyrotated_num += 1
                    if fullyrotated_num % 100 == 0 or fullyrotated_num == 1:
                        #print("Saved file of",mec_name)
                        #print("Fully rotated machine number: %d" % (fullyrotated_num))
                        #print("Number of attemps: %d" % (temp_count))
                        #print("Running Time:", round((time.time() - start_time)/60,2),"mins")
                        np.save(filename_mec,mec_data)
                        np.save(filename_path,path_data)
                        continue
            except:
                #print("weird thing happend at:", temp_count)
                pass

print("Running Time:", round((time.time() - start_time)/60,2),"mins")
