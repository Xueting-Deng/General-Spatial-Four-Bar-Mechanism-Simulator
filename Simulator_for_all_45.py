import numpy as np
import copy
import math
import time
import sympy as sp

from Simulator_Functions import *
from mec_topo import *
import sys

# please input the mec_label or mec_name, take reference from mec_typo for 1-dof mec
mec_label = 22
test_ind = True

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
    

############################## Change ######################################
start = int(sys.argv[1])
end   = int(sys.argv[2])
# start = 0
# end = 1000
filename_mec = f"/gpfs/scratch/xudeng/{mec_name}_mec_test_new.npy"
filename_path = f"/gpfs/scratch/xudeng/{mec_name}_path_test_new.npy"
#     print(filename_mec)
#     print(filename_path)
print(mec_name)
#filename_mec = f"/gpfs/scratch/xudeng/{mec_name}_mec_{start}.npy"
#filename_path = f"/gpfs/scratch/xudeng/{mec_name}_path_{start}.npy"
Joint4 = np.load('/gpfs/home/xudeng/Joint2_7.npy')
Joint1 = np.load('/gpfs/home/xudeng/Joint3_7.npy')
Joint2 = np.load('/gpfs/home/xudeng/Joint4_7.npy')
# Joint4 = np.load('C:\SBU-3\Jupyter-Research\RSCR\\0 Data Generation\Joint2_7.npy')
# Joint1 = np.load('C:\SBU-3\Jupyter-Research\RSCR\\0 Data Generation\Joint3_7.npy')
# Joint2 = np.load('C:\SBU-3\Jupyter-Research\RSCR\\0 Data Generation\Joint4_7.npy')


angle = np.array([30,90,150])*np.pi/180


rotational_axis = []
rotational_axis_vet = [] # To find a value for u5 that is vertical to u6 and u6's projection on xz plane
for ang1 in angle:
    uy = np.cos(ang1)
    for ang2 in angle:
        ux = np.sin(ang1)*np.cos(ang2)
        uz = np.sin(ang1)*np.sin(ang2)
        rotational_axis.append(np.array([ux,uy,uz]))
        vet = np.cross(np.array([ux,uy,uz]),np.array([1,0,0]))
        vet = vet/np.linalg.norm(vet)
        rotational_axis_vet.append(vet)
rotational_axis = np.array(rotational_axis, dtype=np.float32)
rotational_axis_vet = np.array(rotational_axis_vet, dtype=np.float32)



# # Calculate the configuration

mec_num = 45*28
temp_count = 0
step = 360 # this step decides how to devide 2pi
mec_data = np.zeros((mec_num*27,6,3), dtype=np.float32) # only store j1, j2, j3, u0, u1, coupler
path_data = np.zeros((mec_num*27,step,3), dtype=np.float32)

############################## Start Generating ##################################
start_time = time.time()
fullyrotated_num = 0
#print("Did you change rigidbody mesh function?")
for m in range(start,end):
    # print('M:',m)
    if test_ind:
        for k in range(len(rotational_axis)):
            # print("K:",k)
            if test_ind:
                for n in range(len(rotational_axis)):
                    #print("m:",m)
                    #print("k:",k)
                    #print("n:",n)
                    temp_count += 1
                    #print("\r",temp_count,end='')
                    Save_inx = update_jc_ja(mec_name, JC, JA, Joint1, Joint2, Joint4, rotational_axis, rotational_axis_vet, m, k, n)
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

                    try:
                        #Find the initial unknowns and initial guess
                        Unknown_paras, Initial_guess =find_real_unknowns_and_initial_guess(Ground_link, Ground_inx, Actuate_link, JT, UN, JA, JC, Bi_links, Updated_actuated_link)
                        Updated_constraint = assign_value_for_actuator(Final_constraint, UN_actuator_inx, Updated_actuated_link)

                        Initial_para, condition = solve_constraint_equations(Updated_constraint, Unknown_paras, Initial_guess)
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
                            Updated_constraint = assign_value_for_actuator(Final_constraint, UN_actuator_inx, Updated_actuated_link)
                            temp_para, condition = solve_constraint_equations(Updated_constraint, Unknown_paras, Step_para[i])
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
                            if fullyrotated_num % 10 == 0 or fullyrotated_num == 1:
                                print("Fully rotated mechine number: %d" % (fullyrotated_num))
                            if fullyrotated_num == 45:
                                print("Saved file of",mec_name)
                                print("Fully rotated machine number: %d" % (fullyrotated_num))
                                print("Number of attemps: %d" % (temp_count))
                                print("Running Time:", round((time.time() - start_time)/60,2),"mins")
                                # np.save(filename_mec,mec_data)
                                # np.save(filename_path,path_data)
                                test_ind = False
                            #continue
                                break

                # except Exception as e:
                    except:
                        print("weird thing happend at:", temp_count)
                #     print(e)
                        pass
            else:
                break
    else:
            break
########## run end #############
# print("--------------The End of %d---------------------" % (start))
# print("Running Time:", round((time.time() - start_time)/60,2),"mins")
# print("Fully rotated machine number: %d" % (fullyrotated_num))
# print("Number of attemps: %d" % (temp_count))
# np.save(filename_mec,mec_data)
# np.save(filename_path,path_data)
# print('-------------------------------------------------------------')