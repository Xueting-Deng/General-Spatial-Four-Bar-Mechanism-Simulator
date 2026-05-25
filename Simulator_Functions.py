import numpy as np
import copy
from scipy.optimize import root
import math
import time
import sympy as sp

# # Notes:
# 1. when p or c joint exist in link l_012, j1 needs to be w joint. If j0 is p/c joint, then p/c joint belongs to this link. Else if j2 is p/c joint, then p/c joint belongs to the other link that connected to this link

# Function of quaternion
def QuaterTimes(q1, q2):
    a, b, c, d = q1;
    e, f, g, h = q2;
    w = a*e - b*f - c*g - d*h;
    x = b*e + a*f - d*g + c*h;
    y = c*e + d*f + a*g - b*h;
    z = d*e - c*f + b*g + a*h;
    
    return np.array([w, x, y, z])

def QuaterConj(q):
    w, x, y, z = q;
    return np.array([w, -x, -y, -z])

# Rodrigues Rotation Formula
# Work the same with quaternion

def rotate_vector(v, k, angle): #k is the axis of rotation, v is the vector to rotate
    
    # Calculate terms of Rodrigues' formula
    v_rot = (v * np.cos(angle) + 
             np.cross(k, v) * np.sin(angle) + 
             k * np.dot(k, v) * (1 - np.cos(angle)))
    
    return v_rot

def regular_angle(angle):
    angle = angle if abs(angle-1) > 1e-6 else 1
    angle = angle if abs(angle+1) > 1e-6 else -1
    angle = (angle if angle > -np.pi else angle+np.pi) if angle < np.pi else angle-np.pi;#make sure the ang is in (0,pi)
    
    return angle

# make sure each rotational axis is unit vector
def unit_JA(ja):
    ja = copy.deepcopy(ja)
    for i in range(len(ja)):
        for j in range(len(ja[i])):
            #print(ja[i,j])
            if all(ja[i,j] == np.array([0,0,0])):
                continue
            else:
                ja[i,j] = ja[i,j]/np.linalg.norm(ja[i,j])
    return ja

def RigidbodyMesh_old (ref_pos_1, ref_pos_2, ref_axis):
    
    mesh_data = np.zeros((27,3), dtype=np.float64);
    
    center = 0.5*(ref_pos_1-ref_pos_2) + ref_pos_2;
    #x axis
    x_rigid = np.cross((ref_pos_2 - ref_pos_1),ref_axis);
    x_rigid = x_rigid/np.linalg.norm(x_rigid);
    #y axis
    y_rigid = ref_pos_2 - ref_pos_1;
    y_rigid = y_rigid/np.linalg.norm(y_rigid);
    #z axis
    z_rigid = np.cross(x_rigid,y_rigid);
    z_rigid = z_rigid/np.linalg.norm(z_rigid);

    #mesh, unit = 2
    mesh_data[0,:] = center; #[0,0,0]

    mesh_data[1,:] = center + 2*z_rigid; #[0,0,1]
    mesh_data[2,:] = center - 2*z_rigid; #[0,0,-1]
    
    mesh_data[3,:] = center + 2*x_rigid; #[1,0,0]
    mesh_data[4,:] = center - 2*x_rigid; #[-1,0,0]
    
    mesh_data[5,:] = center + 2*y_rigid; #[0,1,0]
    mesh_data[6,:] = center - 2*y_rigid; #[0,-1,0]
    
    
    mesh_data[7,:] = center + 2*x_rigid + 2*y_rigid; #[1,1,0]
    mesh_data[8,:] = center + 2*x_rigid - 2*y_rigid; #[1,-1,0]
    mesh_data[9,:] = center - 2*x_rigid + 2*y_rigid; #[-1,1,0]
    mesh_data[10,:] = center - 2*x_rigid - 2*y_rigid; #[-1,-1,0]
    
    mesh_data[11,:] = center + 2*x_rigid + 2*z_rigid; #[1,0,1]
    mesh_data[12,:] = center + 2*x_rigid - 2*z_rigid; #[1,0,-1]
    mesh_data[13,:] = center - 2*x_rigid + 2*z_rigid; #[-1,0,1]
    mesh_data[14,:] = center - 2*x_rigid - 2*z_rigid; #[-1,0,-1]
    
    mesh_data[15,:] = center + 2*y_rigid + 2*z_rigid; #[0,1,1]
    mesh_data[16,:] = center + 2*y_rigid - 2*z_rigid; #[0,1,-1]
    mesh_data[17,:] = center - 2*y_rigid + 2*z_rigid; #[0,-1,1]
    mesh_data[18,:] = center - 2*y_rigid - 2*z_rigid; #[0,-1,-1]
    
    mesh_data[19,:] = center + 2*x_rigid + 2*y_rigid + 2*z_rigid; #[1,1,1]
    mesh_data[20,:] = center + 2*x_rigid + 2*y_rigid - 2*z_rigid; #[1,1,-1]
    mesh_data[21,:] = center + 2*x_rigid - 2*y_rigid + 2*z_rigid; #[1,-1,1]
    mesh_data[22,:] = center - 2*x_rigid + 2*y_rigid + 2*z_rigid; #[-1,1,1]
    mesh_data[23,:] = center + 2*x_rigid - 2*y_rigid - 2*z_rigid; #[1,-1,-1]
    mesh_data[24,:] = center - 2*x_rigid - 2*y_rigid + 2*z_rigid; #[-1,-1,1]
    mesh_data[25,:] = center - 2*x_rigid + 2*y_rigid - 2*z_rigid; #[-1,1,-1]
    mesh_data[26,:] = center - 2*x_rigid - 2*y_rigid - 2*z_rigid; #[-1,-1,-1]
    
    return mesh_data

def RigidbodyMesh (ref_pos_1, ref_pos_2, ref_axis):
    
    mesh_data = np.zeros((27,3), dtype=np.float64);
    bar_len = np.linalg.norm(ref_pos_2-ref_pos_1)
    #bar_len = 2
    
    center = 0.5*(ref_pos_1-ref_pos_2) + ref_pos_2;
    #x axis
    x_rigid = np.cross((ref_pos_2 - ref_pos_1),ref_axis);
    x_rigid = x_rigid/np.linalg.norm(x_rigid);
    #y axis
    y_rigid = ref_pos_2 - ref_pos_1;
    y_rigid = y_rigid/np.linalg.norm(y_rigid);
    #z axis
    z_rigid = np.cross(x_rigid,y_rigid);
    z_rigid = z_rigid/np.linalg.norm(z_rigid);

    #mesh, unit = 2
    mesh_data[0,:] = center; #[0,0,0]

    mesh_data[1,:] = center + bar_len*z_rigid; #[0,0,1]
    mesh_data[2,:] = center - bar_len*z_rigid; #[0,0,-1]
    
    mesh_data[3,:] = center + bar_len*x_rigid; #[1,0,0]
    mesh_data[4,:] = center - bar_len*x_rigid; #[-1,0,0]
    
    mesh_data[5,:] = center + bar_len*y_rigid; #[0,1,0]
    mesh_data[6,:] = center - bar_len*y_rigid; #[0,-1,0]
    
    
    mesh_data[7,:] = center + bar_len*x_rigid + bar_len*y_rigid; #[1,1,0]
    mesh_data[8,:] = center + bar_len*x_rigid - bar_len*y_rigid; #[1,-1,0]
    mesh_data[9,:] = center - bar_len*x_rigid + bar_len*y_rigid; #[-1,1,0]
    mesh_data[10,:] = center - bar_len*x_rigid - bar_len*y_rigid; #[-1,-1,0]
    
    mesh_data[11,:] = center + bar_len*x_rigid + bar_len*z_rigid; #[1,0,1]
    mesh_data[12,:] = center + bar_len*x_rigid - bar_len*z_rigid; #[1,0,-1]
    mesh_data[13,:] = center - bar_len*x_rigid + bar_len*z_rigid; #[-1,0,1]
    mesh_data[14,:] = center - bar_len*x_rigid - bar_len*z_rigid; #[-1,0,-1]
    
    mesh_data[15,:] = center + bar_len*y_rigid + bar_len*z_rigid; #[0,1,1]
    mesh_data[16,:] = center + bar_len*y_rigid - bar_len*z_rigid; #[0,1,-1]
    mesh_data[17,:] = center - bar_len*y_rigid + bar_len*z_rigid; #[0,-1,1]
    mesh_data[18,:] = center - bar_len*y_rigid - bar_len*z_rigid; #[0,-1,-1]
    
    mesh_data[19,:] = center + bar_len*x_rigid + bar_len*y_rigid + bar_len*z_rigid; #[1,1,1]
    mesh_data[20,:] = center + bar_len*x_rigid + bar_len*y_rigid - bar_len*z_rigid; #[1,1,-1]
    mesh_data[21,:] = center + bar_len*x_rigid - bar_len*y_rigid + bar_len*z_rigid; #[1,-1,1]
    mesh_data[22,:] = center - bar_len*x_rigid + bar_len*y_rigid + bar_len*z_rigid; #[-1,1,1]
    mesh_data[23,:] = center + bar_len*x_rigid - bar_len*y_rigid - bar_len*z_rigid; #[1,-1,-1]
    mesh_data[24,:] = center - bar_len*x_rigid - bar_len*y_rigid + bar_len*z_rigid; #[-1,-1,1]
    mesh_data[25,:] = center - bar_len*x_rigid + bar_len*y_rigid - bar_len*z_rigid; #[-1,1,-1]
    mesh_data[26,:] = center - bar_len*x_rigid - bar_len*y_rigid - bar_len*z_rigid; #[-1,-1,-1]
    
    return mesh_data
# # Calculate the configuration

# Simply check if this mechanism is 1-dof
# Only suitable for usual four-bar!!!
def simple_dof_check(jt, links):
    n = 6*(len(links)-1)
    for joint in jt:
        if joint == 0 or joint == 1:
            n = n-5
        elif joint == 2 or joint == 3:
            n = n-4
        elif joint == 4 or joint == 5:
            n = n-3
        elif joint == 6:
            n = n-6
        else:
            print('wrong type')
            break
    if n == 1:
        print("The DOF of this spatial mechanism is 1")
    else:
        print("The DOF of this spatial mechanism is not 1, this simulator can't solve this mechanism")

# Find both binary and tribary links
def find_link(jj, jt, actuator):
    bi_links = []
    tri_links = []
    ground_link = []
    
    ground_inx = []
    actuated_inx = []
    temp_tri = [0,0,0]
    
    for i in range(len(jj[0])):
        "Below needs to be fixed, for tri_link L012, joint[1] needs to be W joint, and joint[2] needs to be C/P joint"
        if jt[i] == 6:
            temp_tri[1] = i
            temp_inx = 0
            for j in range(len(jj[0])):
                if jj[i,j] == 1 or jj[i,j]==2:
                    temp_tri[temp_inx] = j
                    temp_inx = 2
#             if jt[temp_tri[0]] == 1 or jt[temp_tri[0]] == 3: # make sure tri_link[2] is the special joint
#                 tri_links.append((temp_tri[2],temp_tri[1],temp_tri[0]))
#             else:
#                 tri_links.append((temp_tri[0],temp_tri[1],temp_tri[2]))
            "warning: when both sides are c or p joint, there might be problems" 
            tri_links.append((temp_tri[0],temp_tri[1],temp_tri[2]))      
        else:
            for j in range(i, len(jj[0])):
                if jj[i,j] == 1 and jt[j] != 6:
                    bi_links.append((i,j))
                elif jj[i,j] == 2:
                    ground_link.append(i)
    bi_links = np.array(bi_links, dtype=int)
    tri_links = np.array(tri_links, dtype=int)
    ground_link = np.array(ground_link, dtype=int)
    
    for i in range(len(bi_links)):
        if bi_links[i,0] == actuator:
            if jj[bi_links[i,1],bi_links[i,1]]!= 2:
                actuated_inx.append(i)
                actuated_link = np.array([actuator,bi_links[i,1]], dtype=int)
        elif bi_links[i,1] == actuator:
            if jj[bi_links[i,0],bi_links[i,0]]!= 2:
                actuated_inx.append(i)
                actuated_link = np.array([actuator,bi_links[i,0]], dtype=int)
                
    for i in range(len(tri_links)):
        if tri_links[i,0] == actuator:
            actuated_inx.append(i)
            actuated_link = np.array([actuator,tri_links[i,1],tri_links[i,2]], dtype=int)
        elif tri_links[i,2] == actuator:
            actuated_inx.append(i)
            actuated_link = np.array([actuator,tri_links[i,1],tri_links[i,0]], dtype=int)
    
    for i in range(len(bi_links)):
        if all(bi_links[i]==ground_link):
            ground_inx.append(i)
    ground_inx = np.array(ground_inx, dtype=int)
    
    
    return bi_links, tri_links, ground_link, ground_inx, actuated_link, actuated_inx

# Find LJ matrix
def find_link_joint_table(jj, bi_links, tri_links):
    lj = np.zeros((len(bi_links)+len(tri_links),len(jj)))
    
    for i in range(len(bi_links)):
        for j in range(len(bi_links[i])):
            lj[i,bi_links[i,j]] = 1

    for i in range(len(bi_links),len(bi_links)+len(tri_links)):
        for j in range(len(tri_links[i-len(bi_links)])):
            lj[i,tri_links[i-len(bi_links),j]] = 1
        
    
    lj = lj.astype(int) 
    return lj

def generate_unknown_table (jj):
    columns = len(jj)
    row = 9
    unknown_params = sp.symbols(f'x0:{columns*row}')
    
    un_table = np.empty((columns, row), dtype=object)
    for i in range(columns):
        for j in range(row):
            un_table[i,j] = unknown_params[i*9+j]
    
    return un_table

def compute_angle(jc_1, jc_2, ja_1, ja_2):
    no_angle = np.array([0,0,0])
    angle_1 = 0
    angle_2 = 0
    angle_1_2 = 0
    
    
    # the angle of the 1st joint and the link
    if (all(ja_1 == no_angle)):
        angle_1 = 100   #100 here means the angle doesn't exist!!
        angle_1_2 = 100 # if angle_1 doesn't exist, the angle between two joints doesn't exist
    else:
        angle_1 = np.dot(ja_1, jc_2 - jc_1)/np.linalg.norm(jc_2 - jc_1)
        angle_1 = regular_angle(angle_1)
        angle_1 = np.arccos(angle_1)
        
    # the angle of the 2nd joint and the link
    if (all(ja_2 == no_angle)):
        angle_2 = 100   #100 here means the angle doesn't exist!!
        angle_1_2 = 100 # if angle_2 doesn't exist, the angle between two joints doesn't exist
    else:
        angle_2 = np.dot(ja_2, jc_1 - jc_2)/np.linalg.norm(jc_2 - jc_1)
        #print("angle_2",angle_2)
        angle_2 = regular_angle(angle_2)
        angle_2 = np.arccos(angle_2)
        #print("angle_2",angle_2)
        
    # the angle between two joint
    if angle_1_2 != 100:
        angle_1_2 = np.dot(ja_1, ja_2)
        angle_1_2 = regular_angle(angle_1_2)
        angle_1_2 = np.arccos(angle_1_2)
         
    if angle_1 == 100:
        angle_1 = np.pi/2
    if angle_2 == 100:
        angle_2 = np.pi/2
    if angle_1_2 == 100:
        angle_1_2 = np.pi/2
    
    return(np.array([angle_1, angle_2, angle_1_2]))
    
def compute_bi_link_length_angle (lj, jt, bi_links, jc, ja):
    u_trigger = np.zeros(len(jt))
    bi_link_len = np.zeros(len(bi_links))
    bi_link_ang = np.zeros((len(bi_links),3))
    
    for i in range(len(bi_links)):
        #print("The link:",i+1)
        inx0, inx1 = bi_links[i,0], bi_links[i,1]
        bi_link_len[i] = np.linalg.norm(jc[inx0]-jc[inx1])
        
        if jt[inx0] != 2 and jt[inx1] !=2: # both joints are not u joint
            bi_link_ang[i] = compute_angle(jc[inx0], jc[inx1], ja[inx0,0], ja[inx1,0])
            #print(compute_angle(jc[inx0], jc[inx1], ja[inx0,0], ja[inx1,0]))
#             print('no u')
#             print(i)
#             print(bi_link_ang)
            
        elif jt[inx0] == 2 and jt[inx1] !=2:# j1 is u joint and j2 is not
            if u_trigger[inx0] == 0: # u joint's first angle always belongs to the first link when calculate in order
                bi_link_ang[i] = compute_angle(jc[inx0], jc[inx1], ja[inx0,0], ja[inx1,0])
                u_trigger[inx0] = 1 # mark this u joint's first angle is calculted          
            else:
                #u_trigger[inx0] == 1 # use u joint's second angle for the second link
                bi_link_ang[i] = compute_angle(jc[inx0], jc[inx1], ja[inx0,1], ja[inx1,0])
#             print('j1 u, j2 not')
#             print(i)
#             print(bi_link_ang)
            
        elif jt[inx0] != 2 and jt[inx1] ==2:# j1 is not u joint and j2 is u joint
            if u_trigger[inx1] == 0: # u joint's first angle always belongs to the first link when calculate in order
                bi_link_ang[i] = compute_angle(jc[inx0], jc[inx1], ja[inx0,0], ja[inx1,0])
                u_trigger[inx1] = 1 # mark this u joint's first angle is calculted
            else:
                #u_trigger[inx1] == 1 # use u joint's second angle for the second link
                bi_link_ang[i] = compute_angle(jc[inx0], jc[inx1], ja[inx0,0], ja[inx1,1])
#             print('j1 not j2 u')
#             print(i)
#             print(bi_link_ang)
            
        else: #both j1 and j2 are u joints
            if u_trigger[inx0] == 0 and u_trigger[inx1] == 0:
                bi_link_ang[i] = compute_angle(jc[inx0], jc[inx1], ja[inx0,0], ja[inx1,0])
                u_trigger[inx0] = 1
                u_trigger[inx1] = 1
            elif u_trigger[inx0] == 1 and u_trigger[inx1] == 0:
                bi_link_ang[i] = compute_angle(jc[inx0], jc[inx1], ja[inx0,1], ja[inx1,0])
                u_trigger[inx1] = 1
            elif u_trigger[inx0] == 0 and u_trigger[inx1] == 1:
                bi_link_ang[i] = compute_angle(jc[inx0], jc[inx1], ja[inx0,0], ja[inx1,1])
                u_trigger[inx0] = 1
            else:
                bi_link_ang[i] = compute_angle(jc[inx0], jc[inx1], ja[inx0,1], ja[inx1,1])
#             print('both u')
#             print(i)
#             print(bi_link_ang)            
#     print(u_trigger)
#     print(bi_link_len)
#     print(bi_link_ang)
        
    return bi_link_len, bi_link_ang

def build_bi_link_constraints (un, jt, bi_link, bi_link_len, bi_link_ang, ground_link):
    bi_links = copy.deepcopy(bi_link)
    u_trigger = np.zeros(len(jt))
    F = []
    
#     for i in range(len(bi_links)): # Ground link only uses the first joint of u
#         if (bi_links[i] == ground_link).all():
#             bi_links = np.delete(bi_link,i,0)
#             if ground_link[0] == 2:
#                 u_trigger[ground_link[0]] = 1
#             if ground_link[1] == 2:
#                 u_trigger[ground_link[1]] = 1
    
    for i in range(len(bi_links)):
        #print("link:",i+1)
        inx0, inx1 = bi_links[i,0], bi_links[i,1]
        #f = (un[inx0,0:3]-un[inx1,0:3]).dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]**2 # length constraint
        f = sp.sqrt((un[inx0,0:3]-un[inx1,0:3]).dot(un[inx0,0:3]-un[inx1,0:3])) - bi_link_len[i] # length constraint
        F.append(f)
        #print(f)
        
        if jt[inx0] != 2 and jt[inx1] !=2: # both joints are not u joint
            f = un[inx0,3:6].dot(un[inx1,0:3]-un[inx0,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,0])
            F.append(f)
            #print(f)
            f = un[inx1,3:6].dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,1])
            F.append(f)
            #print(f)
            f = un[inx0,3:6].dot(un[inx1,3:6]) - np.cos(bi_link_ang[i,2])
            F.append(f)
            #print(f)
            
        elif jt[inx0] == 2 and jt[inx1] !=2:# j1 is u joint and j2 is not
            if u_trigger[inx0] == 0: # u joint's first angle always belongs to the first link when calculate in order
                f = un[inx0,3:6].dot(un[inx1,0:3]-un[inx0,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,0])
                F.append(f)
                #print(f)
                f = un[inx1,3:6].dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,1])
                F.append(f)
                #print(f)
                f = un[inx0,3:6].dot(un[inx1,3:6]) - np.cos(bi_link_ang[i,2])
                F.append(f)
                #print(f)
                u_trigger[inx0] = 1 # mark this u joint's first angle is calculted          
            else:
                f = un[inx0,6:9].dot(un[inx1,0:3]-un[inx0,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,0])
                F.append(f)
                #print(f)
                f = un[inx1,3:6].dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,1])
                F.append(f)
                #print(f)
                f = un[inx0,6:9].dot(un[inx1,3:6]) - np.cos(bi_link_ang[i,2])
                F.append(f)
                #print(f)
            
        elif jt[inx0] != 2 and jt[inx1] ==2:# j1 is not u joint and j2 is u joint
            if u_trigger[inx1] == 0: # u joint's first angle always belongs to the first link when calculate in order
                f = un[inx0,3:6].dot(un[inx1,0:3]-un[inx0,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,0])
                F.append(f)
                #print(f)
                f = un[inx1,3:6].dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,1])
                F.append(f)
                #print(f)
                f = un[inx0,3:6].dot(un[inx1,3:6]) - np.cos(bi_link_ang[i,2])
                F.append(f)
                #print(f)
                u_trigger[inx1] = 1 # mark this u joint's first angle is calculted
            else:
                f = un[inx0,3:6].dot(un[inx1,0:3]-un[inx0,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,0])
                F.append(f)
                #print(f)
                f = un[inx1,6:9].dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,1])
                F.append(f)
                #print(f)
                f = un[inx0,3:6].dot(un[inx1,6:9]) - np.cos(bi_link_ang[i,2])
                F.append(f)
                #print(f)

            
        else: #both j1 and j2 are u joints
            if u_trigger[inx0] == 0 and u_trigger[inx1] == 0:
                f = un[inx0,3:6].dot(un[inx1,0:3]-un[inx0,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,0])
                F.append(f)
                f = un[inx1,3:6].dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,1])
                F.append(f)
                f = un[inx0,3:6].dot(un[inx1,3:6]) - np.cos(bi_link_ang[i,2])
                F.append(f)
                u_trigger[inx0] = 1
                u_trigger[inx1] = 1
            elif u_trigger[inx0] == 1 and u_trigger[inx1] == 0:
                f = un[inx0,6:9].dot(un[inx1,0:3]-un[inx0,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,0])
                F.append(f)
                f = un[inx1,3:6].dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,1])
                F.append(f)
                f = un[inx0,6:9].dot(un[inx1,3:6]) - np.cos(bi_link_ang[i,2])
                F.append(f)                
                u_trigger[inx1] = 1
            elif u_trigger[inx0] == 0 and u_trigger[inx1] == 1:
                f = un[inx0,3:6].dot(un[inx1,0:3]-un[inx0,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,0])
                F.append(f)
                f = un[inx1,6:9].dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,1])
                F.append(f)
                f = un[inx0,3:6].dot(un[inx1,6:9]) - np.cos(bi_link_ang[i,2])
                F.append(f)
                u_trigger[inx0] = 1
            else:
                f = un[inx0,6:9].dot(un[inx1,0:3]-un[inx0,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,0])
                F.append(f)
                f = un[inx1,6:9].dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,1])
                F.append(f)
                f = un[inx0,6:9].dot(un[inx1,6:9]) - np.cos(bi_link_ang[i,2])
                F.append(f)
        #print(F)
    return F

def compute_tri_link_length_angle (lj, jt, tri_links, jc, ja):
    tri_link_len = np.zeros(len(tri_links))
    tri_link_ang = np.zeros((len(tri_links),3))
    for i in range(len(tri_links)):
        inx0, inx1, inx2 = tri_links[i,0],tri_links[i,1],tri_links[i,2] #(l_012, 1 is always w joint)
        tri_link_len[i] = np.linalg.norm(jc[inx0]-jc[inx1]) # link length of l01 # this means joint[inx2] has to be the main C/P joint
        
        if jt[inx0] == 2: # if j0 is u joint
            angle_1 = np.dot(ja[inx0,1], jc[inx1] - jc[inx0])/np.linalg.norm(jc[inx1]-jc[inx0])            
            angle_2 = np.dot(ja[inx0,1], ja[inx2,0])
        else:
            angle_1 = np.dot(ja[inx0,0], jc[inx1] - jc[inx0])/np.linalg.norm(jc[inx1]-jc[inx0])
            angle_2 = np.dot(ja[inx0,0], ja[inx2,0])
        
        angle_3 = np.dot(jc[inx0] - jc[inx1], ja[inx2,0])/np.linalg.norm(jc[inx0] - jc[inx1])
        angle_3 = regular_angle(angle_3)
        angle_3 = np.arccos(angle_3)         
        
        angle_1 = regular_angle(angle_1)
        angle_1 = np.arccos(angle_1) 
        angle_2 = regular_angle(angle_2)
        angle_2 = np.arccos(angle_2) 
        
        tri_link_ang[i,0] = angle_1
        tri_link_ang[i,1] = angle_2
        tri_link_ang[i,2] = angle_3
        
    return tri_link_len, tri_link_ang

def build_tri_link_constraints (un, jt, tri_links, tri_link_len, tri_link_ang):
    F=[]
    p_inx = 0
    for i in range(len(tri_links)):
        inx0, inx1, inx2 = tri_links[i,0],tri_links[i,1],tri_links[i,2] #(l_012, 0 is normal joint,1 is always w joint, 2 is c/p)
        
        #f = (un[inx0,0:3]-un[inx1,0:3]).dot(un[inx0,0:3]-un[inx1,0:3]) - tri_link_len[i]**2 # length constraint
        f = sp.sqrt((un[inx0,0:3]-un[inx1,0:3]).dot(un[inx0,0:3]-un[inx1,0:3])) - tri_link_len[i] # length constraint
        F.append(f)
        
        if jt[inx0] == 2: # if j0 is u joint
            f = un[inx0,6:9].dot(un[inx1,0:3]-un[inx0,0:3]) - tri_link_len[i]*np.cos(tri_link_ang[i,0]) # angle_1
            F.append(f)
            f = un[inx0,6:9].dot(un[inx2,3:6]) - np.cos(tri_link_ang[i,1]) # angle_2
            F.append(f)

        else:
            f = un[inx0,3:6].dot(un[inx1,0:3]-un[inx0,0:3]) - tri_link_len[i]*np.cos(tri_link_ang[i,0]) # angle_1
            F.append(f)
            f = un[inx0,3:6].dot(un[inx2,3:6]) - np.cos(tri_link_ang[i,1]) # angle_2
            F.append(f)
        
        f = (un[inx0,0:3]-un[inx1,0:3]).dot(un[inx2,3:6]) - tri_link_len[i]*np.cos(tri_link_ang[i,2])
        F.append(f)
        
        # Parallel Constraints
        u12 = (un[inx2,0:3] - un[inx1,0:3])
        u2 = un[inx2,3:6]
        #f = u12[0]*u2[0] + u12[1]*u2[1] + u12[2]*u2[2] - sp.sqrt(u12.dot(u12))*sp.sqrt(u2.dot(u2))
        #f = u12.dot(u12)*u2.dot(u2) - u2.dot(u12)*u2.dot(u12)
        f = (u12[2]*u2[0] - u12[0]*u2[2])*(u12[2]*u2[0] - u12[0]*u2[2]) + (u12[1]*u2[2] - u12[2]*u2[1])*(u12[1]*u2[2] - u12[2]*u2[1]) + (u12[0]*u2[1] - u12[1]*u2[0])*(u12[0]*u2[1] - u12[1]*u2[0])
        F.append(f)
        f = (u12[2]*u2[0] - u12[0]*u2[2])*(u12[2]*u2[0] - u12[0]*u2[2]) - (u12[1]*u2[2] - u12[2]*u2[1])*(u12[1]*u2[2] - u12[2]*u2[1]) + (u12[0]*u2[1] - u12[1]*u2[0])*(u12[0]*u2[1] - u12[1]*u2[0])
        F.append(f)
            
#         #print("u2:",u2)
#         #print("u12:",u12)
#         cross_1 = u12[1]*u2[2] - u12[2]*u2[1]
#         cross_2 = u12[2]*u2[0] - u12[0]*u2[2]
#         cross_3 = u12[0]*u2[1] - u12[1]*u2[0]
 
    return F

def compute_type_angle(jt, ja):
    type_angle = np.zeros(len(jt))
    for i in range(len(jt)):
        if jt[i] == 2: # if the joint type is U, calculate the joint angle between two axis
            angle = np.dot(ja[i,0], ja[i,1])
            angle = regular_angle(angle)
            angle = np.arccos(angle)
            type_angle[i] = angle
        else:
            type_angle[i] = 100 # 100 here means the angle doesn't exist
            
    return type_angle

def build_type_angle_constraints (un, type_angle):
    F = []
    for i in range(len(type_angle)):
        if type_angle[i] != 100:
            f = un[i,3:6].dot(un[i,6:9]) - np.cos(type_angle[i])
            F.append(f)
    return F

def build_unit_angle_constraints (un, jt):
    F = []
    for i in range(len(jt)):
        if jt[i] == 2: # U joint
            f = (un[i,3:6].dot(un[i,3:6])) - 1
            F.append(f)
            f =(un[i,6:9].dot(un[i,6:9])) - 1
            F.append(f)  
        elif jt[i] == 6 or jt[i] == 4: # W joint
            continue
        else:
            f = (un[i,3:6].dot(un[i,3:6])) - 1
            F.append(f)     
    return F

def compute_p_angle(jt, bi_links, tri_links, jc, ja):
    
    Angle_p = []
    for i in range(len(tri_links)):
        inx0, inx1, inx2 = tri_links[i,0],tri_links[i,1],tri_links[i,2] #(l_012, 1 is always w joint)

#         if jt[inx0] == 1: # if j0 is p joint, p belongs to this link, and there must exist another tri_link               
#             for j in range(len(tri_links)):
#                 if j!= i:
#                     if inx0 == tri_links[j,2]:
#                         inx_t0, inx_t1, inx_t2 = tri_links[j,0], tri_links[j,1], tri_links[j,2]
#                         angle_p = np.dot(jc[inx_t0] - jc[inx_t1], jc[inx0] - jc[inx1])/(np.linalg.norm(jc[inx0]-jc[inx1])*np.linalg.norm(jc[inx_t0] - jc[inx_t1]))
#                         angle_p = regular_angle(angle_p)
#                         angle_p = np.arccos(angle_p)
#                         Angle_p.append(angle_p)                          
                              
        
        if jt[inx2] == 1: # if j2 is p joint
            #print("inx2 is p")
            for j in range(len(bi_links)):
                if inx2 == bi_links[j,0] or inx2 == bi_links[j,1]:
                    #print("bilink,",bi_links[j])
                    inx_b1, inx_b2 = bi_links[j,0], bi_links[j,1]
                    angle_p = np.dot(jc[inx_b2] - jc[inx_b1], jc[inx0] - jc[inx1])/(np.linalg.norm(jc[inx0]-jc[inx1])*np.linalg.norm(jc[inx_b2] - jc[inx_b1]))
                    angle_p = regular_angle(angle_p)
                    angle_p = np.arccos(angle_p)
                    Angle_p.append(angle_p)
            for j in range(len(tri_links)):
                if j!= i:
                    if inx2 == tri_links[j,0]:
                        #print("tri_links,",tri_links[j])
                        inx_t0, inx_t1, inx_t2 = tri_links[j,0], tri_links[j,1], tri_links[j,2]
                        angle_p = np.dot(jc[inx_t0] - jc[inx_t1], jc[inx0] - jc[inx1])/(np.linalg.norm(jc[inx0]-jc[inx1])*np.linalg.norm(jc[inx_t0] - jc[inx_t1]))
                        angle_p = regular_angle(angle_p)
                        angle_p = np.arccos(angle_p)
                        Angle_p.append(angle_p)                        
    Angle_p = np.array(Angle_p)
        
        
    return Angle_p

"warning: where does joint P belongs to, and how to choose the constraint hasn't been completely solved. The code works for"
"1-dof four-bar, but could have issue on other cases."
"Read the note"
def build_p_constraints(jt, bi_links, bi_link_len, tri_links, tri_link_len, jc, ja, angle_p, un):
    F = []
    inx_p = 0
    
    for i in range(len(tri_links)):
        inx0, inx1, inx2 = tri_links[i,0],tri_links[i,1],tri_links[i,2] #(l_012, 1 is always w joint)
        
#         if jt[inx0] == 1: # if j0 is p joint, p belongs to this link, and there must exist another tri_link                       
#             for j in range(len(tri_links)):
#                 if j!= i:
#                     if inx0 == tri_links[j,2]:
#                         inx_t0, inx_t1, inx_t2 = tri_links[j,0], tri_links[j,1], tri_links[j,2]
#                         f = (un[inx_t0,0:3]-un[inx_t1,0:3]).dot(un[inx0,0:3]-un[inx1,0:3]) - tri_link_len[i]*tri_link_len[j]*np.cos(angle_p[inx_p])
#                         F.append(f)
#                         inx_p += 1        
        
        if jt[inx2] == 1: # if j2 is p joint, then constraint p 
            for j in range(len(bi_links)):
                if inx2 == bi_links[j,0] or inx2 == bi_links[j,1]:
                    inx_b1, inx_b2 = bi_links[j,0], bi_links[j,1]
                    f = (un[inx_b2,0:3]-un[inx_b1,0:3]).dot(un[inx0,0:3]-un[inx1,0:3]) - tri_link_len[i]*bi_link_len[j]*np.cos(angle_p[inx_p])
                    F.append(f)
                    inx_p += 1
            for j in range(len(tri_links)):
                if j!= i:
                    if inx2 == tri_links[j,0]:
                        inx_t0, inx_t1, inx_t2 = tri_links[j,0], tri_links[j,1], tri_links[j,2]
                        f = (un[inx_t0,0:3]-un[inx_t1,0:3]).dot(un[inx0,0:3]-un[inx1,0:3]) - tri_link_len[i]*tri_link_len[j]*np.cos(angle_p[inx_p])
                        F.append(f)
                        inx_p += 1 
                    
    return F

def assign_values(symbols, values):
    """
    Assign values to symbolic parameters.

    Args:
        symbols (list or tuple): A list or tuple of symbolic parameters (e.g., [x0, x1, x2]).
        values (list or tuple): A list or tuple of values to substitute into the symbols.

    Returns:
        dict: A dictionary mapping symbols to their values.
    """
#     if len(symbols) != len(values):
#         raise ValueError("The number of symbols and values must match.")
    
    # Create a dictionary mapping symbols to values
    substitutions = {symbol: value for symbol, value in zip(symbols, values)}
    #print(type(substitutions[symbols[0]]))
    
    return substitutions

def update_actuated_link (actuate_link, actuator, actuate_inx, ja, jc, jt, bi_link_len, tri_link_len, phi):
    
    # using quaternion to calculate new position
    #if actuate_link[0] == actuator: 
    actuate_vector = (jc[actuate_link[1],:] - jc[actuate_link[0],:])/np.linalg.norm(jc[actuate_link[1],:] - jc[actuate_link[0],:])
    actuate_vector =  np.append(0, actuate_vector)
    #else:
    #    print('Actuator format wrong. Please let the first element of the actuate_link be where the actuator located on')
    
    if len(actuate_link) == 2:
        link_len = bi_link_len[actuate_inx[0]]
    elif len(actuate_link) == 3:
        link_len = tri_link_len[actuate_inx[0]]
    
    theta = np.pi*phi/360; #thera is half of phi  
    actuator_axis = ja[actuator,0]
    rotation_quater = np.array([np.cos(theta), np.sin(theta)*actuator_axis[0], np.sin(theta)*actuator_axis[1], np.sin(theta)*actuator_axis[2]]) #the quaternion of rotation
    new_actuated_vec = QuaterTimes(QuaterTimes(rotation_quater,actuate_vector),QuaterConj(rotation_quater))
    new_pos = new_actuated_vec[1:4]*link_len + jc[actuator] # new position of the joint connected to the actuated joint
    new_axis = np.array([0,0,0,0], dtype=np.float64)
    
    #if jt[actuate_link[1]] ==  2: # if it's u joint
        #new_axis = rotate_vector(ja[actuate_link[1],0], actuator_axis, theta*2) # new vector of the joint connected to the actuated joint if it has
    new_axis = QuaterTimes(QuaterTimes(rotation_quater,np.append(0,ja[actuate_link[1],0])),QuaterConj(rotation_quater))
        #????????updated_actuated_link = np.vstack((new_pos,new_axis[1:]))
        #除了S joint以外，R, P, U, C的new_axis都应该被更新
    
    if len(actuate_link) == 2:
        updated_actuated_link = np.vstack((new_pos,new_axis[1:]))
    elif len(actuate_link) == 3:
        new_axis_2 = QuaterTimes(QuaterTimes(rotation_quater,np.append(0,ja[actuate_link[2],0])),QuaterConj(rotation_quater))
        updated_actuated_link = np.vstack((new_pos,new_axis[1:],new_axis_2[1:]))
    else:
        print("something wrong with the update_actuated_link")
    
    #print(updated_actuated_link)
    return updated_actuated_link
    
def assign_value_for_knowns(ground_link, ground_inx, actuate_link, jt, ori_un, ja, jc, bi_links):#, updated_actuated_link):
    un = copy.deepcopy(ori_un)
    # assign 0 to non-rotational axis for non-u joints
    for i in range(len(jt)):
        if jt[i] != 2:
            un[i,6:9] = np.array([0,0,0])
        if jt[i] == 6 or jt[i] == 4:
            un[i,3:6] = np.array([0,0,0])
        
    # assign original value to the fixed joints
    for i in ground_link:
        un[i,0:3] = jc[i]
        if jt[i] == 2: # judge if use the first or the second joint of u
            for j in range(len(bi_links)):
                if i == bi_links[j,0] or i == bi_links[j,1]:
                    if ground_inx <= j: # U还没有被用过，ground用u1，u1是known 
                        un[i,3:6] = ja[i,0]
                        #print(i, 'th ground u1, bi_link', j)
                    else: # U被用过了，ground用u2，u2是known
                        un[i,6:9] = ja[i,1]
                        #print(i, 'th ground u1, bi_link', j)
                    break
            #un[i,3:6] = ja[i,0]
        else:
            un[i,3:6] = ja[i,0]
            un[i,6:9] = ja[i,1]
        
    # assign updated actuated link location and axis
#     un[actuate_link[1], 0:3] = updated_actuated_link[0]
#     un[actuate_link[1], 3:6] = updated_actuated_link[1]
                    
    return un

def remove_unvalid_constraint (initial_constraint):
    removal = []
    for i in range(len(initial_constraint)):
        #print("----------------------")
        #print(i)
        #print(type(initial_constraint[i]))
        #print(initial_constraint[i])
        if isinstance(initial_constraint[i],np.float64) or isinstance(initial_constraint[i],int) or isinstance(initial_constraint[i],float):
            if abs(initial_constraint[i]) < 1e-4:
                removal.append(i)
                #print("remove1:",initial_constraint[i])
            else:
                #print(abs(initial_constraint[i]))
                print('The', i, "th constraint is calculatable but not zero, the constraint equations have error 1")
        elif isinstance(initial_constraint[i],sp.core.numbers.Zero):# or isinstance(initial_constraint[i],sp.core.numbers.NaN):
            #print("remove2:",initial_constraint[i])
            removal.append(i)
        elif isinstance(initial_constraint[i],sp.core.numbers.Float):
            if abs(initial_constraint[i]) < 1e-4:
                removal.append(i)
                #print("remove3:",initial_constraint[i])
            else:
                print('The', i, "th constraint is calculatable but not zero, the constraint equations have error 2")
#         elif isinstance(initial_constraint[i], sp.core.add.Add) == False and isinstance(initial_constraint[i], sp.core.mul.Mul) == False:
#             removal.append(i)  
#             print("remove2:",initial_constraint[i])
    final_constraint = np.delete(initial_constraint,removal)
    return final_constraint

def find_real_unknowns_and_initial_guess(ground_link, ground_inx, actuate_link, jt, un, ja, jc, bi_links, updated_actuated_link):
    real_un = []
    initial_guess = []
    #p_inx = 0
    
    #c_inx = None
    for i in range(len(jt)):
        if (i == ground_link[0]) or (i == ground_link[1]): # skip ground link except the axis of u joint
            if jt[i] == 2: 
                #real_un.append(un[i,6:9])
                #initial_guess.append(ja[i,1])
                for j in range(len(bi_links)):
                    if bi_links[j,0] == i or bi_links[j,1] == i: 
                        if ground_inx > j: #u已被用过，地杆占角u_2， u_1未知
                            real_un.append(un[i,3:6])
                            initial_guess.append(ja[i,0])
                        else: # u在地杆内，或地杆后，地杆占角u_1， u_2未知
                            real_un.append(un[i,6:9]) #u2未知
                            initial_guess.append(ja[i,1])
                        break
#             elif jt[i] == 3 or jt[i] == 1: #c or p
#                 c_inx = 100 # find the inx of C joint's axis's y value

            else:
                continue
            #print(real_un)
                
        elif i == actuate_link[1]: #skip the actuated link except the axis of u joint
            if jt[i] == 2:
                real_un.append(un[i,6:9])
                initial_guess.append(ja[i,1]) # update intial guess with the updated actuated joint it's a
            else:
                continue
        elif (len(actuate_link) == 3) and (i == actuate_link[2]): # if the actuated link is a tri_link, then the unkown 
            # will be the position of c or p joint while the c/p joint's axis will be known
            "it's risky here because if len(actuate_link)=2, then actuate_link[2] won't exist. But it seems python will calculate"
            "the part before and first..."
            #if i == actuate_link[2]:
            #print("skip:",i)
            real_un.append(un[i,0:3])
            initial_guess.append(jc[i])
        
        elif jt[i] == 2: #u
            real_un.append(un[i])
            initial_guess.append(jc[i])
            initial_guess.append(ja[i,0])
            initial_guess.append(ja[i,1])
        elif jt[i] == 6 or jt[i] == 4: # w and s
            real_un.append(un[i,0:3])
            initial_guess.append(jc[i])
        else:
            #print("the len is:",real_un[0])
            #print("the len is:",len(real_un[0]))
#             if jt[i] == 3 or jt[i] == 1 : #c or p
#                 p_inx += 1
            real_un.append(un[i,0:6])
            initial_guess.append(jc[i])
            initial_guess.append(ja[i,0])
            
#     for j in range(p_inx):
#         real_un.append(np.array([un[len(jt),j]]))
#         initial_guess.append(np.array([1]))
    
    real_un = np.concatenate(real_un)
    initial_guess = np.concatenate(initial_guess)
    return real_un, initial_guess

# def solve_constraint_equations(constraint_equations_sym, unknown_params, guess):
    
#     def constraint_eqs(upara):
#         #print(np.array([eq(*upara) for eq in constraint_equations]))
#         return np.array([eq(*upara) for eq in constraint_equations])
    
#     constraint_equations = [sp.lambdify(unknown_params, eq) for eq in constraint_equations_sym]
#     #print(constraint_eqs)

#     result = root(constraint_eqs, guess, tol = 1e-5)#,method='krylov')
#     #print(result)
#     #print("-----------------")

    
#     if result.success:
#         solution = result.x
# #         print(result.fun)
# #         print("-----------------")
# #         print(solution)
#         return solution, True
#     else:
#         mes = result.message;
#         #print(result.fun)
#         #print(mes)
        
        # return None, False    


def find_actuator_inx (un, actuate_link):
    if len(actuate_link) == 2:
        un_actuator_inx = np.hstack((un[actuate_link[1], 0:3], un[actuate_link[1], 3:6]))
    elif len(actuate_link) == 3:
        un_actuator_inx = np.hstack((un[actuate_link[1], 0:3], un[actuate_link[1], 3:6], un[actuate_link[2], 3:6]))
    else:
        print("something wrong with the actuator")
        
    return un_actuator_inx

def assign_value_for_actuator(final_constraint, un_actuator_inx, updated_actuated_link):
    updated_constraint = []
    values = assign_values(un_actuator_inx,updated_actuated_link.flatten())
    #print(values)
    for i in range(len(final_constraint)):
        updated_constraint.append(final_constraint[i].xreplace(values))

    updated_constraint = np.array(updated_constraint)
    updated_constraint = remove_unvalid_constraint(updated_constraint)
    return updated_constraint

#################################### Plot and Normalization Functions #####################################################
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import directed_hausdorff
from scipy.interpolate import splprep, splev

def plotLine3D(P1,P2, color='gray',linewidth=2,alpha=1):
    ax.plot3D([P1[0],P2[0]], [P1[1],P2[1]], [P1[2],P2[2]],linewidth=linewidth, color=color, alpha=alpha)

def plotMec(jt, jc, ja, bi_links, tri_links, ground_link, coupler, scale=0.5):
    for i in range(len(jt)):
        if jt[i] == 0 or jt[i] == 1: # R, P, C joint
            ax.plot3D([jc[i,0], jc[i,0]+scale*ja[i,0,0]], [jc[i,1], jc[i,1]+scale*ja[i,0,1]], [jc[i,2], jc[i,2]+scale*ja[i,0,2]],linewidth=2, color='blue') 
            ax.plot3D([jc[i,0], jc[i,0]-scale*ja[i,0,0]], [jc[i,1], jc[i,1]-scale*ja[i,0,1]], [jc[i,2], jc[i,2]-scale*ja[i,0,2]],linewidth=2, color='blue')  

            ax.plot3D([jc[i,0], jc[i,0]+0.3*scale*ja[i,0,0]], [jc[i,1], jc[i,1]+0.3*scale*ja[i,0,1]], [jc[i,2], jc[i,2]+0.3*scale*ja[i,0,2]],linewidth=6, color='blue') 
            ax.plot3D([jc[i,0], jc[i,0]-0.3*scale*ja[i,0,0]], [jc[i,1], jc[i,1]-0.3*scale*ja[i,0,1]], [jc[i,2], jc[i,2]-0.3*scale*ja[i,0,2]],linewidth=6, color='blue')  

        elif jt[i] == 3: # C joint
            ax.plot3D([jc[i,0], jc[i,0]+scale*ja[i,0,0]], [jc[i,1], jc[i,1]+scale*ja[i,0,1]], [jc[i,2], jc[i,2]+scale*ja[i,0,2]],linewidth=2, color='magenta') 
            ax.plot3D([jc[i,0], jc[i,0]-scale*ja[i,0,0]], [jc[i,1], jc[i,1]-scale*ja[i,0,1]], [jc[i,2], jc[i,2]-scale*ja[i,0,2]],linewidth=2, color='magenta')  

            ax.plot3D([jc[i,0], jc[i,0]+0.3*scale*ja[i,0,0]], [jc[i,1], jc[i,1]+0.3*scale*ja[i,0,1]], [jc[i,2], jc[i,2]+0.3*scale*ja[i,0,2]],linewidth=6, color='magenta') 
            ax.plot3D([jc[i,0], jc[i,0]-0.3*scale*ja[i,0,0]], [jc[i,1], jc[i,1]-0.3*scale*ja[i,0,1]], [jc[i,2], jc[i,2]-0.3*scale*ja[i,0,2]],linewidth=6, color='magenta')  

        elif jt[i] == 2: # U joint
            ax.plot3D([jc[i,0], jc[i,0]+scale*ja[i,0,0]], [jc[i,1], jc[i,1]+scale*ja[i,0,1]], [jc[i,2], jc[i,2]+scale*ja[i,0,2]],linewidth=2, color='red') 
            ax.plot3D([jc[i,0], jc[i,0]-scale*ja[i,0,0]], [jc[i,1], jc[i,1]-scale*ja[i,0,1]], [jc[i,2], jc[i,2]-scale*ja[i,0,2]],linewidth=2, color='red')  
            
            ax.plot3D([jc[i,0], jc[i,0]+scale*ja[i,1,0]], [jc[i,1], jc[i,1]+scale*ja[i,1,1]], [jc[i,2], jc[i,2]+scale*ja[i,1,2]],linewidth=2, color='orange') 
            ax.plot3D([jc[i,0], jc[i,0]-scale*ja[i,1,0]], [jc[i,1], jc[i,1]-scale*ja[i,1,1]], [jc[i,2], jc[i,2]-scale*ja[i,1,2]],linewidth=2, color='orange')  

            
            ax.plot3D([jc[i,0], jc[i,0]+0.3*scale*ja[i,0,0]], [jc[i,1], jc[i,1]+0.3*scale*ja[i,0,1]], [jc[i,2], jc[i,2]+0.3*scale*ja[i,0,2]],linewidth=6, color='red') 
            ax.plot3D([jc[i,0], jc[i,0]-0.3*scale*ja[i,0,0]], [jc[i,1], jc[i,1]-0.3*scale*ja[i,0,1]], [jc[i,2], jc[i,2]-0.3*scale*ja[i,0,2]],linewidth=6, color='red')  
            
            ax.plot3D([jc[i,0], jc[i,0]+0.3*scale*ja[i,1,0]], [jc[i,1], jc[i,1]+0.3*scale*ja[i,1,1]], [jc[i,2], jc[i,2]+0.3*scale*ja[i,1,2]],linewidth=6, color='orange') 
            ax.plot3D([jc[i,0], jc[i,0]-0.3*scale*ja[i,1,0]], [jc[i,1], jc[i,1]-0.3*scale*ja[i,1,1]], [jc[i,2], jc[i,2]-0.3*scale*ja[i,1,2]],linewidth=6, color='orange')  

        elif jt[i] == 4: # S joint
            ax.scatter3D(jc[i,0],jc[i,1],jc[i,2],marker='o',linewidths=3, color='green')
        
        elif jt[i] == 6: # W joint
            ax.scatter3D(jc[i,0],jc[i,1],jc[i,2],marker='o',linewidths=2, color='black')
    
    for i in range(len(bi_links)):
        inx0, inx1 = bi_links[i,0], bi_links[i,1]
        plotLine3D(jc[inx0],jc[inx1], color='green')

    for i in range(len(tri_links)):
        inx0, inx1, inx2 = tri_links[i,0], tri_links[i,1], tri_links[i,2]
        plotLine3D(jc[inx0],jc[inx1], color='cyan')
        plotLine3D(jc[inx1],jc[inx2], color='cyan')
    
    ground_inx0, ground_inx1 = ground_link[0], ground_link[1]
    plotLine3D(jc[ground_inx0],jc[ground_inx1], color='black',linewidth=3)
    ax.scatter3D(coupler[0],coupler[1],coupler[2],marker='o',linewidths=0.1, color='orange')
    
def plotPath(Pts, limit=5, color = 'gray',linestyle = 'line1', dot_size = 1):

        
    if linestyle == 'line1':
        xline=Pts[:,0]
        yline=Pts[:,1]
        zline=Pts[:,2]
        xline = np.append(xline,xline[0])
        yline = np.append(yline,yline[0])
        zline = np.append(zline,zline[0])
        ax.plot3D(xline, yline, zline, color)
    elif linestyle == 'point1':
        ax.scatter3D(Pts[:,0], Pts[:,1], Pts[:,2], c=color, s=dot_size)
        #ax.scatter3D(Pts[0,0], Pts[0,1], Pts[0,2], c='Green', s=dot_size*5)
        #ax.scatter3D(Pts[-1,0], Pts[-1,1], Pts[-1,2], c='Blue', s=dot_size*5)
    elif linestyle == 'point2':
        ax.scatter3D(Pts[0,:], Pts[1,:], Pts[2,:], c=color, s=dot_size)
        #ax.scatter3D(Pts[0,0], Pts[1,0], Pts[2,0], c='Green', s=dot_size*5)
        #ax.scatter3D(Pts[0,-1], Pts[1,-1], Pts[2,-1], c='Blue', s=dot_size*5)
    elif linestyle == 'line2':
        xline=Pts[0,:]
        yline=Pts[1,:]
        zline=Pts[2,:]
        xline = np.append(xline,xline[0])
        yline = np.append(yline,yline[0])
        zline = np.append(zline,zline[0])
        ax.plot3D(xline, yline, zline, color)
    else:
        print('wrong linesyle')
    
    #ax.auto_scale_xyz([-limit, limit], [-limit, limit], [-limit, limit])
    #plt.tight_layout()

def IfRemove(path, th=1):
    condition = False
    dist = np.sum(np.power(np.subtract(path[0,:],path[-1,:]),2))
    if dist > th:
        condition = True
    return condition

def RemoveOpen (path_set, mec_set):
    num = len(path_set)
    num2 = len(mec_set)
    
    if num != num2:
        print('Datasets have different length')
    
    inx = []
    for i in range(num):
        con = IfRemove(path_set[i])
        if con == True:
            inx.append(i)
            
    path_set = np.delete(path_set,inx, 0)
    mec_set = np.delete(mec_set,inx, 0)
    
    return path_set, mec_set

def scale_mec(omec):
    mec = copy.deepcopy(omec)
    mec[0] = omec[0]/3
    mec[1] = omec[1]/3
    mec[2] = omec[2]/3
    mec[5] = omec[5]/3

    return mec

#Evaluate the smooth factor of bspline
def Coordinate_meansquared(y1,y2):
    MSE = 0;
    sub = np.subtract(y1,y2);
    l = len(sub[0]);
    for i in range(l-1):
        diff = math.sqrt(sub[0,i]**2 + sub[1,i]**2 + sub[2,i]**2)
        MSE += diff
    MSE = MSE/l
    return MSE

def Path_Interpolate_Open(path, num_pts=100, smooth=0.01):
    px = path[:,0]
    py = path[:,1]
    pz = path[:,2]
    
    tck, u = splprep([px,py,pz],s=smooth) #splprep 给出插值结果y，splev用前者给出的新插值函数y，用输入的新节点x_new找到对应y值

    u_fine = np.linspace(0,1,num_pts)
    new_points = splev(u_fine,tck)
    
    return new_points

def Path_Interpolate_Closed(path, num_pts=100, smooth=0.01):
    px = path[:,0]
    py = path[:,1]
    pz = path[:,2]
    
    px = np.append(px,px[0]);
    py = np.append(py,py[0]);
    pz = np.append(pz,pz[0]);
    
    tck, u = splprep([px,py,pz],s=smooth) #splprep 给出插值结果y，splev用前者给出的新插值函数y，用输入的新节点x_new找到对应y值
    u_fine = np.linspace(0,1,num_pts)
    new_points = splev(u_fine,tck)
    
    return new_points


def reflect_data_3d(path):
    x = path[0,:]
    y = path[1,:]
    z = path[2,:]
    
    x_r = (x - np.mean(x))**3
    y_r = (y - np.mean(y))**3
    z_r = (z - np.mean(z))**3
    
    m11 = np.sum(x_r) / len(x_r) #x_e
    m12 = np.sum(y_r) / len(y_r) #y_e
    m13 = np.sum(z_r) / len(z_r) #z_e
    
    signm11 = np.sign(m11)
    signm12 = np.sign(m12)
    signm13 = np.sign(m13)
    
#     if np.abs(signm11) < 1e-5:
#         signm11 = 1
#     if np.abs(signm12) < 1e-5:
#         signm12 = 1
#     if np.abs(signm13) < 1e-5:
#         signm13 = 1
    
    reflectionMat = np.array(
    [[signm11, 0, 0],
     [0, signm12, 0],
     [0, 0, signm13]
    ])
    
    if np.abs(m11) < np.abs(m13) and np.abs(m13) < np.abs(m12):
        reflectionMat = np.matmul(np.array([[1,0,0],[0,0,1],[0,1,0]]), reflectionMat)
    elif np.abs(m12) < np.abs(m11) and np.abs(m11) < np.abs(m13):
        reflectionMat = np.matmul(np.array([[0,1,0],[1,0,0],[0,0,1]]), reflectionMat)
    elif np.abs(m12) < np.abs(m13) and np.abs(m13) < np.abs(m11):
        reflectionMat = np.matmul(np.array([[0,1,0],[0,0,1],[1,0,0]]), reflectionMat)
    elif np.abs(m13) < np.abs(m11) and np.abs(m11) < np.abs(m12):
        reflectionMat = np.matmul(np.array([[0,0,1],[1,0,0],[0,1,0]]), reflectionMat)
    elif np.abs(m13) < np.abs(m12) and np.abs(m12) < np.abs(m11):
        reflectionMat = np.matmul(np.array([[0,0,1],[0,1,0],[1,0,0]]), reflectionMat)
        
    reflected_curve = np.matmul(reflectionMat,path)
    
    return reflected_curve, reflectionMat

# Calculate the total arc length of a 3D curve
def calcArcLen(path):
    ArcLen=0
    path_t = path.T
    for j in range(0,len(path_t)-1):
        d=np.linalg.norm(path_t[j]-path_t[j+1])
        ArcLen+=d
    ArcLen += np.linalg.norm(path_t[-1]-path_t[0])
    return ArcLen

# Normolize path with reflection
def normalizePath(path):
    #Translation
    trans = np.mean(path,axis=1)
    trans_path = np.transpose(np.transpose(path) - trans) #tans_path size: 3*num
                            
    #Rotation
    cov_mat = np.cov(trans_path) # calculate covariance matrix
    eig_val,eig_vec = np.linalg.eigh(cov_mat)
    if np.linalg.det(eig_vec) <= 0:
        eig_vec = -eig_vec
    rot = np.linalg.inv(eig_vec)
    rot_path = np.matmul(rot,trans_path)
    #eigen = np.linalg.svd(trans_path.T)[2] #V^T, each row is a eigenvector
    #rot_path2 = np.matmul(np.linalg.inv(eigen.T),trans_path)
    
    #reflection
    reflec_path, reflect_matrix = reflect_data_3d(rot_path)
    re_rot = np.matmul(reflect_matrix,rot)
    
    #Scaling
    scale =calcArcLen(reflec_path)
    #scale = scale1/100
    #scale = 1
    
    norm_path=reflec_path/scale

    return norm_path, trans, re_rot, scale

# Normalize Mechanism
def normalizeMec(mec, trans, rot, scale):
    # Translation
    trans_mec = copy.deepcopy(mec)
    trans_mec[0,:] = mec[0,:]-trans
    trans_mec[1,:] = mec[1,:]-trans
    trans_mec[2,:] = mec[2,:]-trans
    trans_mec[5,:] = mec[5,:]-trans
    trans_mec = np.transpose(trans_mec)    
    
    #Rotation
    rot_mec = np.transpose(np.matmul(rot, trans_mec))
    
    #Scaling
    rot_mec[0,:] = rot_mec[0,:]/scale
    rot_mec[1,:] = rot_mec[1,:]/scale
    rot_mec[2,:] = rot_mec[2,:]/scale
    rot_mec[5,:] = rot_mec[5,:]/scale  
    
    return rot_mec

def normalizeJ0(jc0, ja0, trans, rot, scale):
    trans_jc0 = jc0-trans
    trans_jc0 = np.transpose(trans_jc0)
    rot_jc0 = np.transpose(np.matmul(rot, trans_jc0))
    new_jc0 = rot_jc0/scale
    
    new_ja0 = np.transpose(np.matmul(rot, ja0))
    return new_jc0, new_ja0
####################################  After Training Simulation Functions  #####################################################
def Open_curve_simulation (jj, jt, jc, ja, actuator, step, Root_accuracy = 1e-5):
    def fun(u, p):
        # u: 1D array for unknowns (x27..x32), p: 1D array for params (x9,x10,x11)
        vals = (*u, *p)
        return np.asarray(F_num(*vals), dtype=float).ravel()
    # Calculate the configuration
    Bi_links, Tri_links, Ground_link, Ground_inx, Actuate_link, Actuate_inx = find_link(jj, jt, actuator)
    LJ = find_link_joint_table(jj, Bi_links, Tri_links)
    UN = generate_unknown_table(jj)
    Bi_link_len, Bi_link_ang= compute_bi_link_length_angle (LJ, jt, Bi_links, jc, ja)
    Tri_link_len, Tri_link_ang = compute_tri_link_length_angle (LJ, jt, Tri_links, jc, ja)
    Spacial_p_ang = compute_p_angle(jt, Bi_links, Tri_links, jc, ja)
    Type_angle = compute_type_angle(jt, ja)

    # Find Constraint Eqs
    Updated_UN = assign_value_for_knowns(Ground_link, Ground_inx, Actuate_link, jt, UN, ja, jc, Bi_links) #update UN
    UN_actuator_inx = find_actuator_inx (UN, Actuate_link)

    Bi_constraint = build_bi_link_constraints (Updated_UN, jt, Bi_links, Bi_link_len, Bi_link_ang, Ground_link)
    Tri_constraint= build_tri_link_constraints (Updated_UN, jt, Tri_links, Tri_link_len, Tri_link_ang)
    Type_constraint = build_type_angle_constraints (Updated_UN, Type_angle)
    Unit_constraint = build_unit_angle_constraints (Updated_UN, jt)
    Special_p_constraint = build_p_constraints(jt, Bi_links, Bi_link_len, Tri_links, Tri_link_len, jc, ja, Spacial_p_ang, Updated_UN)
    Final_constraint = remove_unvalid_constraint(np.array(Bi_constraint + Tri_constraint+ Type_constraint + Unit_constraint + Special_p_constraint))

    "calculate forward from the initial point"
    phi = 0
    Updated_actuated_link = update_actuated_link(Actuate_link, actuator, Actuate_inx, ja, jc, jt, Bi_link_len, Tri_link_len,phi)
    #Find the initial unknowns and initial guess
    Unknown_paras, Initial_guess =find_real_unknowns_and_initial_guess(Ground_link, Ground_inx, Actuate_link, jt, UN, ja, jc, Bi_links, Updated_actuated_link)
    arglist = [*Unknown_paras, *UN_actuator_inx]
    F_filtered = sp.Matrix([f for f in Final_constraint if f.free_symbols & set(Unknown_paras)])
    F_num = sp.lambdify(arglist, F_filtered, modules = 'numpy') #change the constraint equations from symbolic to numerical
    F_wrapped = lambda u, p: fun(u, p)
    result_temp = root(F_wrapped, Initial_guess, args=(Updated_actuated_link.flatten()), tol = Root_accuracy)
    Initial_para, condition = result_temp.x, result_temp.success

    if condition == False:
        print("Doesn't have initial solution")
        return None
    
    "calculate forward from the initial point"
    phi = 0
    Step_para = np.zeros((step,len(Unknown_paras))) # store every position of point_4
    Step_actuated_link = np.zeros((step,len(Updated_actuated_link),len(Updated_actuated_link[0])))
    Step_para[0,:] = Initial_para
    Step_actuated_link[0,:] = Updated_actuated_link


    for i in range(step-1):
        phi = phi + 360/step
        Updated_actuated_link = update_actuated_link(Actuate_link, actuator, Actuate_inx, ja, jc, jt, Bi_link_len, Tri_link_len,phi)
        result_temp = root(F_wrapped, Step_para[i], args=(Updated_actuated_link.flatten()), tol = Root_accuracy)
        temp_para, condition = result_temp.x, result_temp.success   

        if condition == True:
            #print('-------------------------')
            Step_para[i+1,:] = temp_para
            Step_actuated_link[i+1,:] = Updated_actuated_link
        else:
            break
            
    "calculate backward from the initial point"
    phi = 0
    #print("step1:",i)
    step2 = step-i
    Step_para_2= np.zeros((step,len(Unknown_paras))) # store every position of point_4
    Step_actuated_link_2 = np.zeros((step,len(Updated_actuated_link),len(Updated_actuated_link[0])))
    
    Step_para_2[0,:] = Initial_para
    Step_actuated_link_2[0,:] = Updated_actuated_link
    
    for i in range(step2-1):
        phi = phi - 360/step
        Updated_actuated_link = update_actuated_link(Actuate_link, actuator, Actuate_inx, ja, jc, jt, Bi_link_len, Tri_link_len,phi)
        result_temp = root(F_wrapped, Step_para_2[i], args=(Updated_actuated_link.flatten()), tol = Root_accuracy)
        temp_para, condition = result_temp.x, result_temp.success   

        if condition == True:
            #print('-------------------------')
            Step_para_2[i+1,:] = temp_para
            Step_actuated_link_2[i+1,:] = Updated_actuated_link
        else:
            break
            
    "flip the second round of calculation and final combine"    
    zero_inx = 0
    for i in range(len(Step_para_2)-1):
        if abs(np.sum(Step_actuated_link_2[i,0] - np.array([0,0,0])))< 1e-5 and  abs(np.sum(Step_actuated_link_2[i+1,0] - np.array([0,0,0])))< 1e-5:
            zero_inx = i-1
            break
        zero_inx = i

    Step_para_2 = Step_para_2[0:zero_inx,:]
    Step_actuated_link_2 = Step_actuated_link_2[0:zero_inx,:]
    Step_para_2_r = np.flip(Step_para_2,0)
    Step_actuated_link_2_r = np.flip(Step_actuated_link_2,0)

    
    zero_inx = 0
    for i in range(len(Step_para)-1):
        if abs(np.sum(Step_actuated_link[i,0] - np.array([0,0,0])))< 1e-5 and  abs(np.sum(Step_actuated_link[i+1,0] - np.array([0,0,0])))< 1e-5:
            zero_inx = i-1
            break
        zero_inx = i

    Step_para = Step_para[0:zero_inx,:]
    Step_actuated_link = Step_actuated_link[0:zero_inx,:]

    Final_step_para = np.concatenate((Step_para_2_r, Step_para), axis=0)
    Final_step_actuated_link = np.concatenate((Step_actuated_link_2_r, Step_actuated_link), axis=0)
    
    return Final_step_para, Final_step_actuated_link


def Coupler_in_Rigidbody (ref_pos_1, ref_pos_2, ref_axis, coupler, frame1):
    
    center_f2 = 0.5*(ref_pos_1-ref_pos_2) + ref_pos_2;
    #x axis
    x_rigid = np.cross((ref_pos_2 - ref_pos_1),ref_axis);
    x_rigid = x_rigid/np.linalg.norm(x_rigid);
    #y axis
    y_rigid = ref_pos_2 - ref_pos_1;
    y_rigid = y_rigid/np.linalg.norm(y_rigid);
    #z axis
    z_rigid = np.cross(x_rigid,y_rigid);
    z_rigid = z_rigid/np.linalg.norm(z_rigid);
    
    frame2 = np.array([x_rigid,y_rigid,z_rigid])
    
    rotation_f1_f2 = find_rotation_matrix(frame1,frame2)
    
    coupler_f2 = np.zeros((1,3), dtype=np.float32)
    coupler_f2 = np.transpose(np.dot(np.linalg.inv(rotation_f1_f2), (coupler-center_f2).T))
    #print("joint 6 in f2:", coupler_f2)

    return coupler_f2

def Coupler_Path_Find(ref_pos_1, ref_pos_2, ref_axis, coupler_relative_pos, frame1):
    
    #calculate frame 2
    center_f2 = 0.5*(ref_pos_1-ref_pos_2) + ref_pos_2;
    #x axis
    x_rigid = np.cross((ref_pos_2 - ref_pos_1),ref_axis);
    x_rigid = x_rigid/np.linalg.norm(x_rigid);
    #y axis
    y_rigid = ref_pos_2 - ref_pos_1;
    y_rigid = y_rigid/np.linalg.norm(y_rigid);
    #z axis
    z_rigid = np.cross(x_rigid,y_rigid);
    z_rigid = z_rigid/np.linalg.norm(z_rigid);
    
    frame2 = np.array([x_rigid,y_rigid,z_rigid])
    
    rotation_f1_f2 = find_rotation_matrix(frame1,frame2)

    coupler_f1 = np.zeros((1,3), dtype=np.float32)
    coupler_f1 = np.transpose(np.dot(rotation_f1_f2, (coupler_relative_pos).T))+center_f2
    
    return coupler_f1

"Find rotation matrix between two frames"
"P2 = R*P1, R = P2*inv(P1), Pi means P in frame i"
def find_rotation_matrix(axis_vectors_F1, axis_vectors_F2):
    axis_vectors_F1 = axis_vectors_F1.T
    axis_vectors_F2 = axis_vectors_F2.T
    R = np.dot(axis_vectors_F2, np.linalg.inv(axis_vectors_F1))
    return R
####################################  Animation  Functions  ##############################################################

# from vpython import *
def animation_setting (jt, jc, ja, bi_links, tri_links):
    joint = []
    for i in range(len(jt)):
        if jt[i] == 0: # R joint
            joint.append(cylinder(pos=vector(jc[i,0],jc[i,1],jc[i,2]), axis=vector(0.5*ja[i,0,0],0.5*ja[i,0,1],0.5*ja[i,0,2]),radius=0.3,color=color.blue))
            joint.append(cylinder(pos=vector(jc[i,0],jc[i,1],jc[i,2]), axis=vector(-0.5*ja[i,0,0],-0.5*ja[i,0,1],-0.5*ja[i,0,2]),radius=0.3,color=color.blue))

        elif jt[i] == 1: # P joint
            joint.append(box(pos=vector(jc[i,0],jc[i,1],jc[i,2]), axis=vector(0.5*ja[i,0,0],0.5*ja[i,0,1],0.5*ja[i,0,2]), length=1, height=0.5, width=0.5, color=color.black))
            joint.append(box(pos=vector(jc[i,0],jc[i,1],jc[i,2]), axis=vector(0.5*ja[i,0,0],0.5*ja[i,0,1],0.5*ja[i,0,2]),length=1, height=0.5, width=0.5, color=color.black))

        elif jt[i] == 3: # C joint
            joint.append(cylinder(pos=vector(jc[i,0],jc[i,1],jc[i,2]), axis=vector(0.3*ja[i,0,0],0.3*ja[i,0,1],0.3*ja[i,0,2]),radius=0.3,color=color.magenta))
            joint.append(cylinder(pos=vector(jc[i,0],jc[i,1],jc[i,2]), axis=vector(-0.3*ja[i,0,0],-0.3*ja[i,0,1],-0.3*ja[i,0,2]),radius=0.3,color=color.magenta))
        
        elif jt[i] == 2: # U joint
            joint.append(cylinder(pos=vector(jc[i,0],jc[i,1],jc[i,2]), axis=vector(0.5*ja[i,0,0],0.5*ja[i,0,1],0.5*ja[i,0,2]),radius=0.3,color=color.red))
            joint.append(cylinder(pos=vector(jc[i,0],jc[i,1],jc[i,2]), axis=vector(-0.5*ja[i,0,0],-0.5*ja[i,0,1],-0.5*ja[i,0,2]),radius=0.3,color=color.red))
            
            joint.append(cylinder(pos=vector(jc[i,0],jc[i,1],jc[i,2]), axis=vector(0.5*ja[i,1,0],0.5*ja[i,1,1],0.5*ja[i,1,2]),radius=0.3,color=color.orange))
            joint.append(cylinder(pos=vector(jc[i,0],jc[i,1],jc[i,2]), axis=vector(-0.5*ja[i,1,0],-0.5*ja[i,1,1],-0.5*ja[i,1,2]),radius=0.3,color=color.orange))

        elif jt[i] == 4: # S joint
            joint.append(sphere(pos=vector(jc[i,0],jc[i,1],jc[i,2]),radius=0.3,color=color.green))
        elif jt[i] == 6: # W joint
            joint.append(sphere(pos=vector(jc[i,0],jc[i,1],jc[i,2]),radius=0.15,color=color.black))
            
    link = []
    for i in range(len(bi_links)):
        inx0, inx1 = bi_links[i,0], bi_links[i,1]
        temp_vector = -vector(jc[inx0,0],jc[inx0,1],jc[inx0,2]) + vector(jc[inx1,0],jc[inx1,1],jc[inx1,2])
        link.append(cylinder(pos=vector(jc[inx0,0],jc[inx0,1],jc[inx0,2]), axis=temp_vector, radius=0.15, color=color.purple))

    for i in range(len(tri_links)):
        inx0, inx1, inx2 = tri_links[i,0], tri_links[i,1], tri_links[i,2]
        temp_vector = -vector(jc[inx0,0],jc[inx0,1],jc[inx0,2]) + vector(jc[inx1,0],jc[inx1,1],jc[inx1,2])
        link.append(cylinder(pos=vector(jc[inx0,0],jc[inx0,1],jc[inx0,2]), axis=temp_vector, radius=0.15, color=color.cyan))
        
        if jt[inx2] == 3: # C joint

            #temp_vector = -vector(jc[inx2,0],jc[inx2,1],jc[inx2,2]) + vector(jc[inx1,0],jc[inx1,1],jc[inx1,2])
            temp_vector = 2*(vector(ja[inx2,0,0], ja[inx2,0,1], ja[inx2,0,2]))
            link.append(cylinder(pos=vector(jc[inx2,0],jc[inx2,1],jc[inx2,2]), axis=temp_vector*3, radius=0.15, color=color.cyan))
            link.append(cylinder(pos=vector(jc[inx2,0],jc[inx2,1],jc[inx2,2]), axis=-temp_vector*3, radius=0.15, color=color.cyan))
        elif jt[inx2] == 1: # P joint

            #temp_vector = -vector(jc[inx2,0],jc[inx2,1],jc[inx2,2]) + vector(jc[inx1,0],jc[inx1,1],jc[inx1,2])
            temp_vector = 2*(vector(ja[inx2,0,0], ja[inx2,0,1], ja[inx2,0,2]))
            link.append(box(pos=vector(jc[inx2,0],jc[inx2,1],jc[inx2,2]), axis=temp_vector, length=10, height=0.25, width=0.25, color=color.cyan))
            link.append(box(pos=vector(jc[inx2,0],jc[inx2,1],jc[inx2,2]), axis=temp_vector, length=10, height=0.25, width=0.25, color=color.cyan))
             
        
    return joint, link

def find_unknowns_inx(ground_link, ground_inx, actuate_link, jt, bi_links):
    joint_inx = [] # [joint inx]
    axis_inx = [] # [joint inx, axis inx]
    for i in range(len(jt)):
        if (i == ground_link[0]) or (i == ground_link[1]): # skip ground link except the axis of u joint
            if jt[i] == 2:
                for j in range(len(bi_links)):
                    if bi_links[j,0] == i or bi_links[j,1] == i: 
                        if ground_inx > j: #u已被用过，地杆占角u_2， u_1未知
                            axis_inx.append(np.array([i,0]))
                        else:
                            axis_inx.append(np.array([i,1]))
                        break

            else:
                continue
                
        elif i == actuate_link[1]: #skip the actuated link except the axis of u joint
            if jt[i] == 2:
                axis_inx.append(np.array([i,1]))
            else:
                continue
                
        elif (len(actuate_link) == 3) and (i == actuate_link[2]): # if the actuated link is a tri_link, then the unkown will be the position of c or p joint
            #while the c/p joint's axis will be known
            "it's risky here because if len(actuate_link)=2, then actuate_link[2] won't exist. But it seems python will calculate"
            "the part before and first..."
            #if i == actuate_link[2]:
            #print("skip:",i)
            joint_inx.append(i)
                           
        elif jt[i] == 2: #u
            joint_inx.append(i)
            axis_inx.append(np.array([i,0]))
            axis_inx.append(np.array([i,1]))

        elif jt[i] == 6 or jt[i] == 4: # w or s
            joint_inx.append(i)
        else:
            joint_inx.append(i)
            axis_inx.append(np.array([i,0]))
            
    joint_inx = np.array(joint_inx)
    axis_inx = np.array(axis_inx)
    return joint_inx, axis_inx

def animation_update_jc_ja (actuate_link, an_update_joint, an_update_axis, jc, ja, step_para, step_actuated_link):
    an_updated_jc = copy.deepcopy(jc)
    an_updated_ja = copy.deepcopy(ja)
    inx = 0
    
    joint_inx = 0
    for i in range(len(jc)): # go through all joints one by one
        if joint_inx < len(an_update_joint):
            if an_update_joint[joint_inx] == i: # if joint i needs to be updated
                an_updated_jc[i] = step_para[inx:inx+3]
                joint_inx += 1
                inx += 3
        for j in range(len(an_update_axis)):
            if an_update_axis[j,0] == i:
                if an_update_axis[j,1] == 0:
                    an_updated_ja[i,0] = step_para[inx:inx+3]
                    inx += 3
                if an_update_axis[j,1] == 1:
                    an_updated_ja[i,1] = step_para[inx:inx+3]
                    inx += 3
                    
    an_updated_jc[actuate_link[1]] = step_actuated_link[0]
    an_updated_ja[actuate_link[1],0] = step_actuated_link[1]
    if len(actuate_link)==3:
        an_updated_ja[actuate_link[2],0] = step_actuated_link[2]

    return an_updated_jc, an_updated_ja

def animation_update (an_joint, an_link, jt, up_jc, up_ja, bi_links, tri_links):
    joint_inx = 0
    for i in range(len(jt)):
        if jt[i] == 0 or jt[i] == 1: # R, P joint
            an_joint[joint_inx].pos = vector(up_jc[i,0],up_jc[i,1],up_jc[i,2])
            an_joint[joint_inx].axis = vector(0.5*up_ja[i,0,0],0.5*up_ja[i,0,1],0.5*up_ja[i,0,2])
            joint_inx += 1
            an_joint[joint_inx].pos = vector(up_jc[i,0],up_jc[i,1],up_jc[i,2])
            an_joint[joint_inx].axis = vector(-0.5*up_ja[i,0,0],-0.5*up_ja[i,0,1],-0.5*up_ja[i,0,2])
            joint_inx += 1
#             print("r or p", joint_inx)
            
        elif jt[i] == 1: # P joint
            an_joint[joint_inx].pos = vector(up_jc[i,0],up_jc[i,1],up_jc[i,2])
            an_joint[joint_inx].axis = vector(0.5*up_ja[i,0,0],0.5*up_ja[i,0,1],0.5*up_ja[i,0,2])
            an_joint[joint_inx].length = 1
            joint_inx += 1
            an_joint[joint_inx].pos = vector(up_jc[i,0],up_jc[i,1],up_jc[i,2])
            an_joint[joint_inx].axis = vector(-0.5*up_ja[i,0,0],-0.5*up_ja[i,0,1],-0.5*up_ja[i,0,2])
            an_joint[joint_inx].length = 1
            joint_inx += 1
        
        elif jt[i] == 2: # U joint
            an_joint[joint_inx].pos = vector(up_jc[i,0],up_jc[i,1],up_jc[i,2])
            an_joint[joint_inx].axis = vector(0.5*up_ja[i,0,0],0.5*up_ja[i,0,1],0.5*up_ja[i,0,2])            
            joint_inx += 1
            an_joint[joint_inx].pos = vector(up_jc[i,0],up_jc[i,1],up_jc[i,2])
            an_joint[joint_inx].axis = vector(-0.5*up_ja[i,0,0],-0.5*up_ja[i,0,1],-0.5*up_ja[i,0,2])            
            joint_inx += 1
#             print("u1", joint_inx)
            
            an_joint[joint_inx].pos = vector(up_jc[i,0],up_jc[i,1],up_jc[i,2])
            an_joint[joint_inx].axis = vector(0.5*up_ja[i,1,0],0.5*up_ja[i,1,1],0.5*up_ja[i,1,2])            
            joint_inx += 1       
            an_joint[joint_inx].pos = vector(up_jc[i,0],up_jc[i,1],up_jc[i,2])
            an_joint[joint_inx].axis = vector(-0.5*up_ja[i,1,0],-0.5*up_ja[i,1,1],-0.5*up_ja[i,1,2])            
            joint_inx += 1 
#             print("u2", joint_inx)
            
        elif jt[i] == 3: # C joint
            an_joint[joint_inx].pos = vector(up_jc[i,0],up_jc[i,1],up_jc[i,2])
            an_joint[joint_inx].axis = vector(0.3*up_ja[i,0,0],0.3*up_ja[i,0,1],0.3*up_ja[i,0,2])
            joint_inx += 1
            an_joint[joint_inx].pos =vector(up_jc[i,0],up_jc[i,1],up_jc[i,2])
            an_joint[joint_inx].axis = vector(-0.3*up_ja[i,0,0],-0.3*up_ja[i,0,1],-0.3*up_ja[i,0,2])
            joint_inx += 1
#             print("c", joint_inx)

        elif jt[i] == 4: # S joint
            an_joint[joint_inx].pos = vector(up_jc[i,0],up_jc[i,1],up_jc[i,2])
            joint_inx += 1 
#             print("s", joint_inx)
        elif jt[i] == 6: # W joint
            an_joint[joint_inx].pos = vector(up_jc[i,0],up_jc[i,1],up_jc[i,2])
            joint_inx += 1 
#             print("w", joint_inx)
            
    link_inx = 0
    for i in range(len(bi_links)):
        inx0, inx1 = bi_links[i,0], bi_links[i,1]
        temp_vector = -vector(up_jc[inx0,0],up_jc[inx0,1],up_jc[inx0,2]) + vector(up_jc[inx1,0],up_jc[inx1,1],up_jc[inx1,2])
        an_link[link_inx].pos = vector(up_jc[inx0,0],up_jc[inx0,1],up_jc[inx0,2])
        an_link[link_inx].axis = temp_vector
        link_inx += 1
        
    for i in range(len(tri_links)):
        inx0, inx1, inx2 = tri_links[i,0], tri_links[i,1], tri_links[i,2]
        temp_vector = -vector(up_jc[inx0,0],up_jc[inx0,1],up_jc[inx0,2]) + vector(up_jc[inx1,0],up_jc[inx1,1],up_jc[inx1,2])
        an_link[link_inx].pos = vector(up_jc[inx0,0],up_jc[inx0,1],up_jc[inx0,2])
        an_link[link_inx].axis = temp_vector
        link_inx += 1
        
        if jt[inx2] == 3: # C joint
        #temp_vector = -vector(up_jc[inx2,0],up_jc[inx2,1],up_jc[inx2,2]) + vector(up_jc[inx1,0],up_jc[inx1,1],up_jc[inx1,2])
            temp_vector = 2*(vector(up_ja[inx2,0,0], up_ja[inx2,0,1], up_ja[inx2,0,2]))
            an_link[link_inx].pos = vector(up_jc[inx1,0],up_jc[inx1,1],up_jc[inx1,2])
            an_link[link_inx].axis = 3*temp_vector
            link_inx += 1

            an_link[link_inx].pos = vector(up_jc[inx1,0],up_jc[inx1,1],up_jc[inx1,2])
            an_link[link_inx].axis = -3*temp_vector
            link_inx += 1
        
        elif jt[inx2] == 1: # P joint
        #temp_vector = -vector(up_jc[inx2,0],up_jc[inx2,1],up_jc[inx2,2]) + vector(up_jc[inx1,0],up_jc[inx1,1],up_jc[inx1,2])
            temp_vector = 2*(vector(up_ja[inx2,0,0], up_ja[inx2,0,1], up_ja[inx2,0,2]))
            an_link[link_inx].pos = vector(up_jc[inx1,0],up_jc[inx1,1],up_jc[inx1,2])
            an_link[link_inx].axis = temp_vector
            an_link[link_inx].length = 20
            link_inx += 1

            an_link[link_inx].pos = vector(up_jc[inx1,0],up_jc[inx1,1],up_jc[inx1,2])
            an_link[link_inx].axis = -temp_vector
            an_link[link_inx].length = 20
            link_inx += 1
    return an_joint, an_link
        

