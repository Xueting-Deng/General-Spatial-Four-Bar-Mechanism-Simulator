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

def RigidbodyMesh (ref_pos_1, ref_pos_2, ref_axis):
    
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

# R-0, P-1, U-2, C-3, S-4, E (plane)-5, W (welded) -6

# Find both binary and tribary links
def find_link(jj, jt, actuator):
    bi_links = []
    tri_links = []
    ground_link = []
    
    ground_inx = []
    actuated_inx = []
    temp_tri = [0,0,0]
    
    for i in range(len(jj[0])):
        if jt[i] == 6:
            temp_tri[1] = i
            temp_inx = 0
            for j in range(len(jj[0])):
                if jj[i,j] == 1 or jj[i,j]==2:
                    temp_tri[temp_inx] = j
                    temp_inx = 2
            tri_links.append(temp_tri)
        
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
        angle_1 = angle_1 if abs(angle_1-1) > 1e-6 else 1
        angle_1 = (angle_1 if angle_1 > -np.pi else angle_1+np.pi) if angle_1 < np.pi else angle_1-np.pi;#make sure the ang is in (0,pi)
        angle_1 = np.arccos(angle_1)
        
    # the angle of the 2nd joint and the link
    if (all(ja_2 == no_angle)):
        angle_2 = 100   #100 here means the angle doesn't exist!!
        angle_1_2 = 100 # if angle_2 doesn't exist, the angle between two joints doesn't exist
    else:
        angle_2 = np.dot(ja_2, jc_1 - jc_2)/np.linalg.norm(jc_2 - jc_1)
        angle_2 = angle_2 if abs(angle_2-1) > 1e-6 else 1
        angle_2 = (angle_2 if angle_2 > -np.pi else angle_2+np.pi) if angle_2 < np.pi else angle_2-np.pi;#make sure the ang is in (0,pi)
        angle_2 = np.arccos(angle_2)
        
    # the angle between two joint
    if angle_1_2 != 100:
        angle_1_2 = np.dot(ja_1, ja_2)
        angle_1_2 = angle_1_2 if abs(angle_1_2-1) > 1e-6 else 1
        angle_1_2 = (angle_1_2 if angle_1_2 > -np.pi else angle_1_2+np.pi) if angle_1_2 < np.pi else angle_1_2-np.pi;#make sure the ang is in (0,pi)
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
        inx0, inx1 = bi_links[i,0], bi_links[i,1]
        bi_link_len[i] = np.linalg.norm(jc[inx0]-jc[inx1])
        
        if jt[inx0] != 2 and jt[inx1] !=2: # both joints are not u joint
            bi_link_ang[i] = compute_angle(jc[inx0], jc[inx1], ja[inx0,0], ja[inx1,0])
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
        inx0, inx1 = bi_links[i,0], bi_links[i,1]
        #f = (un[inx0,0:3]-un[inx1,0:3]).dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]**2 # length constraint
        f = sp.sqrt((un[inx0,0:3]-un[inx1,0:3]).dot(un[inx0,0:3]-un[inx1,0:3])) - bi_link_len[i] # length constraint
        F.append(f)
        
        if jt[inx0] != 2 and jt[inx1] !=2: # both joints are not u joint
            f = un[inx0,3:6].dot(un[inx1,0:3]-un[inx0,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,0])
            F.append(f)
            f = un[inx1,3:6].dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,1])
            F.append(f)
            f = un[inx0,3:6].dot(un[inx1,3:6]) - np.cos(bi_link_ang[i,2])
            F.append(f)
            
            
        elif jt[inx0] == 2 and jt[inx1] !=2:# j1 is u joint and j2 is not
            if u_trigger[inx0] == 0: # u joint's first angle always belongs to the first link when calculate in order
                f = un[inx0,3:6].dot(un[inx1,0:3]-un[inx0,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,0])
                F.append(f)
                f = un[inx1,3:6].dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,1])
                F.append(f)
                f = un[inx0,3:6].dot(un[inx1,3:6]) - np.cos(bi_link_ang[i,2])
                F.append(f)
                u_trigger[inx0] = 1 # mark this u joint's first angle is calculted          
            else:
                f = un[inx0,6:9].dot(un[inx1,0:3]-un[inx0,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,0])
                F.append(f)
                f = un[inx1,3:6].dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,1])
                F.append(f)
                f = un[inx0,6:9].dot(un[inx1,3:6]) - np.cos(bi_link_ang[i,2])
                F.append(f)

            
        elif jt[inx0] != 2 and jt[inx1] ==2:# j1 is not u joint and j2 is u joint
            if u_trigger[inx1] == 0: # u joint's first angle always belongs to the first link when calculate in order
                f = un[inx0,3:6].dot(un[inx1,0:3]-un[inx0,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,0])
                F.append(f)
                f = un[inx1,3:6].dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,1])
                F.append(f)
                f = un[inx0,3:6].dot(un[inx1,3:6]) - np.cos(bi_link_ang[i,2])
                F.append(f)
                u_trigger[inx1] = 1 # mark this u joint's first angle is calculted
            else:
                f = un[inx0,3:6].dot(un[inx1,0:3]-un[inx0,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,0])
                F.append(f)
                f = un[inx1,6:9].dot(un[inx0,0:3]-un[inx1,0:3]) - bi_link_len[i]*np.cos(bi_link_ang[i,1])
                F.append(f)
                f = un[inx0,3:6].dot(un[inx1,6:9]) - np.cos(bi_link_ang[i,2])
                F.append(f)

            
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
                f = un[inx0,6:0].dot(un[inx1,6:9]) - np.cos(bi_link_ang[i,2])
                F.append(f)
    return F

def compute_tri_link_length_angle (lj, jt, tri_links, jc, ja):
    tri_link_len = np.zeros(len(tri_links))
    tri_link_ang = np.zeros((len(tri_links),3))
    for i in range(len(tri_links)):
        inx0, inx1, inx2 = tri_links[i,0],tri_links[i,1],tri_links[i,2] #(l_012, 1 is always w joint)
        tri_link_len[i] = np.linalg.norm(jc[inx0]-jc[inx1]) # link length of l01
        
        if jt[inx0] == 2: # if j0 is u joint
            angle_1 = np.dot(ja[inx0,1], jc[inx1] - jc[inx0])/np.linalg.norm(jc[inx1]-jc[inx0])            
            angle_2 = np.dot(ja[inx0,1], ja[inx2,0])
        else:
            angle_1 = np.dot(ja[inx0,0], jc[inx1] - jc[inx0])/np.linalg.norm(jc[inx1]-jc[inx0])
            angle_2 = np.dot(ja[inx0,0], ja[inx2,0])
        
        angle_3 = np.dot(jc[inx0] - jc[inx1], ja[inx2,0])/np.linalg.norm(jc[inx0] - jc[inx1])
        angle_3 = (angle_3 if angle_3 > -np.pi else angle_3+np.pi) if angle_3 < np.pi else angle_3-np.pi;#make sure the ang is in (0,pi)
        angle_3 = np.arccos(angle_3)         
        
        angle_1 = (angle_1 if angle_1 > -np.pi else angle_1+np.pi) if angle_1 < np.pi else angle_1-np.pi;#make sure the ang is in (0,pi)
        angle_1 = np.arccos(angle_1) 
        angle_2 = (angle_2 if angle_2 > -np.pi else angle_2+np.pi) if angle_2 < np.pi else angle_2-np.pi;#make sure the ang is in (0,pi)
        angle_2 = np.arccos(angle_2) 
        
        tri_link_ang[i,0] = angle_1
        tri_link_ang[i,1] = angle_2
        tri_link_ang[i,2] = angle_3
        
    return tri_link_len, tri_link_ang

def build_tri_link_constraints (un, jt, tri_links, tri_link_len, tri_link_ang):
    F=[]
    for i in range(len(tri_links)):
        inx0, inx1, inx2 = tri_links[i,0],tri_links[i,1],tri_links[i,2] #(l_012, 1 is always w joint)
        
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
        u12 = un[inx2,0:3] - un[inx1,0:3]
        u2 = un[inx2,3:6]
        cross_1 = u12[1]*u2[2] - u12[2]*u2[1]
        cross_2 = u12[0]*u2[2] - u12[2]*u2[0]
        cross_3 = u12[0]*u2[1] - u12[1]*u2[0]
        
        #if abs(u2[2]) <= 10e-8:
#         F.append(u12[2])
#         F.append(cross_3)
        #else:
        F.append(cross_1)
        F.append(cross_2)
        
    return F

def compute_type_angle(jt, ja):
    type_angle = np.zeros(len(jt))
    for i in range(len(jt)):
        if jt[i] == 2: # if the joint type is U, calculate the joint angle between two axis
            angle = np.dot(ja[i,0], ja[i,1])
            angle = angle if abs(abs(angle)-1) > 1e-5 else 1
            angle = (angle if angle > -np.pi else angle+np.pi) if angle < np.pi else angle-np.pi;#make sure the ang is in (0,pi)
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
            f = sp.sqrt(un[i,3:6].dot(un[i,3:6])) - 1
            F.append(f)
            f = sp.sqrt(un[i,6:9].dot(un[i,6:9])) - 1
            F.append(f)  
        elif jt[i] == 6: # W joint
            continue
        else:
            f = sp.sqrt(un[i,3:6].dot(un[i,3:6])) - 1
            F.append(f)     
    return F

def compute_p_angle(jt, bi_links, tri_links, jc, ja):
    
    Angle_p = []
    for i in range(len(tri_links)):
        inx0, inx1, inx2 = tri_links[i,0],tri_links[i,1],tri_links[i,2] #(l_012, 1 is always w joint)
#         if jt[inx0] == 1: # if j0 is p joint
#             print("inx0 is p")
#             for j in range(len(bi_links)):
#                 if inx0 == bi_links[j,0] or inx0 == bi_links[j,1]:
#                     print("bilink,",bi_links[j])
#                     inx_b1, inx_b2 = bi_links[j,0], bi_links[j,1]
#                     angle_p = np.dot(jc[inx_b2] - jc[inx_b1], jc[inx2] - jc[inx1])/(np.linalg.norm(jc[inx2]-jc[inx1])*np.linalg.norm(jc[inx_b2] - jc[inx_b1]))
                   
#                     angle_p = (angle_p if angle_p > -np.pi else angle_p+np.pi) if angle_p < np.pi else angle_p-np.pi;#make sure the ang is in (0,pi)
#                     angle_p = np.arccos(angle_p)
#                     Angle_p.append(angle_p) 
        if jt[inx2] == 1: # if j2 is p joint
            #print("inx2 is p")
            for j in range(len(bi_links)):
                if inx2 == bi_links[j,0] or inx2 == bi_links[j,1]:
                    #print("bilink,",bi_links[j])
                    inx_b1, inx_b2 = bi_links[j,0], bi_links[j,1]
                    angle_p = np.dot(jc[inx_b2] - jc[inx_b1], jc[inx0] - jc[inx1])/(np.linalg.norm(jc[inx0]-jc[inx1])*np.linalg.norm(jc[inx_b2] - jc[inx_b1]))
            
                    angle_p = (angle_p if angle_p > -np.pi else angle_p+np.pi) if angle_p < np.pi else angle_p-np.pi;#make sure the ang is in (0,pi)
                    angle_p = np.arccos(angle_p)
                    Angle_p.append(angle_p)
    Angle_p = np.array(Angle_p)
        
        
    return Angle_p

def build_p_constraints(jt, bi_links, bi_link_len, tri_links, tri_link_len, jc, ja, angle_p, un):
    F = []
    inx_p = 0
    
    for i in range(len(tri_links)):
        inx0, inx1, inx2 = tri_links[i,0],tri_links[i,1],tri_links[i,2] #(l_012, 1 is always w joint, 2 is p/c joint)

        if jt[inx2] == 1: # if j2 is p joint
            for j in range(len(bi_links)):
                if inx2 == bi_links[j,0] or inx2 == bi_links[j,1]:
                    inx_b1, inx_b2 = bi_links[j,0], bi_links[j,1]
                    f = (un[inx_b2,0:3]-un[inx_b1,0:3]).dot(un[inx0,0:3]-un[inx1,0:3]) - tri_link_len[i]*bi_link_len[j]*np.cos(angle_p[inx_p])
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

def update_actuated_link (actuate_link, actuator, actuate_inx, ja, jc, jt, bi_link_len, phi):
    
    # using quaternion to calculate new position
    #if actuate_link[0] == actuator: 
    actuate_vector = (jc[actuate_link[1],:] - jc[actuate_link[0],:])/np.linalg.norm(jc[actuate_link[1],:] - jc[actuate_link[0],:])
    actuate_vector =  np.append(0, actuate_vector)
    #else:
    #    print('Actuator format wrong. Please let the first element of the actuate_link be where the actuator located on')

    theta = np.pi*phi/360; #thera is half of phi  
    actuator_axis = ja[actuator,0]
    
    rotation_quater = np.array([np.cos(theta), np.sin(theta)*actuator_axis[0], np.sin(theta)*actuator_axis[1], np.sin(theta)*actuator_axis[2]]) #the quaternion of rotation
    
    new_actuated_vec = QuaterTimes(QuaterTimes(rotation_quater,actuate_vector),QuaterConj(rotation_quater))
    new_pos = new_actuated_vec[1:4]*bi_link_len[actuate_inx[0]] + jc[actuator] # new position of the joint connected to the actuated joint
    new_axis = np.array([0,0,0,0], dtype=np.float64)
    
    #if jt[actuate_link[1]] ==  2: # if it's u joint
        #new_axis = rotate_vector(ja[actuate_link[1],0], actuator_axis, theta*2) # new vector of the joint connected to the actuated joint if it has
    new_axis = QuaterTimes(QuaterTimes(rotation_quater,np.append(0,ja[actuate_link[1],0])),QuaterConj(rotation_quater))
        #????????updated_actuated_link = np.vstack((new_pos,new_axis[1:]))
        #除了S joint以外，R, P, U, C的new_axis都应该被更新
    
    
    updated_actuated_link = np.vstack((new_pos,new_axis[1:]))
    
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
                    if ground_inx < j: # U还没有被用过，ground用u1，u1是known 
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
        if isinstance(initial_constraint[i],np.float64):
            if abs(initial_constraint[i]) < 1e-5:
                removal.append(i)
            else:
                print('The', i, "th constraint is calculatable but not zero, the constraint equations have error")
        elif isinstance(initial_constraint[i], sp.core.add.Add) == False:
            removal.append(i)   
    final_constraint = np.delete(initial_constraint,removal)
    return final_constraint

def find_real_unknowns_and_initial_guess(ground_link, ground_inx, actuate_link, jt, un, ja, jc, bi_links, updated_actuated_link):
    real_un = []
    initial_guess = []
    for i in range(len(jt)):
        if (i == ground_link[0]) or (i == ground_link[1]): # skip ground link except the axis of u joint
            if jt[i] == 2: 
                #real_un.append(un[i,6:9])
                #initial_guess.append(ja[i,1])
                for j in range(len(bi_links)):
                    if bi_links[j,0] == i or bi_links[j,1] == i: 
                        if ground_inx < j: #u没被用过
                            real_un.append(un[i,6:9]) #u2未知
                            initial_guess.append(ja[i,1])
                        else:
                            real_un.append(un[i,3:6])
                            initial_guess.append(ja[i,0])
                        break

            else:
                continue
            #print(real_un)
                
        elif i == actuate_link[1]: #skip the actuated link except the axis of u joint
            if jt[i] == 2:
                real_un.append(un[i,6:9])
                initial_guess.append(ja[i,1]) # update intial guess with the updated actuated joint it's a
            else:
                continue
        
        elif jt[i] == 2: #u
            real_un.append(un[i])
            initial_guess.append(jc[i])
            initial_guess.append(ja[i,0])
            initial_guess.append(ja[i,1])
        elif jt[i] == 6 or jt[i] == 4: # w and s
            real_un.append(un[i,0:3])
            initial_guess.append(jc[i])
        else:
            real_un.append(un[i,0:6])
            initial_guess.append(jc[i])
            initial_guess.append(ja[i,0])
    real_un = np.concatenate(real_un)
    initial_guess = np.concatenate(initial_guess)
    return real_un, initial_guess

def solve_constraint_equations(constraint_equations_sym, unknown_params, guess):
    
    def constraint_eqs(upara):
        #print(np.array([eq(*upara) for eq in constraint_equations]))
        return np.array([eq(*upara) for eq in constraint_equations])
    
    constraint_equations = [sp.lambdify(unknown_params, eq, 'numpy') for eq in constraint_equations_sym]

    result = root(constraint_eqs, guess)#,method='krylov')
    #print(result)

    
    if result.success:
        solution = result.x
        #print(result.fun)
        return solution, True
    else:
        mes = result.message
        #print(result.fun)
        #print(mes)
        
        return None, False    

def find_actuator_inx (un, actuate_link):
    un_actuator_inx = np.hstack((un[actuate_link[1], 0:3], un[actuate_link[1], 3:6]))
    return un_actuator_inx

def assign_value_for_actuator(final_constraint, un_actuator_inx, updated_actuated_link):
    updated_constraint = []
    values = assign_values(un_actuator_inx,updated_actuated_link.flatten())
    for i in range(len(final_constraint)):
        updated_constraint.append(final_constraint[i].xreplace(values))

    updated_constraint = np.array(updated_constraint)
    updated_constraint = remove_unvalid_constraint(updated_constraint)
    return updated_constraint

# Animation
def animation_setting (jt, jc, ja, bi_links, tri_links):
    joint = []
    for i in range(len(jt)):
        if jt[i] == 0: # R joint
            joint.append(cylinder(pos=vector(jc[i,0],jc[i,1],jc[i,2]), axis=vector(0.5*ja[i,0,0],0.5*ja[i,0,1],0.5*ja[i,0,2]),radius=0.3,color=color.blue))
            joint.append(cylinder(pos=vector(jc[i,0],jc[i,1],jc[i,2]), axis=vector(-0.5*ja[i,0,0],-0.5*ja[i,0,1],-0.5*ja[i,0,2]),radius=0.3,color=color.blue))

        elif jt[i] == 1: # P joint
            joint.append(box(pos=vector(jc[i,0],jc[i,1],jc[i,2]), axis=vector(0.5*ja[i,0,0],0.5*ja[i,0,1],0.5*ja[i,0,2]), length=1, height=0.5, width=0.5, color=color.black))
            joint.append(box(pos=vector(jc[i,0],jc[i,1],jc[i,2]), axis=vector(-0.5*ja[i,0,0],-0.5*ja[i,0,1],-0.5*ja[i,0,2]),length=1, height=0.5, width=0.5, color=color.black))

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
            link.append(cylinder(pos=vector(jc[inx2,0],jc[inx2,1],jc[inx2,2]), axis=temp_vector, radius=0.15, color=color.cyan))
            link.append(cylinder(pos=vector(jc[inx2,0],jc[inx2,1],jc[inx2,2]), axis=-temp_vector, radius=0.15, color=color.cyan))
        elif jt[inx2] == 1: # P joint

            #temp_vector = -vector(jc[inx2,0],jc[inx2,1],jc[inx2,2]) + vector(jc[inx1,0],jc[inx1,1],jc[inx1,2])
            temp_vector = 2*(vector(ja[inx2,0,0], ja[inx2,0,1], ja[inx2,0,2]))
            link.append(box(pos=vector(jc[inx2,0],jc[inx2,1],jc[inx2,2]), axis=temp_vector, length=10, height=0.25, width=0.25, color=color.cyan))
            link.append(box(pos=vector(jc[inx2,0],jc[inx2,1],jc[inx2,2]), axis=-temp_vector, length=10, height=0.25, width=0.25, color=color.cyan))
             
        
    return joint, link

def find_unknowns_inx(ground_link, ground_inx, actuate_link, jt, bi_links):
    joint_inx = [] # [joint inx]
    axis_inx = [] # [joint inx, axis inx]
    for i in range(len(jt)):
        if (i == ground_link[0]) or (i == ground_link[1]): # skip ground link except the axis of u joint
            if jt[i] == 2:
                for j in range(len(bi_links)):
                    if bi_links[j,0] == i or bi_links[j,1] == i: 
                        if ground_inx < j: #u没被用过
                            axis_inx.append(np.array([i,1])) #u1已知， u2未知
                        else:
                            axis_inx.append(np.array([i,0]))
                        break

            else:
                continue
                
        elif i == actuate_link[1]: #skip the actuated link except the axis of u joint
            if jt[i] == 2:
                axis_inx.append(np.array([i,1]))
            else:
                continue
        
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
                    
    an_updated_jc[Actuate_link[1]] = step_actuated_link[0]
    an_updated_ja[Actuate_link[1],0] = step_actuated_link[1]

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
            an_joint[joint_inx].pos =vector(up_jc[i,0],up_jc[i,1],up_jc[i,2])
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
            an_link[link_inx].axis = temp_vector
            link_inx += 1

            an_link[link_inx].pos = vector(up_jc[inx1,0],up_jc[inx1,1],up_jc[inx1,2])
            an_link[link_inx].axis = -temp_vector
            link_inx += 1
        
        elif jt[inx2] == 1: # C joint
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
        
 
import numpy as np
import copy
from scipy.optimize import root
import math
import time
import sys

import sympy as sp
#from tqdm.notebook import tqdm

# R-0, P-1, U-2, C-3, S-4, E (plane)-5, W (welded) -6

# JJ (Joint-Joint) matrix represents the connection between joints, 2 means ground, 1 means connceted. 
# It follows the natural order J0, J1, J2, ...


JJ = np.array([
    [2, 1, 0, 0, 1],
    [1, 0, 1, 0, 0],
    [0, 1, 0, 1, 0],
    [0, 0, 1, 0, 1],
    [1, 0, 0, 1, 2],
], dtype=int)

# actuators? - temporarily actuate J0 - temporarily actuate R joint

if (JJ != (np.transpose(JJ)) ).all():
    print ("Joint-to-Joint matrix is wrong.")
    
Actuator = 0 # Put the actuator at J0
JT = np.array([0,2,6,1,4], dtype=int)

JC = np.array([
    [0, 0, 0],#j0
    [0, 2, 0], #j1
    [3, 3, 0], #j2
    [3, 3, 3], #j3
    [5, 1, 6], #j4
], dtype=np.float64)

JA = np.array([
    [[0, 0, 1], [0, 0, 0]], # R
    [[np.sqrt(2)/2, 0, np.sqrt(2)/2], [0, 1, 0]], # u
    [[0, 0, 0], [0, 0, 0]], # w
    [[0,0,1], [0, 0, 0]], # p
    [[0, 0, 0], [0, 0, 0]], # s
], dtype=np.float64)

# Calculate the configuration
Bi_links, Tri_links, Ground_link, Ground_inx, Actuate_link, Actuate_inx = find_link(JJ, JT, Actuator)
LJ = find_link_joint_table(JJ, Bi_links, Tri_links)
UN = generate_unknown_table(JJ)
Bi_link_len, Bi_link_ang= compute_bi_link_length_angle (LJ, JT, Bi_links, JC, JA)
Tri_link_len, Tri_link_ang = compute_tri_link_length_angle (LJ, JT, Tri_links, JC, JA)
Spacial_p_ang = compute_p_angle(JT, Bi_links, Tri_links, JC, JA)
Type_angle = compute_type_angle(JT, JA)
print('Bi_links:\n', Bi_links)
print('Tri_links:\n', Tri_links)
print('Ground_link:', Ground_link)
print('Ground_inx:', Ground_inx)
print('Actuate_link:', Actuate_link)
print('Actuate_inx:', Actuate_inx)
print('Bi_link_len:', Bi_link_len)
print('Bi_link_ang:', Bi_link_ang)
print('Tri_link_len:', Tri_link_len)
print('Tri_link_ang:', Tri_link_ang)
print('Special_P_ang:', Spacial_p_ang)
print('Type_ang:', Type_angle)

#Build and solve constraint equations
Updated_UN = assign_value_for_knowns(Ground_link, Ground_inx, Actuate_link, JT, UN, JA, JC, Bi_links) #update UN
UN_actuator_inx = find_actuator_inx (UN, Actuate_link)

Bi_constraint = build_bi_link_constraints (Updated_UN, JT, Bi_links, Bi_link_len, Bi_link_ang, Ground_link)
Tri_constraint = build_tri_link_constraints (Updated_UN, JT, Tri_links, Tri_link_len, Tri_link_ang)
Type_constraint = build_type_angle_constraints (Updated_UN, Type_angle)
Unit_constraint = build_unit_angle_constraints (Updated_UN, JT)
Special_p_constraint = build_p_constraints(JT, Bi_links, Bi_link_len, Tri_links, Tri_link_len, JC, JA, Spacial_p_ang, Updated_UN)
Final_constraint = remove_unvalid_constraint(np.array(Bi_constraint + Tri_constraint + Type_constraint + Unit_constraint + Special_p_constraint))

# Calculate the initial state
phi = 0
Updated_actuated_link = update_actuated_link(Actuate_link, Actuator, Actuate_inx, JA, JC, JT, Bi_link_len, phi)
Updated_constraint = assign_value_for_actuator(Final_constraint, UN_actuator_inx, Updated_actuated_link)

#Find the initial unknowns and initial guess
Unknown_paras, Initial_guess =find_real_unknowns_and_initial_guess(Ground_link, Ground_inx, Actuate_link, JT, UN, JA, JC, Bi_links, Updated_actuated_link)

Initial_para, condition = solve_constraint_equations(Updated_constraint, Unknown_paras, Initial_guess)
print(condition)

mec_num = 1
step = 360
#mec_data = np.zeros((mec_num*1,len(Unknown_paras),3), dtype=np.float64)
#path_data = np.zeros((mec_num*1,step,3), dtype=np.float64)

Step_para = np.zeros((step,len(Unknown_paras))) # store every position of point_4
Step_actuated_link = np.zeros((step,len(Updated_actuated_link),len(Updated_actuated_link[0])))
Rigid_mesh = np.zeros((27, step, 3))

Step_para[0,:] = Initial_para
Step_actuated_link[0,:] = Updated_actuated_link
Rigid_mesh[:,0,:] = RigidbodyMesh(Step_para[0,0:3],Step_actuated_link[0,0,:],Step_para[0,3:6])

start_time = time.time()
phi = 0

for i in range(step-1):
    phi = phi + 360/step
    Updated_actuated_link = update_actuated_link(Actuate_link, Actuator, Actuate_inx, JA, JC, JT, Bi_link_len, phi)
    Updated_constraint = assign_value_for_actuator(Final_constraint, UN_actuator_inx, Updated_actuated_link)
    Updated_constraint = Updated_constraint

    temp_para, condition = solve_constraint_equations(Updated_constraint, Unknown_paras, Step_para[i])
    if condition == True:
        #print('-------------------------')
        Step_para[i+1,:] = temp_para
        Step_actuated_link[i+1,:] = Updated_actuated_link
        Rigid_mesh[:,i+1,:]  = RigidbodyMesh(Step_para[i+1,0:3],Step_actuated_link[i+1,0,:],Step_para[i+1,3:6])
    else:
        print(i)
        print("no full rotation");
        break
print(i)
print(time.time()-start_time)

# animation
from vpython import *
scene = canvas(width=800,height=500,center=vector(0,0,0),background=color.white);

an_Joint, an_Link = animation_setting (JT, JC, JA, Bi_links, Tri_links)
an_Update_joint, an_Update_axis = find_unknowns_inx(Ground_link, Ground_inx, Actuate_link, JT, Bi_links)

while True:
    for i in range(1,step):
        rate(40)
        an_Update_jc, an_Update_ja = animation_update_jc_ja (Actuate_link, an_Update_joint, an_Update_axis, JC, JA, Step_para[i],Step_actuated_link[i])
        #print(an_Update_jc)
        #print(an_Update_ja)
        animation_update(an_Joint, an_Link, JT, an_Update_jc, an_Update_ja, Bi_links, Tri_links)