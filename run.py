import numpy as np
import math
import os
from fractions import Fraction
from chembondpy import bonder
def operation(Matrix,op_all):
    MT=[]
    for op in op_all:
        ops=op.replace(' ','').split(',')
        for no in ['x','y','z']:
            if no in ops[0]:
                X=ops[0].lstrip().split(no)
        for no in ['x','y','z']:
            if no in ops[1]:
                Y=ops[1].lstrip().split(no)
        for no in ['x','y','z']:
            if no in ops[2]:
                Z=ops[2].lstrip().split(no)
        MM=Matrix.copy()
        matrix1=Matrix.copy()
        if ops[0].replace('-','')[0]=='x':
            if X[0]=='-':
                if X[1]!='':
                    matrix1[:,0]=-matrix1[:,0]+Fraction(X[1])
                else:
                    matrix1[:,0]=-matrix1[:,0]
            else:
                if X[1]!='':
                    matrix1[:,0]=matrix1[:,0]+Fraction(X[1])
                else:
                    matrix1[:,0]=matrix1[:,0]   

        if ops[0].replace('-','')[0]=='y':
            if X[0]=='-':
                if X[1]!='':
                    matrix1[:,0]=-matrix1[:,1]+Fraction(X[1])
                else:
                    matrix1[:,0]=-matrix1[:,1]
             
            else:
                if X[1]!=''  :            
                    matrix1[:,0]=matrix1[:,1]+Fraction(X[1])
                else:
                    matrix1[:,0]=matrix1[:,1]

        if ops[0].replace('-','')[0]=='z':
            if X[0]=='-':
                if X[1]!='':
                    matrix1[:,0]=-matrix1[:,2]+Fraction(X[1])
                else:
                    matrix1[:,0]=-matrix1[:,2]
            else:
                if X[1]!='': 
                    matrix1[:,0]=matrix1[:,2]+Fraction(X[1])
                else:
                    matrix1[:,0]=matrix1[:,2]
        MM[:,0]=matrix1[:,0]
        matrix1=Matrix.copy()
        if ops[1].replace('-','')[0]=='y':

            if Y[0]=='-':
                if Y[1]!='':
                    matrix1[:,1]=-matrix1[:,1]+Fraction(Y[1])
                else:
                    matrix1[:,1]=-matrix1[:,1]
            else:
                if Y[1]!='':    
                    matrix1[:,1]=matrix1[:,1]+Fraction(Y[1])
                else:
                     matrix1[:,1]=matrix1[:,1]

        if ops[1].replace('-','')[0]=='x':
            if Y[0]=='-':
                if Y[1]!='':
                    matrix1[:,1]=-matrix1[:,0]+Fraction(Y[1])
                else:
                    matrix1[:,1]=-matrix1[:,0]
            else:
                if Y[1]!='':    
                    matrix1[:,1]=matrix1[:,0]+Fraction(Y[1])
                else:
                     matrix1[:,1]=matrix1[:,0]

        if ops[1].replace('-','')[0]=='z':
            if Y[0]=='-':
                if Y[1]!='':
                    matrix1[:,1]=-matrix1[:,2]+Fraction(Y[1])
                else:
                    matrix1[:,1]=-matrix1[:,2]
            else:
                if Y[1]!='':    
                    matrix1[:,1]=matrix1[:,2]+Fraction(Y[1])
                else:
                     matrix1[:,1]=matrix1[:,2]

        MM[:,1]=matrix1[:,1]
        matrix1=Matrix.copy()
       
        if ops[2].replace('-','')[0]=='z':
            if Z[0]=='-':
                if Z[1]!='':
                    matrix1[:,2]=-matrix1[:,2]+Fraction(Z[1])
                else:
                     matrix1[:,2]=-matrix1[:,2]
            else:
                if Z[1]!='':
                    matrix1[:,2]=matrix1[:,2]+Fraction(Z[1])
                else:
                    matrix1[:,2]=matrix1[:,2]

        if ops[2].replace('-','')[0]=='x':
            if Z[0]=='-':
                if Z[1]!='':
                    matrix1[:,2]=-matrix1[:,0]+Fraction(Z[1])
                else:
                     matrix1[:,2]=-matrix1[:,0]
            else:
                if Z[1]!='':
                    matrix1[:,2]=matrix1[:,0]+Fraction(Z[1])
                else:
                    matrix1[:,2]=matrix1[:,0]

        if ops[2].replace('-','')[0]=='y':
            if Z[0]=='-':
                if Z[1]!='':
                    matrix1[:,2]=-matrix1[:,1]+Fraction(Z[1])
                else:
                     matrix1[:,2]=-matrix1[:,1]
            else:
                if Z[1]!='':
                    matrix1[:,2]=matrix1[:,1]+Fraction(Z[1])
                else:
                    matrix1[:,2]=matrix1[:,1]

        MM[:,2]=matrix1[:,2]
        MT.append(MM.round(4).tolist())
    return MT


def appendseg(initialcoor,latercoor,trancoor):
    import numpy as np
  
    initialvec1=np.array(initialcoor[1])-np.array(initialcoor[0])+np.array(initialcoor[2])-np.array(initialcoor[0])
    initialvec2=np.cross(np.array(initialcoor[1])-np.array(initialcoor[0]),np.array(initialcoor[2])-np.array(initialcoor[0]))\
                /np.linalg.norm(np.cross(np.array(initialcoor[1])-np.array(initialcoor[0]),np.array(initialcoor[2])-np.array(initialcoor[0])))\
                *(np.linalg.norm(np.array(initialcoor[1])-np.array(initialcoor[0]))+
                 np.linalg.norm( np.array(initialcoor[2])-np.array(initialcoor[0])))
    initialvec3=np.cross(initialvec1,initialvec2)/np.linalg.norm(np.cross(initialvec1,initialvec2))\
    *(np.linalg.norm(initialvec1)+np.linalg.norm(initialvec2))
    initialm=np.matrix(np.vstack((initialvec1,initialvec2,initialvec3))).T

    latervec1=np.array(latercoor[1])-np.array(latercoor[0])+np.array(latercoor[2])-np.array(latercoor[0])
    latervec2=np.cross(np.array(latercoor[1])-np.array(latercoor[0]),np.array(latercoor[2])-np.array(latercoor[0]))\
    /np.linalg.norm(np.cross(np.array(latercoor[1])-np.array(latercoor[0]),np.array(latercoor[2])-np.array(latercoor[0])))\
    *(np.linalg.norm(np.array(latercoor[1])-np.array(latercoor[0]))+np.linalg.norm(np.array(latercoor[2])-np.array(latercoor[0])))
    latervec3=np.cross(latervec1,latervec2)/np.linalg.norm(np.cross(latervec1,latervec2))\
    *(np.linalg.norm(latervec1)+np.linalg.norm(latervec2))
    laterm=np.matrix(np.vstack((latervec1,latervec2,latervec3))).T

    tranm=[np.array(trancoor[i])-np.array(initialcoor[0]) for i in range(len(trancoor))]
    tranm=np.matrix(tranm).T

    m=np.linalg.inv(initialm).dot(tranm)
   
    new=np.array((laterm.dot(m)).T)
    newnew=[(new[i]+np.array(latercoor[0])).tolist() for i in range(len(trancoor))]
    return  newnew

def reflect_frage(new_xyz_c_mat,original_linker):
    linker_in_mof=np.array(new_xyz_c_mat)
    old_linker=np.array(original_linker)           
    lin=0
    vectors=0
    while vectors==0:
        one_length=np.around((linker_in_mof[lin,:]-linker_in_mof[lin+1,:]).dot(linker_in_mof[lin,:]-linker_in_mof[lin+1,:])-\
        (old_linker[lin,:]-old_linker[lin+1,:]).dot(old_linker[lin,:]-old_linker[lin+1,:]),1)
        two_length=np.around((linker_in_mof[lin+2,:]-linker_in_mof[lin+1,:]).dot(linker_in_mof[lin+2,:]-linker_in_mof[lin+1,:])-\
        (old_linker[lin+2,:]-old_linker[lin+1,:]).dot(old_linker[lin+2,:]-old_linker[lin+1,:]),1)
        if one_length==0.0 and two_length==0.0:
            vectors=1
            initialcoor=[original_linker[lin,:],original_linker[lin+1,:],original_linker[lin+2,:]]
            latercoor=[new_xyz_c_mat[lin,:],new_xyz_c_mat[lin+1,:],new_xyz_c_mat[lin+2,:]]
            trancoor=original_linker
            new_xyz_c_mat_maped=np.array(appendseg(initialcoor,latercoor,trancoor))
        else:
            lin=lin+1
            if lin>linker_in_mof.shape[0]-3:
                break
    return new_xyz_c_mat_maped

def rotation(Angle,new_xyz_c_mat0,center,VV):
    VV=np.array(VV)
    new_xyz_c_mat=new_xyz_c_mat0.T
    a,b,c=center[0],center[1],center[2]
    length=(VV.dot(VV))**0.5
    u,v,w=VV[0]/length,VV[1]/length,VV[2]/length
    angle=math.radians(Angle)

    matrix=np.matrix([[u*u+(v*v+w*w)*math.cos(angle),u*v*(1-math.cos(angle))-w*math.sin(angle),u*w*(1-math.cos(angle))+v*math.sin(angle),(a*(v*v+w*w)-u*(b*v+c*w))*(1-math.cos(angle))+(b*w-c*v)*math.sin(angle)],\
[u*v*(1-math.cos(angle))+w*math.sin(angle),v*v+(u*u+w*w)*math.cos(angle),v*w*(1-math.cos(angle))-u*math.sin(angle),(b*(u*u+w*w)-v*(a*u+c*w))*(1-math.cos(angle))+(c*u-a*w)*math.sin(angle)],\
[u*w*(1-math.cos(angle))-v*math.sin(angle),v*w*(1-math.cos(angle))+u*math.sin(angle),w*w+(u*u+v*v)*math.cos(angle),(c*(u*u+v*v)-w*(a*u+b*v))*(1-math.cos(angle))+(a*v-b*u)*math.sin(angle)],\
[0,0,0,1]])
    add_line=np.ones(new_xyz_c_mat.shape[1])
    result=matrix.dot(np.array(np.insert(new_xyz_c_mat,3,values=add_line,axis=0)))
    res=result.round(6)
    return res.T[:,0:3]
def get_xyz_c(new_fractonal_c,a,b,c,A1,B1,C1):
    A=math.radians(A1)
    B=math.radians(B1)
    C=math.radians(C1)
    matrix=np.matrix([[a,b*math.cos(C),c*math.cos(B)],[0,b*math.sin(C),c*(math.cos(A)-math.cos(B)*math.cos(C))/math.sin(C)],[0,0,c*np.sqrt(math.sin(B)*math.sin(B)-((math.cos(A)-math.cos(B)*math.cos(C))/math.sin(C))**2)]])
    xyz=np.matrix(new_fractonal_c)
    m=xyz.T
    ff=(matrix.dot(m)).T.round(4)
    return ff

def get_frac_c(new_fractonal_c,a,b,c,A1,B1,C1):
    A=math.radians(A1)
    B=math.radians(B1)
    C=math.radians(C1)
    matrix=np.matrix([[a,b*math.cos(C),c*math.cos(B)],[0,b*math.sin(C),c*(math.cos(A)-math.cos(B)*math.cos(C))/math.sin(C)],[0,0,c*np.sqrt(math.sin(B)*math.sin(B)-((math.cos(A)-math.cos(B)*math.cos(C))/math.sin(C))**2)]])
    M=np.linalg.inv(matrix)
    xyz=np.matrix(new_fractonal_c)
    m=xyz.T
    fraction=M.dot((m))
    Frac=fraction.T.round(4)
    return Frac

def check_distance_of_cell(CIFfile,unitcell):
    a,b,c,A1,B1,C1=unitcell[0],unitcell[1],unitcell[2],unitcell[3],unitcell[4],unitcell[5]
    data0=np.loadtxt(CIFfile+'.cif',skiprows=19,usecols=[2,3,4])
    data=get_xyz_c(data0,a,b,c,A1,B1,C1)
    error=0
    for i in range(data.shape[0]):
        for j in range(i):
            VV=data[i,:]-data[j,:]
            distance=(VV.dot(VV))**0.5
            if 0<distance<0.8:              
                error=1
                break
        break
    return error

def check_coll(CIFfile,unitcell):
    a,b,c,A1,B1,C1=unitcell[0],unitcell[1],unitcell[2],unitcell[3],unitcell[4],unitcell[5]
    data=np.loadtxt(CIFfile+'.cif',skiprows=18,usecols=[2,3,4])
    #data=get_xyz_c(data0,a,b,c,A1,B1,C1)
    with open(CIFfile+'.cif') as cif_data:
        cif_datas=cif_data.readlines()
    Label=[i.split()[0] for i in cif_datas[18:]]
    C_N=[]
    wrong_bonds=[]
    for i in range(0,180):
        for j in range(180,194):
            VV=np.array(data)[i,:]-np.array(data)[j,:]
            for ii in range(3):
                 while abs(VV[ii])>0.5:
                    if VV[ii]<0:
                        VV[ii]=VV[ii]+1
                    else:
                        VV[ii]=VV[ii]-1
            distance=((abs(VV[0])*a)**2+(abs(VV[1])*b)**2+(abs(VV[2])*c)**2)**0.5

            if 0<distance<0.8: 
                wrong_bonds.append(1)
            else:
                num=bonder([Label[i],Label[j]],distance,[])
                if num==1 and [Label[i],Label[j]] !=['N','C']:
                    wrong_bonds.append(1)
                if num==1 and [Label[i],Label[j]] ==['N','C']:
                    C_N.append(1) 
    return sum(wrong_bonds)+abs(sum(C_N)-2)

def cal_distance_to_determined_C(all_N_atoms):

    all_distance=[]
    index=[] 
    all_connections=[] 
    for nn in range(len(all_N_atoms[0])):
        for i in range(1,len(all_N_atoms)):
            for j in range(len(all_N_atoms[i])):
                VV=np.array(all_N_atoms[0][nn])-np.array(all_N_atoms[i][j])
                length=abs((VV.dot(VV))**0.5)
                if 14.5<length<16.7:
                    all_connections.append([all_N_atoms[0][nn],all_N_atoms[i][j]])
                    index.append(str(nn)+str(i)+str(j))######### first number is the N number of center nodes,second is the number of nodes ,last is the N number of node
                    all_distance.append(length)

    return all_connections,index

#################################################################################
def write_CIF(out_name,matrixs,Label,cell_value):
    if 1>0:
        with open(out_name+'.cif','w') as W:
            W.write('_cell_length_a    '+str(cell_value[0])+'\n'+'_cell_length_b    '+str(cell_value[1])+'\n'+'_cell_length_c    '+str(cell_value[2])+'\n')
            W.write('_cell_angle_alpha    '+str(cell_value[3])+'\n'+'_cell_angle_beta    '+str(cell_value[4])+'\n'+'_cell_angle_gamma    '+str(cell_value[5])+'\n')
            W.write("_symmetry_space_group_name_H-M		'P1'"+'\n'+'_symmetry_Int_Tables_number		1'+'\n'+'_symmetry_cell_setting		Monoclinic'+'\n')
            W.write('loop_'+'\n'+'_symmetry_equiv_pos_as_xyz'+'\n'+"'+x,+y,+z'"+'\n'+'loop_'+'\n'+'_atom_site_label'+'\n'+'_atom_site_type_symbol'+'\n'+'_atom_site_fract_x'+'\n'+'_atom_site_fract_y'+'\n'+'_atom_site_fract_z'+'\n')
            for matrix in matrixs[0:4]:
                for i in range(len(Label)):  
                    W.write(Label[i]+' '+Label[i]+' '+str(matrix[i][0])+' '+str(matrix[i][1])+' '+str(matrix[i][2])+'\n')

def write_linker_CIF(out_name,matrixs,Label,cell_value):
    if 1>0:
        with open(out_name+'.cif','a') as W:
            for matrix in matrixs[0:8]:
                for i in range(len(Label)):  
                    W.write(Label[i]+' '+Label[i]+' '+str(matrix[i][0])+' '+str(matrix[i][1])+' '+str(matrix[i][2])+'\n')

def cal_fractional_min_dis(pointsA,C_coordinates,a,b,c):#####  node target atoms, all_atoms list,return minx istance
    dis_alls=[]
    for nnn1 in C_coordinates[0:]:
        VV=np.array(nnn1)-np.array(pointsA)
        for i in range(3):
            while abs(VV[i])>0.5:
                if VV[i]<0:
                    VV[i]=VV[i]+1
                else:
                    VV[i]=VV[i]-1
        length=((abs(VV[0])*a)**2+(abs(VV[1])*b)**2+(abs(VV[2])*c)**2)**0.5
        dis_alls.append(abs(length))
    AA_dis=sorted(dis_alls)
    NNNN_dis=min(dis_alls)
    return NNNN_dis

def check_N_N_dis(N_atoms, unitcell):
   # print (N_atoms)
    a,b,c=unitcell[0],unitcell[1],unitcell[2]
    wrong_NN=0
    all_N_N_dis=[]
    for i in range(4):
        min_dis_NN=cal_fractional_min_dis(N_atoms[i], N_atoms[4:16], a, b, c)
       # print (i,min_dis_NN,66666666666,N_atoms[i])
        all_N_N_dis.append(min_dis_NN)
    if min(all_N_N_dis)<2:
        wrong_NN=1
    return wrong_NN

def cal_fractional_lin_dis(pointsA,C_coordinates,a,b,c):#####  linker target atoms, all_atoms list,return minx istance
    dis_alls=[]
    all_vectors=[]
    for nnn1 in C_coordinates[0:]:
        VV=np.array(nnn1)-np.array(pointsA)
        for i in range(3):
            while abs(VV[i])>0.5:
                if VV[i]<0:
                    VV[i]=VV[i]+1
                else:
                    VV[i]=VV[i]-1
        length=((abs(VV[0])*a)**2+(abs(VV[1])*b)**2+(abs(VV[2])*c)**2)**0.5
        dis_alls.append(abs(length))
        all_vectors.append([VV[0]*a,VV[1]*b,VV[2]*c])
    NNNN_dis=min(dis_alls)
    
    return NNNN_dis,all_vectors[dis_alls.index(NNNN_dis)],C_coordinates[dis_alls.index(NNNN_dis)] ##### return min C-N, vector C-N,and the N coordinates


def check_node_overate(all_atoms_frac,unitcell):#############3 node
    over=0
    num_of_symmetry,num_of_node=16,45
    a,b,c,A1,B1,C1=unitcell[0],unitcell[1],unitcell[2],unitcell[3],unitcell[4],unitcell[5]
    all_N_atoms=[]
    for i in range(num_of_symmetry):
        for j in range(1,5):
            all_N_atoms.append(all_atoms_frac[i][j])  

    all_C_N_atoms=[]
    for i in range(num_of_symmetry):
        for j in range(5,9):
            all_C_N_atoms.append(all_atoms_frac[i][j])  

    four_min_dis=[] 
    for i in range(4):
        dis=cal_fractional_min_dis(all_N_atoms[i],all_N_atoms[4:],a,b,c)  
        four_min_dis.append(dis)
    if len(four_min_dis)>0:
        if sum(four_min_dis)/len(four_min_dis)<0.5:
            over=1
        #print (four_min_dis)
       # print ('all nitrogen dis')
    return over,all_N_atoms,sum(four_min_dis)/len(four_min_dis),all_C_N_atoms

def calculated_angle(V1,V2):

    L1=(V1.dot(V1))**0.5
    L2=(V2.dot(V2))**0.5
    #print L1,L2
    cos_angle=V1.dot(V2)/(L1*L2)
    if cos_angle>1:
        cos_angle=1    
    angle_=np.arccos(cos_angle)
    angle=angle_*360/2/np.pi
   # print angle
    return angle
def cal_C_H_vector(a,b,c,A1,B1,C1,C_atoms,H_atoms):
    H_C=[]
    unit=[a,b,c]
    VV=np.array(H_atoms)-np.array(C_atoms)
    for i in range(3):
        while abs(VV[i])>0.5:
            if VV[i]<0:
                VV[i]=VV[i]+1
            else:
                VV[i]=VV[i]-1
    H_C=[VV[0]*a,VV[1]*b,VV[2]*c]
    #print (H_C,4444)
    return np.array(H_C)


def check_linker_bonds(all_C_atoms,all_N_atoms,unitcell,all_H_atoms,all_C_N_atoms):###############linker
    bond_write=0
    a,b,c,A1,B1,C1=unitcell[0],unitcell[1],unitcell[2],unitcell[3],unitcell[4],unitcell[5]

    four_min_dis=[] 
    
    write=[]
    all_angles=[]
    all_angles_nodes=[]
    for i in range(len(all_C_atoms)):
        all_dis=[]
        dis,N_C_vector,N_coordinates=cal_fractional_lin_dis(all_C_atoms[i],all_N_atoms,a,b,c) 

        dis2,C_C_vector,C_N_coordinates=cal_fractional_lin_dis(N_coordinates,all_C_N_atoms,a,b,c)########## node C-C vector

        H_C=cal_C_H_vector(a,b,c,A1,B1,C1,all_C_atoms[i],all_H_atoms[i])
        H_C_N=calculated_angle(np.array(N_C_vector),np.array(H_C)) ###########3 angle-linker
        C_C_N=180-calculated_angle(np.array(N_C_vector),np.array(C_C_vector)) ########### angle-node


        all_angles.append(H_C_N)
        all_angles_nodes.append(C_C_N)
        all_dis.append(dis) 

        write.append(abs(dis-1.273))
        four_min_dis.append(abs(dis-1.273))
    if max(four_min_dis)<=0.15:
        bond_write=1   
    
    angle_error=abs(sum(all_angles)/len(all_angles)-122.439)/122.439
    angle_error2=abs(sum(all_angles_nodes)/len(all_angles_nodes)-120.414)/120.414
    
    return bond_write,sum(four_min_dis)/len(four_min_dis),angle_error,angle_error2


def get_node(XYZ,angles,unitcell,first_atom):########## node xyz coordinates, rotated angle, one node fractional coordinates
    r1,r2,r3=angles[0],angles[1],angles[2]
    center=[np.sum(XYZ[:,i])/XYZ.shape[0] for i in range(XYZ.shape[1])] #rotation
    #r1,r2,r3=-82.9,-59.5,-58.5
    XYZ=rotation(r1,XYZ,center,[1,0,0])
    XYZ=rotation(r2,XYZ,center,[0,1,0])
    XYZ=rotation(r3,XYZ,center,[0,0,1])
    a,b,c,A1,B1,C1=unitcell[0],unitcell[1],unitcell[2],unitcell[3],unitcell[4],unitcell[5]
    matrix=get_frac_c(XYZ,a,b,c,A1,B1,C1)
    center=[np.sum(matrix[:,i])/matrix.shape[0] for i in range(matrix.shape[1])] #rotation
    #np.savetxt(name+'one_.cif',matrix,'%.4f')
    x0,y0,z0=first_atom[0],first_atom[1],first_atom[2]
    mx=x0-matrix[0,0]
    my=y0-matrix[0,1]
    mz=z0-matrix[0,2] 
    #mx,my,mz=0,0,0
    matrix[:,0]=matrix[:,0]+mx ######### move
    matrix[:,1]=matrix[:,1]+my ######### move
    matrix[:,2]=matrix[:,2]+mz ######### move

    MT=operation(matrix,[' x, y, z','-x+1/2, -y, z+1/2','-y+3/4, x+1/4, z+1/4','y+3/4, -x+3/4, z+3/4','x+1/2, y+1/2, z+1/2','-x+1, -y+1/2, z+1','-y+5/4, x+3/4, z+3/4','y+5/4, -x+5/4, z+5/4','-x, -y, -z','x-1/2, y, -z-1/2','y-3/4, -x-1/4, -z-1/4','-y-3/4, x-3/4, -z-3/4','-x+1/2, -y+1/2, -z+1/2','x, y+1/2, -z','y-1/4, -x+1/4, -z+1/4','-y-1/4, x-1/4, -z-1/4'])############## get all node coordinates

    return MT

def get_linker(XYZ,angles,unitcell,first_atom):##########linker xyz coordinates, rotated angle, one node fractional coordinates
    r1,r2,r3=angles[0],angles[1],angles[2]
    center=[np.sum(XYZ[:,i])/XYZ.shape[0] for i in range(XYZ.shape[1])] #rotation
    #r1,r2,r3=-82.9,-59.5,-58.5
    XYZ=rotation(r1,XYZ,center,[1,0,0])
    XYZ=rotation(r2,XYZ,center,[0,1,0])
    XYZ=rotation(r3,XYZ,center,[0,0,1])
    a,b,c,A1,B1,C1=unitcell[0],unitcell[1],unitcell[2],unitcell[3],unitcell[4],unitcell[5]
    matrix=get_frac_c(XYZ,a,b,c,A1,B1,C1)
    center=[np.sum(matrix[:,i])/matrix.shape[0] for i in range(matrix.shape[1])] #rotation
    x0,y0,z0=first_atom[0],first_atom[1],first_atom[2]
    mx=x0-center[0]
    my=y0-center[1]
    mz=z0-center[2]
    #mx,my,mz=0,0,0
    matrix[:,0]=matrix[:,0]+mx ######### move
    matrix[:,1]=matrix[:,1]+my ######### move
    matrix[:,2]=matrix[:,2]+mz ######### move

    MT=operation(matrix,[' x, y, z','-x+1/2, -y, z+1/2','-y+3/4, x+1/4, z+1/4','y+3/4, -x+3/4, z+3/4','x+1/2, y+1/2, z+1/2','-x+1, -y+1/2, z+1','-y+5/4, x+3/4, z+3/4','y+5/4, -x+5/4, z+5/4','-x, -y, -z','x-1/2, y, -z-1/2','y-3/4, -x-1/4, -z-1/4','-y-3/4, x-3/4, -z-3/4','-x+1/2, -y+1/2, -z+1/2','x, y+1/2, -z','y-1/4, -x+1/4, -z+1/4','-y-1/4, x-1/4, -z-1/4'])############## get all node coordinates

    return MT

def constructed_CIF(parameters,node_name,linker_name,NName):

    angles=[parameters[0],parameters[1],parameters[2]]
    angles_linker=[parameters[3],parameters[4],parameters[5]]

    unitcell=[20.150,20.150,8.850,90,90.0,90.0] ########3 unitcells
    node_xyz=np.loadtxt(node_name,skiprows=2,usecols=[1,2,3])  ######### node file
    with open(node_name) as Xs:
        Xss=Xs.readlines()[2:]
    node_label=[i.split()[0] for i in Xss]   

    MT=get_node(node_xyz,angles,unitcell,[0.5,0.75,0.125])
    overate,all_N_atoms,ave_dis,all_C_N_atoms=check_node_overate(MT,unitcell) ########### if overate,all N,averdis
    #print (overate,ave_dis)
    out_name='./CIF/'+NName
    write_CIF(out_name,MT,node_label,unitcell)
    wrong_NN=0
    wrong_NN=check_N_N_dis(all_N_atoms,unitcell)
    with open('all_records.txt','a') as Rt:
        Rt.write(out_name+' '+str(wrong_NN)+' '+str(overate)+'\n')

    if overate==1 and wrong_NN==0:         ######################### calculated
        linker_xyz=np.loadtxt(linker_name,skiprows=2,usecols=[1,2,3])  ######### linker file
        with open(linker_name) as Xs:
            Xss=Xs.readlines()[2:]
        linker_label=[i.split()[0] for i in Xss]    
        
        MT_linker=get_linker(linker_xyz,angles_linker,unitcell,[0.5,0,0])
        all_C_atoms=MT_linker[0][0:2]
        all_H_atoms=MT_linker[0][2:4]##########3 in fact this is C atoms
        bond_write,two_C,angles_error,angles_error2=check_linker_bonds(all_C_atoms,all_N_atoms,unitcell,all_H_atoms,all_C_N_atoms)########### if write C-N bonds, bond distance

        write_linker_CIF(out_name,MT_linker,linker_label,unitcell)
        rewards=1-1/(ave_dis+two_C+angles_error+angles_error2+1)
        check_bond=check_coll(out_name,unitcell)
        R_value=rewards+check_bond#+R/100############################

        with open('result.txt','a') as RR:
            RR.write(out_name[6:]+' '+str(round(ave_dis,2))+' '+str(round(two_C,2))+' '+str(round(angles_error,2))+' '+str(round(angles_error2,2))+' '+str(round(rewards,2))+' '+str(round(R_value,2))+'\n')
    else:
        R_value=50+ave_dis
        os.system('rm '+out_name+'.cif')
    return round(R_value,2)
       
import time

def demo_func(x1, x2, x3, x4,x5,x6):
    node_name='MS_draw_node.xyz'
    linker_name='MS_draw_linker.xyz'
    PP=[x1,x2,x3,x4,x5,x6]
    CIF_name=str(PP[0])+'_'+str(PP[1])+'_'+str(PP[2])+'_'+str(PP[3])+'_'+str(PP[4])+'_'+str(PP[5])
    R_value=constructed_CIF(PP,node_name,linker_name,CIF_name)
    return R_value
import time
from sko.PSO import PSO
tim1=time.time()
pso = PSO(func=demo_func, n_dim=6, pop=2000, max_iter=100, lb=[0,0,0, 0, 0,0], ub=[180,180,360, 180, 180,360])
tim2=time.time()
fitness = pso.run()
import matplotlib.pyplot as plt
plt.switch_backend('agg')
with open('inf.txt','a') as inf:
    inf.write(str(tim2-tim1)+'\n')
    for i in pso.gbest_y_hist:
        inf.write(str(i)+'\n')
print (pso.gbest_y_hist)
plt.plot(pso.gbest_y_hist)
plt.savefig('1.png')
plt.close()


