# -*- coding: utf-8 -*-

import numpy as np
import math
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error
import scipy as sp


Ligand_Atom = ['C','N','O','S','P','F','Cl','Br','I']

Year = '2013'
pre = '../data/' + Year + '/'

f1 = open(pre + 'name/train_data_' + Year + '.txt')
pre_train_data = f1.readlines()
train_data = eval(pre_train_data[0])
f1.close()

f1 = open(pre + 'name/test_data_' + Year + '.txt')
pre_test_data = f1.readlines()
test_data = eval(pre_test_data[0])
f1.close()

f1 = open(pre + 'name/all_data_' + Year + '.txt')
pre_all_data = f1.readlines()
all_data = eval(pre_all_data[0])
f1.close()

def get_index(a,b):
    t = len(b)
    if a=='Cl':
        return 6
    if a=='CL':
        return 6
    if a=='Br':
        return 7
    if a=='BR':
        return 7
    
    for i in range(t):
        if a[0]==b[i]:
            return i
    return -1


def extract_coordinate_of_ligand(start,end):
    #t1 = len(all_data)
    for i in range(start,end):
        print('process {0}-th '.format(i))
        name = all_data[i]
        
        ligand = {}
        for ii in range(9):
            ligand[Ligand_Atom[ii]] = []
            
        t2 = pre + 'refined/' + name + '/' + name + '_ligand.mol2'
        f2 = open(t2,'r')
        contents = f2.readlines()
        t3 = len(contents)
        start = 0
        end = 0
        for jj in range(t3):
            if contents[jj][0:13]=='@<TRIPOS>ATOM':
                start = jj + 1
                continue
            if contents[jj][0:13]=='@<TRIPOS>BOND':
                end = jj - 1
                break
        for kk in range(start,end+1):
            if contents[kk][8:17]=='thiophene':
                print('thiophene',kk)
            atom = contents[kk][47:49]
            atom = atom.strip()
            index = get_index(atom,Ligand_Atom)
            if index==-1:
                continue
            else:
                    
                ligand[Ligand_Atom[index]].append(contents[kk][17:46])
        f2.close()
        
        #print(ligand)
        # C
        t = len(ligand['C'])
        temp = np.zeros((t,3))
        for tt in range(t):
            temp[tt][0] = float(ligand['C'][tt][0:9])
            temp[tt][1] = float(ligand['C'][tt][9:19])
            temp[tt][2] = float(ligand['C'][tt][19:29])
        #print(temp)
        filename = pre + 'ligand_coordinate/' + name + '_C_coordinate.csv'
        np.savetxt(filename,temp,delimiter=',')
        
        # N
        t = len(ligand['N'])
        temp = np.zeros((t,3))
        for tt in range(t):
            temp[tt][0] = float(ligand['N'][tt][0:9])
            temp[tt][1] = float(ligand['N'][tt][9:19])
            temp[tt][2] = float(ligand['N'][tt][19:29])
        #print(temp)
        filename = pre + 'ligand_coordinate/' + name + '_N_coordinate.csv'
        np.savetxt(filename,temp,delimiter=',')
        
        # O
        t = len(ligand['O'])
        temp = np.zeros((t,3))
        for tt in range(t):
            temp[tt][0] = float(ligand['O'][tt][0:9])
            temp[tt][1] = float(ligand['O'][tt][9:19])
            temp[tt][2] = float(ligand['O'][tt][19:29])
        #print(temp)
        filename = pre + 'ligand_coordinate/' + name + '_O_coordinate.csv'
        np.savetxt(filename,temp,delimiter=',')
        
        # C,N
        t1 = len(ligand['C'])
        t2 = len(ligand['N'])
        temp = np.zeros((t1+t2,3))
        cc = 0
        for tt in range(t1):
            temp[cc][0] = float(ligand['C'][tt][0:9])
            temp[cc][1] = float(ligand['C'][tt][9:19])
            temp[cc][2] = float(ligand['C'][tt][19:29])
            cc = cc + 1
        for tt in range(t2):
            temp[cc][0] = float(ligand['N'][tt][0:9])
            temp[cc][1] = float(ligand['N'][tt][9:19])
            temp[cc][2] = float(ligand['N'][tt][19:29])
            cc = cc + 1
        #print(temp)
        filename = pre + 'ligand_coordinate/' + name + '_CN_coordinate.csv'
        np.savetxt(filename,temp,delimiter=',')
        
        # C,O
        t1 = len(ligand['C'])
        t2 = len(ligand['O'])
        temp = np.zeros((t1+t2,3))
        cc = 0
        for tt in range(t1):
            temp[cc][0] = float(ligand['C'][tt][0:9])
            temp[cc][1] = float(ligand['C'][tt][9:19])
            temp[cc][2] = float(ligand['C'][tt][19:29])
            cc = cc + 1
        for tt in range(t2):
            temp[cc][0] = float(ligand['O'][tt][0:9])
            temp[cc][1] = float(ligand['O'][tt][9:19])
            temp[cc][2] = float(ligand['O'][tt][19:29])
            cc = cc + 1
        #print(temp)
        filename = pre + 'ligand_coordinate/' + name + '_CO_coordinate.csv'
        np.savetxt(filename,temp,delimiter=',')
        
        # C,N,O
        t1 = len(ligand['C'])
        t2 = len(ligand['N'])
        t3 = len(ligand['O'])
        temp = np.zeros((t1+t2+t3,3))
        cc = 0
        for tt in range(t1):
            temp[cc][0] = float(ligand['C'][tt][0:9])
            temp[cc][1] = float(ligand['C'][tt][9:19])
            temp[cc][2] = float(ligand['C'][tt][19:29])
            cc = cc + 1
        for tt in range(t2):
            temp[cc][0] = float(ligand['N'][tt][0:9])
            temp[cc][1] = float(ligand['N'][tt][9:19])
            temp[cc][2] = float(ligand['N'][tt][19:29])
            cc = cc + 1
        for tt in range(t3):
            temp[cc][0] = float(ligand['O'][tt][0:9])
            temp[cc][1] = float(ligand['O'][tt][9:19])
            temp[cc][2] = float(ligand['O'][tt][19:29])
            cc = cc + 1
        
        filename = pre + 'ligand_coordinate/' + name + '_CNO_coordinate.csv'
        np.savetxt(filename,temp,delimiter=',')
        
        # C,N,O,S,P,F,Cl,Br,I
        t = 0
        for tt in range(9):
            t = t + len(ligand[Ligand_Atom[tt]])
        temp = np.zeros((t,3))
        cc = 0
        for tt in range(9):
            for ttt in range(len(ligand[Ligand_Atom[tt]])):
                temp[cc][0] = float(ligand[Ligand_Atom[tt]][ttt][0:9])
                temp[cc][1] = float(ligand[Ligand_Atom[tt]][ttt][9:19])
                temp[cc][2] = float(ligand[Ligand_Atom[tt]][ttt][19:29])
                cc = cc + 1
        #print(temp)
        filename = pre + 'ligand_coordinate/' + name + '_all_coordinate.csv'
        np.savetxt(filename,temp,delimiter=',')


def distance_of_two_points(p1,p2):
    temp = pow(p1[0]-p2[0],2) + pow(p1[1]-p2[1],2) + pow(p1[2]-p2[2],2)
    res = pow(temp,0.5)
    return res

def get_complex(N,filtration,typ):
    name = all_data[N]
    filename = pre + 'ligand_coordinate/' + name + '_' + typ + '_coordinate.csv'
    point = np.loadtxt(filename,delimiter=',')
    t = point.shape
    if len(t)==1:
        return [ [0,0,0,0] ]
    #  create the distance list, then sort
    dis_list = []
    for i in range(t[0]):
        for j in range(i+1,t[0]):
            dis = distance_of_two_points(point[i], point[j])
            if dis<=filtration:
                dis_list.append([ i,j,dis ])
    dis_list = sorted(dis_list,key=lambda x:(x[2]))
    
    # create filtered neighbourhood complex
    simplices = []
    count = 0
    neighbor = []
    is_edge_matrix = np.zeros((t[0],t[0]))
    is_triangle_matrix = np.zeros((t[0],t[0],t[0]))
    
    # set neighbor as empty
    for i in range(t[0]):
        neighbor.append([])
        
    # add all points as 0-simplices
    for i in range(t[0]):
        temp = [ count, 0, 0, i ]
        simplices.append(temp)
        count = count + 1
    
    # add higher dimensional simplices
    for i in range(len(dis_list)):
        one = dis_list[i][0]
        two = dis_list[i][1]
        dis = dis_list[i][2]
        
        # one
        number_had = len(neighbor[one])
        if number_had==0:
            neighbor[one].append(two)
        else:
            
            # edge
            for ii in range(len(neighbor[one])):
                one2 = neighbor[one][ii]
                m = min(two,one2)
                M = max(two,one2)
                if is_edge_matrix[m][M]==0:
                    is_edge_matrix[m][M] = 1
                    temp = [count, dis, 1, m, M]
                    simplices.append(temp)
                    count = count + 1
            
            # triangle
            for ii in range(len(neighbor[one])):
                one2 = neighbor[one][ii]
                for jj in range(ii+1,len(neighbor[one])):
                    one3 =  neighbor[one][jj]
                    m = min(two,one2,one3)
                    M = max(two,one2,one3)
                    mid = two + one2 + one3 - m - M
                    if is_triangle_matrix[m][mid][M]==0:
                        is_triangle_matrix[m][mid][M] = 1
                        temp = [count, dis, 2, m, mid, M]
                        simplices.append(temp)
                        count = count + 1
            
            neighbor[one].append(two)
            
        
        # two
        number_had = len(neighbor[two])
        if number_had==0:
            neighbor[two].append(one)
        else:
            
            # edge
            for ii in range(len(neighbor[two])):
                one2 = neighbor[two][ii]
                m = min(one,one2)
                M = max(one,one2)
                if is_edge_matrix[m][M]==0:
                    is_edge_matrix[m][M] = 1
                    temp = [count, dis, 1, m, M]
                    simplices.append(temp)
                    count = count + 1
            
            # triangle
            for ii in range(len(neighbor[two])):
                one2 = neighbor[two][ii]
                for jj in range(ii+1,len(neighbor[two])):
                    one3 =  neighbor[two][jj]
                    m = min(one,one2,one3)
                    M = max(one,one2,one3)
                    mid = one + one2 + one3 - m - M
                    if is_triangle_matrix[m][mid][M]==0:
                        is_triangle_matrix[m][mid][M] = 1
                        temp = [count, dis, 2, m, mid, M]
                        simplices.append(temp)
                        count = count + 1
            
            neighbor[two].append(one)
    return simplices




def get_point_index(point,points):
    for i in range(len(points)):
        if point==points[i]:
            return i

def get_edge_index(p1,p2,edges):
    for i in range(len(edges)):
        if (p1==edges[i][0])&(p2==edges[i][1]):
            return i
    

def eigenvalue0_to_file(typ,simplices,name,filtration,grid):
    #print('process {0}-{1} combination of {2}'.format(P,L,name))
    filename = pre + 'ligand_eigenvalue_' + str(filtration) + '_zero/'
    
    
    if len(simplices)==0:
        # no complex, use -1 in the first position as a signal 
        filename1 = filename + name + '_' + typ + '_eigenvalue_0D.txt'
        res = [-1]
        f = open(filename1,'w')
        f.writelines(str(res))
        f.close()
        
        return
        
    
    #get 0-dimension laplacian
    
    number0 = int(filtration/grid)
    
    eigenvalue0 = [1] # have complex, use 1 in the first position as a signal
    for i in range(number0 + 1):
        # get eigenvalue for each filtra0 value with grid from 2 to filtration
        filtra0 = i * grid
        
        points = []
        edges = []
        for r in range(len(simplices)):
            if simplices[r][1]<=filtra0:
                if simplices[r][2]==0:
                    points.append(simplices[r][3])
                elif simplices[r][2]==1:
                    edges.append([ simplices[r][3] , simplices[r][4] ])
                
            else:
                break
        
        row = len(points)
        column = len(edges)
        
        if column==0:
            # only have points, no edges
            res = []
            for ii in range(row):
                res.append(0)
            eigenvalue0.append(res)
            
        else:
            zero_boundary = np.zeros((row,column))
            for j in range(column):
                one = edges[j][0]
                two = edges[j][1]
                index1 = get_point_index(one, points)
                index2 = get_point_index(two, points)
                zero_boundary[index1][j] = -1
                zero_boundary[index2][j] = 1
            Laplacian = np.dot( zero_boundary, zero_boundary.T )
            values = np.linalg.eigvalsh(Laplacian)
            res = []
            for iii in range(len(values)):
                res.append( values[iii] )
            eigenvalue0.append(res)
            #print(eigenvalue0)
    
    
    filename1 = filename + name + '_' + typ + '_eigenvalue_0D.txt'
    f = open(filename1,'w')
    f.writelines(str(eigenvalue0))
    f.close()


def eigenvalue1_to_file(typ,simplices,name,filtration,grid):
    #print('process {0}-{1} combination of {2}'.format(P,L,name))
    filename = pre + 'ligand_eigenvalue_' + str(filtration) + '_'  + 'one/'
    
    
    if len(simplices)==0:
        # no complex, use -1 in the first position as a signal 
        filename1 = filename + name + typ + '_eigenvalue_1D.txt'
        res = [-1]
        f = open(filename1,'w')
        f.writelines(str(res))
        f.close()
        return
        
    
    #get 1-dimension laplacian
    number1 = int(filtration/grid)
    eigenvalue1 = [1] # have complex, use 1 in the first position as a signal
    for i in range(number1 + 1):
        # get eigenvalue for each filtra0 value with grid from 2 to filtration
        filtra0 = i * grid
        points = []
        edges = []
        triangles = []
        for r in range(len(simplices)):
            if simplices[r][1]<=filtra0:
                if simplices[r][2]==0:
                    points.append(simplices[r][3])
                elif simplices[r][2]==1:
                    edges.append([ simplices[r][3] , simplices[r][4] ])
                elif simplices[r][2]==2:
                    triangles.append([ simplices[r][3], simplices[r][4], simplices[r][5] ])
                    
                
            else:
                break
        
        N0 = len(points)
        N1 = len(edges)
        N2 = len(triangles)
        
        
        if N1==0:
            # only have points, no edges
            res = [-1]
            eigenvalue1.append(res)
        elif (N1>0)&(N2==0):
            one_boundary = np.zeros((N0,N1))
            for j in range(N1):
                one = edges[j][0]
                two = edges[j][1]
                index1 = get_point_index(one, points)
                index2 = get_point_index(two, points)
                one_boundary[index1][j] = -1
                one_boundary[index2][j] = 1
            Laplacian = np.dot( one_boundary.T, one_boundary )
            values = np.linalg.eigvalsh(Laplacian)
            res = []
            for iii in range(len(values)):
                res.append( values[iii] )
            eigenvalue1.append(res)
        elif (N1>0)&(N2>0):
            
            one_boundary = np.zeros((N0,N1))
            for j in range(N1):
                one = edges[j][0]
                two = edges[j][1]
                index1 = get_point_index(one, points)
                index2 = get_point_index(two, points)
                one_boundary[index1][j] = -1
                one_boundary[index2][j] = 1
            L1 = np.dot( one_boundary.T, one_boundary )
            
            
            two_boundary = np.zeros((N1,N2))
            for j in range(N2):
                one = triangles[j][0]
                two = triangles[j][1]
                three = triangles[j][2]
                index1 = get_edge_index(one, two, edges)
                index2 = get_edge_index(one, three, edges)
                index3 = get_edge_index(two, three, edges)
                two_boundary[index1][j] = 1
                two_boundary[index2][j] = -1
                two_boundary[index3][j] = 1
                
            L2 = np.dot( two_boundary, two_boundary.T)
            Laplacian = L1 + L2
            values = np.linalg.eigvalsh(Laplacian)
            res = []
            for iii in range(len(values)):
                res.append( values[iii] )
            eigenvalue1.append(res)
    
    
    filename1 = filename + name + typ + '_eigenvalue_1D.txt'
    f = open(filename1,'w')
    f.writelines(str(eigenvalue1))
    f.close()





def eigenvalue_to_file(start,end,filtration0,filtration1,grid):
    
    for i in range(start,end):
        print(i)
        name = all_data[i]
        for typ in ['C','N','O','CN','CO','CNO','all']:
            simplices = get_complex(i,filtration0,typ)
            eigenvalue0_to_file(typ,simplices,name,filtration0,grid)
            eigenvalue1_to_file(typ,simplices,name,filtration1,grid)
        
def train_feature_to_file(start,end,filtration0,filtration1,grid):
    row = end - start
    N = 10
    number0 = int (filtration0/grid )
    number1 = int (filtration1/grid)
    add_value = int(10 * grid)
    column = 7 * number0 * N + 7*number1*N
    
    feature_matrix = np.zeros((row,column))
    #pre0 = 'E:\\singapore\\pdb_bind\\' + Year + '\\' + 'max_cocycle_bar\\eigenvalue_' + str(cutoff) + '_' + str(filtration0) + '_' +  'zero\\' 
    #pre = 'E:\\singapore\\pdb_bind\\' + Year + '\\' + 'max_cocycle_bar\\eigenvalue_' + str(cutoff) + '_' + str(filtration0) + '_' + str(filtration1) + '\\' 
    pre_zero = pre + 'ligand_eigenvalue_' + str(filtration0) + '_zero/'
    pre_one = pre + 'ligand_eigenvalue_' + str(filtration1) + '_one/'
    
    
    for i in range(start,end):
        print(i)
        name = train_data[i]
        count = 0
        for typ in ['C','N','O','CN','CO','CNO','all']:
            filename0 = pre_zero + name + '_' + typ + '_' + 'eigenvalue_0D.txt'
            f0 = open(filename0)
            pre_eigenvalue0 = f0.readlines()
            eigenvalue0 = eval(pre_eigenvalue0[0])
            f0.close()
            
            filename1 = pre_one + name + typ + '_' + 'eigenvalue_1D.txt'
            f1 = open(filename1)
            pre_eigenvalue1 = f1.readlines()
            eigenvalue1 = eval(pre_eigenvalue1[0])
            f1.close()
            
            
                
                
                
            if eigenvalue0[0]==-1:
                for ii in range(number0):
                    for iii in range(N):
                        feature_matrix[i-start][count] = 0
                        count = count + 1
                 
            else:
                #number0 = 2
                for ii in range(1,number0+1,add_value):
                    value = []
                    all_value = []
                    c0 = 0
                    for iii in range(len(eigenvalue0[ii])):
                        v = eigenvalue0[ii][iii]
                        if v<=0.000000001:
                            c0 = c0 + 1
                            all_value.append(0)
                        else:
                            value.append(v)
                            all_value.append(v)
                        #print(value)
                    feature_matrix[i-start][count] = c0
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_max(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_min(value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_mean(all_value)
                    count = count + 1
                        
                        
                    feature_matrix[i-start][count] = get_std(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_sum(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_generalized_mean_graph_energy(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_spectral_moment(all_value,-1)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_spectral_moment(all_value,2)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_spectral_moment(all_value,3)
                    count = count + 1
                    ###################################################################
             
            if eigenvalue1[0]==-1:
                for ii in range(number1):
                    for iii in range(N):
                        feature_matrix[i-start][count] = 0
                        count = count + 1
                 
            else:
                #number0 = 2
                for ii in range(1,number1+1,add_value):
                    value = []
                    all_value = []
                    c0 = 0
                    for iii in range(len(eigenvalue1[ii])):
                        v = eigenvalue1[ii][iii]
                        if v<=0.000000001:
                            c0 = c0 + 1
                            all_value.append(0)
                        else:
                            value.append(v)
                            all_value.append(v)
                        #print(value)
                    feature_matrix[i-start][count] = c0
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_max(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_min(value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_mean(all_value)
                    count = count + 1
                        
                        
                    feature_matrix[i-start][count] = get_std(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_sum(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_generalized_mean_graph_energy(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_spectral_moment(all_value,-1)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_spectral_moment(all_value,2)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_spectral_moment(all_value,3)
                    count = count + 1
            
                        
    filename = pre + 'pocket_feature/ligand_eigenvalue01_train.csv'
    np.savetxt(filename,feature_matrix,delimiter=',')
    


def test_feature_to_file(start,end,filtration0,filtration1,grid):
    row = end - start
    N = 10
    number0 = int (filtration0/grid )
    number1 = int (filtration1/grid)
    add_value = int(10 * grid)
    column = 7 * number0 * N + 7*number1*N
    
    feature_matrix = np.zeros((row,column))
    #pre0 = 'E:\\singapore\\pdb_bind\\' + Year + '\\' + 'max_cocycle_bar\\eigenvalue_' + str(cutoff) + '_' + str(filtration0) + '_' +  'zero\\' 
    #pre = 'E:\\singapore\\pdb_bind\\' + Year + '\\' + 'max_cocycle_bar\\eigenvalue_' + str(cutoff) + '_' + str(filtration0) + '_' + str(filtration1) + '\\' 
    pre_zero = pre + 'ligand_eigenvalue_' + str(filtration0) + '_zero/'
    pre_one = pre + 'ligand_eigenvalue_' + str(filtration1) + '_one/'
    
    
    for i in range(start,end):
        print(i)
        name = test_data[i]
        count = 0
        for typ in ['C','N','O','CN','CO','CNO','all']:
            filename0 = pre_zero + name + '_' + typ + '_' + 'eigenvalue_0D.txt'
            f0 = open(filename0)
            pre_eigenvalue0 = f0.readlines()
            eigenvalue0 = eval(pre_eigenvalue0[0])
            f0.close()
            
            filename1 = pre_one + name + typ + '_' + 'eigenvalue_1D.txt'
            f1 = open(filename1)
            pre_eigenvalue1 = f1.readlines()
            eigenvalue1 = eval(pre_eigenvalue1[0])
            f1.close()
            
            
                
                
                
            if eigenvalue0[0]==-1:
                for ii in range(number0):
                    for iii in range(N):
                        feature_matrix[i-start][count] = 0
                        count = count + 1
                 
            else:
                #number0 = 2
                for ii in range(1,number0+1,add_value):
                    value = []
                    all_value = []
                    c0 = 0
                    for iii in range(len(eigenvalue0[ii])):
                        v = eigenvalue0[ii][iii]
                        if v<=0.000000001:
                            c0 = c0 + 1
                            all_value.append(0)
                        else:
                            value.append(v)
                            all_value.append(v)
                        #print(value)
                    feature_matrix[i-start][count] = c0
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_max(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_min(value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_mean(all_value)
                    count = count + 1
                        
                        
                    feature_matrix[i-start][count] = get_std(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_sum(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_generalized_mean_graph_energy(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_spectral_moment(all_value,-1)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_spectral_moment(all_value,2)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_spectral_moment(all_value,3)
                    count = count + 1
                    ###################################################################
             
            if eigenvalue1[0]==-1:
                for ii in range(number1):
                    for iii in range(N):
                        feature_matrix[i-start][count] = 0
                        count = count + 1
                 
            else:
                #number0 = 2
                for ii in range(1,number1+1,add_value):
                    value = []
                    all_value = []
                    c0 = 0
                    for iii in range(len(eigenvalue1[ii])):
                        v = eigenvalue1[ii][iii]
                        if v<=0.000000001:
                            c0 = c0 + 1
                            all_value.append(0)
                        else:
                            value.append(v)
                            all_value.append(v)
                        #print(value)
                    feature_matrix[i-start][count] = c0
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_max(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_min(value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_mean(all_value)
                    count = count + 1
                        
                        
                    feature_matrix[i-start][count] = get_std(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_sum(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_generalized_mean_graph_energy(all_value)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_spectral_moment(all_value,-1)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_spectral_moment(all_value,2)
                    count = count + 1
                        
                    feature_matrix[i-start][count] = get_spectral_moment(all_value,3)
                    count = count + 1
            
                        
    filename = pre + 'pocket_feature/ligand_eigenvalue01_test.csv'
    np.savetxt(filename,feature_matrix,delimiter=',')
    




def get_combined_feature_dis_and_charge(typ):
    pre = '../data/' + Year + '/pocket_feature/'
    filename1 = pre + 'dis_' + typ + '.csv'
    filename2 = pre + 'charge_' + typ + '.csv'
    m1 = np.loadtxt(filename1,delimiter=',')
    m2 = np.loadtxt(filename2,delimiter=',')
    t1 = m1.shape
    t2 = m2.shape
    m = np.zeros((t1[0],t1[1]+t2[1] ))
    
    for i in range(t1[0]):
        
        m[i][ 0 : t1[1] ] = m1[i,: ]
        m[i][ t1[1] : t1[1]+t2[1] ] = m2[i,:]
            
    filename3 = pre + 'dis_charge_' + typ + '.csv'
    np.savetxt(filename3,m,delimiter=',')
    
def get_combined_feature_dis_and_charge_and_ligand(typ):
    pre = '../data/' + Year + '/pocket_feature/'
    filename1 = pre + 'dis_charge_' + typ + '.csv'
    filename2 = pre + 'ligand_eigenvalue01_' + typ + '.csv'
    m1 = np.loadtxt(filename1,delimiter=',')
    m2 = np.loadtxt(filename2,delimiter=',')
    t1 = m1.shape
    t2 = m2.shape
    m = np.zeros((t1[0],t1[1]+t2[1] ))
    
    for i in range(t1[0]):
        
        m[i][ 0 : t1[1] ] = m1[i,: ]
        m[i][ t1[1] : t1[1]+t2[1] ] = m2[i,:]
            
    filename3 = pre + 'dis_charge_ligand_' + typ + '.csv'
    np.savetxt(filename3,m,delimiter=',')






############################################################################################################
# machine_learning algorithm starts.
    
def gradient_boosting(X_train,Y_train,X_test,Y_test):
    params={'n_estimators': 40000, 'max_depth': 6, 'min_samples_split': 2,
                'learning_rate': 0.001, 'loss': 'ls','max_features':'sqrt','subsample':0.7}
    regr = GradientBoostingRegressor(**params)
    regr.fit(X_train,Y_train)
    pearson_coorelation = sp.stats.pearsonr(Y_test,regr.predict(X_test))
    mse1 = mean_squared_error(Y_test, regr.predict(X_test))
    mse2 = pow(mse1,0.5)
    #mse3 = mse2/0.7335
    mse3 = mse2
    return [pearson_coorelation[0],mse3]


def get_pearson_correlation(typ):
    pre = '../data/' + Year
    feature_matrix_of_train = np.loadtxt( pre + '/pocket_feature/' + typ + '_train.csv',delimiter=',' )
    target_matrix_of_train = np.loadtxt( pre + '/pocket_feature/' + 'target_matrix_of_train.csv',delimiter=',' )
    feature_matrix_of_test = np.loadtxt( pre + '/pocket_feature/' + typ + '_test.csv',delimiter=',' )
    target_matrix_of_test = np.loadtxt( pre + '/pocket_feature/' +  'target_matrix_of_test.csv',delimiter=',' )
    
    number = 10
    P = np.zeros((number,1))
    M = np.zeros((number,1)) 
    print(feature_matrix_of_train.shape,feature_matrix_of_test.shape)
    for i in range(number):
        [P[i][0],M[i][0]] = gradient_boosting(feature_matrix_of_train,target_matrix_of_train,feature_matrix_of_test,target_matrix_of_test)
        print(P[i])
    median_p = np.median(P)
    median_m = np.median(M)
    print('for data ' + Year + ', 10 results for ' + typ + '-model are:')
    print(P)
    print('median pearson correlation values are')
    print(median_p)
    print('median mean squared error values are')
    print(median_m)
    
    
############################################################################################################
# machine_learning algorithm ends.
    
    
extract_coordinate_of_ligand(0,2959)
eigenvalue_to_file(0,2959,10,5,0.1)
train_feature_to_file(0,2764,10,5,0.1)
test_feature_to_file(0,195,10,5,0.1)

get_combined_feature_dis_and_charge('train')
get_combined_feature_dis_and_charge('test')
get_pearson_correlation('dis_charge','')

get_combined_feature_dis_and_charge_and_ligand('train')
get_combined_feature_dis_and_charge_and_ligand('test')
get_pearson_correlation('dis_charge_ligand','')

