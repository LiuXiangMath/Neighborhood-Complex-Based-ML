# -*- coding: utf-8 -*-

import numpy as np
import math
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error
import scipy as sp

Year = '2013'


def get_combined_feature_dis_and_charge(typ):
    pre = '../data/' + Year + '/pocket_feature/'
    filename1 = pre + 'laplacian_' + typ + '_10_10_0.1.csv'
    filename2 = pre + 'charge_laplacian_' + typ + '_50_10_10_0.02.csv'
    m1 = np.loadtxt(filename1,delimiter=',')
    m2 = np.loadtxt(filename2,delimiter=',')
    t1 = m1.shape
    t2 = m2.shape
    m = np.zeros((t1[0],t1[1]+t2[1] ))
    
    for i in range(t1[0]):
        
        m[i][ 0 : t1[1] ] = m1[i,: ]
        m[i][ t1[1] : t1[1]+t2[1] ] = m2[i,:]
            
    filename3 = pre + 'mix_dis_and_charge_laplacian_' + typ + '_50_10_10_0.02.csv'
    np.savetxt(filename3,m,delimiter=',')
    
def get_combined_feature_dis_and_charge_and_ligand(typ):
    pre = '../data/' + Year + '/pocket_feature/'
    filename1 = pre + 'mix_dis_and_charge_laplacian_' + typ + '_50_10_10_0.02.csv'
    filename2 = pre + '10_ligand_' + typ + '.csv'
    m1 = np.loadtxt(filename1,delimiter=',')
    m2 = np.loadtxt(filename2,delimiter=',')
    t1 = m1.shape
    t2 = m2.shape
    m = np.zeros((t1[0],t1[1]+t2[1] ))
    
    for i in range(t1[0]):
        
        m[i][ 0 : t1[1] ] = m1[i,: ]
        m[i][ t1[1] : t1[1]+t2[1] ] = m2[i,:]
            
    filename3 = pre + 'mix_dis_and_charge_and_ligand_laplacian_' + typ + '_50_10_10_0.02.csv'
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


def get_pearson_correlation(typ,pref):
    pre = '../data/' + Year
    feature_matrix_of_train = np.loadtxt( pre + '/pocket_feature/' + pref +'mix_' + typ + '_laplacian_train_50_10_10_0.02.csv',delimiter=',' )
    target_matrix_of_train = np.loadtxt( pre + '/pocket_feature/' + 'target_matrix_of_train.csv',delimiter=',' )
    feature_matrix_of_test = np.loadtxt( pre + '/pocket_feature/' + pref + 'mix_' + typ + '_laplacian_test_50_10_10_0.02.csv',delimiter=',' )
    target_matrix_of_test = np.loadtxt( pre + '/pocket_feature/' +  'target_matrix_of_test.csv',delimiter=',' )
    
    number = 10
    P = np.zeros((number,1))
    M = np.zeros((number,1)) 
    
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
    
get_combined_feature_distance_and_charge('train')
get_combined_feature_dis_and_charge('test')
get_pearson_correlation('dis_and_charge','')

get_combined_feature_distance_and_charge_and_ligand('train')
get_combined_feature_dis_and_charge_and_ligand('test')
get_pearson_correlation('dis_and_charge_and_ligand','')




