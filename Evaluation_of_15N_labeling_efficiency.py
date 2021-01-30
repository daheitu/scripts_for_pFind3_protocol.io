# -*- coding: utf-8 -*-
"""
Created on Mon Jul 08 09:19:15 2019

@author: shaoguangcan
"""
import numpy as np
path = 'pQuant.spectra'
def Evaluation_of_15N_labeling_efficiency():
    w = open(r'Evaluation_of_15N_labeling_efficiency.txt','w')
    Title1 = ['Total number','Range [lowest, highest]','Middle50 [25%, 75%]','Mean value','Median value']
    Title2 = ['Number','efficiency']
    w.write('\t'.join(Title1) + '\n')
    n=0
    all_count = []
    All_title_count = []
    all_title_count = []
    with open(path,'r') as f:
        tmpline_cell=[]
        l = f.readlines()

        for i in range(0, len(l)):
            line = l[i].rstrip('\n')
            cell = line.split('\t')
            if 'I,Q,02' in cell[0]:
                n = n+1
                number = ('Value'+ str(n))
                all_title_count.append(number)
                all_title_count.append(cell[1])
                All_title_count.append(all_title_count)
                all_title_count =[]
                all_count.append(cell[1])

        all_count =list(map(float, all_count))
        count_average = np.mean(all_count)
        count_median = np.median(all_count)
        count_50 = np.percentile(all_count, (25, 75), interpolation='midpoint')
        w.write('\t'.join([str(n),str([min(all_count),max(all_count)]),str(count_50),str(count_average),str(count_median)]) +'\n'+'\n')
        w.write('\t'.join(['Number','Labeling efficiency'])+ '\n')
        for x in range(0,len(All_title_count)):
            w.write('\t'.join(All_title_count[x]) + '\n')
        print('Good Job!')

if __name__ == '__main__':
    Evaluation_of_15N_labeling_efficiency()
