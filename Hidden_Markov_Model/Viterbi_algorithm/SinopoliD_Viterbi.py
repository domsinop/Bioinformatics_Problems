import numpy as np
import pandas as pd
import re
from copy import copy
import copy
import glob
import sys
import os

def main():
    split_path = os.path.split(sys.argv[1])
    file_name = sys.argv[1]
    data = np.array(pd.read_csv(file_name, sep = '/n', engine = 'python',header = None))
    #Initial Conditions
    A = np.array([[(2/5),(2/5),(1/5),0],[0,(3/5),(2/5),0],[0,0,(2/5),(3/5)],[(1/2),0,0,(1/2)]])
    B = np.array([[(1/3),(1/6),0,(1/2)],[(1/4),(1/4),(1/4),(1/4)],[0,(3/16),(1/16),(3/4)],[(4/5),(1/10),(1/10),(0)]])
    Pi = np.array([1,0,0,0])
    create_txt_file(get_dna(data),viterbi(get_dna(data),A,B,Pi),split_path)
    
def get_dna(data):
    dna_positions = []
    dna = []
    for i in range(len(data)):
        if(re.search('>',data[i][0])):
            dna_positions.append(i)
    for i in range(len(dna_positions)):
        for j in range(dna_positions[i]+1,dna_positions[i]+2):
            dna.append(np.reshape(np.array(list(data[j][0])),(np.array(list(data[j][0])).size,1)))
    return dna
def viterbi(dna,A,B,Pi):
    scores = []
    paths = []
    for i in range(len(dna)):
        path = []
        alpha = np.zeros((dna[i].shape[0],A.shape[0]))
        keys = np.zeros((dna[i].shape[0],A.shape[0]))
        dna_num = dna_to_num(dna[i])
        #initlize first row
        for j in range(len(Pi)):
            alpha[0,j] = Pi[j]*B[j,int(dna_num[j])]
        #fill in rest
        for k in range(1,alpha.shape[0]):
            for l in range(alpha.shape[1]):
                value = 0
                key = 0
                for z in range(alpha.shape[1]):
                    test_value =  (alpha[k-1,z]*A[z,l]*B[l,int(dna_num[k])])
                    if(test_value > value):
                        value = test_value
                        key = z
                        
                alpha[k,l] = value
                if(value != 0):
                    keys[k,l] = key+1
                else:
                    keys[k,l] = key
             
        #get path
        path = np.zeros(len(dna[i]))
        path[0] = str(np.argmax(alpha[alpha.shape[0]-1,:])+1)
        w = (len(dna[i])) - 1
        for b in range(len(dna[i])-1):
            path[b+1] = str(keys[w,int(path[b]-1)])
            w = w - 1
        path = np.flip(np.array(path))
        paths.append(path)
      
    return paths

def dna_to_num(dna_seq):
    dna_num = np.zeros(len(dna_seq))
    for i in range(len(dna_seq)):
        if(dna_seq[i] == 'A'):
            dna_num[i] = 0
        if(dna_seq[i] == 'C'):
            dna_num[i] = 1
        if(dna_seq[i] == 'G'):
            dna_num[i] = 2
        if(dna_seq[i] == 'T'):
            dna_num[i] = 3
    return dna_num

def create_txt_file(dna,paths,split_path):
    #file = open("Viterbi_Path.txt", "w")
    file = open(split_path[0]+'/'+str.split(split_path[1],'.')[0]+'_Viterbi_'+"results.txt", "w")
    for i in range(len(dna)):
        np.savetxt(file,np.append(np.reshape(['Seq',str(i+1)+':     '],[1,2]),np.transpose(np.array(dna[i])),1),fmt = '%s')
        new = get_list(paths[i])
        np.savetxt(file,np.reshape(['State:     ',' '.join(map(str,new))],[1,2]),fmt ='%s')
        np.savetxt(file,[' '],fmt = '%s')
        np.savetxt(file,[' '],fmt = '%s')
    file.close()


def get_list(dna):
    new = []
    for i in range(len(dna)):
        new.append(str(int(dna[i])))
    return new

if __name__ == "__main__":
    main()