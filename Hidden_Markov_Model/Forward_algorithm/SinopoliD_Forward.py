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
    create_txt_file(get_dna(data),forward(get_dna(data),A,B,Pi),split_path)


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

def forward(dna,A,B,Pi):
    scores = []
    for i in range(len(dna)):
        alpha = np.zeros((dna[i].shape[0],A.shape[0]))
        dna_num = dna_to_num(dna[i])
        #initlize first row
        for j in range(len(Pi)):
            alpha[0,j] = Pi[j]*B[j,int(dna_num[j])]
        #fill in rest
        for k in range(1,alpha.shape[0]):
            for l in range(alpha.shape[1]):
                value = 0
                for z in range(alpha.shape[1]):
                    value = value + (alpha[k-1,z]*A[z,l]*B[l,int(dna_num[k])])
                alpha[k,l] = value
        #get score
        score = 0
        for b in range(alpha.shape[1]):
            score = score + alpha[alpha.shape[0]-1,b]
        scores.append(score)
      
    return scores

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

def create_txt_file(dna,scores,split_path):
    #file = open("forward_scores.txt", "w")
    file = open(split_path[0]+'/'+str.split(split_path[1],'.')[0]+'_Forward_'+"results.txt", "w")
    for i in range(len(dna)):
        np.savetxt(file,np.append(np.reshape(['Seq',str(i+1)+':'],[1,2]),np.transpose(np.array(dna[i])),1),fmt = '%s')
        np.savetxt(file,np.reshape(['PR(Seq'+str(i+1)+'|Model):',scores[i]],[1,2]),fmt ='%s')
        np.savetxt(file,[' '],fmt = '%s')
        np.savetxt(file,[' '],fmt = '%s')
    file.close()

if __name__ == "__main__":
    main()
