import numpy as np
import pandas as pd
import re
import copy
import glob
import sys
import os



# In[48]:



def main():
  split_path = os.path.split(sys.argv[1])
  file_name = sys.argv[1]
  data = np.array(pd.read_csv(file_name, sep = '/n', engine = 'python',header = None))
  dna_1,dna_2 = get_dna(data)
  scores = get_score(dna_1,dna_2)
  sim,bc = get_sim(dna_1,dna_2,scores)
  s = {}
  k = 1
  key = list(bc.keys())[-1]
  s[1] = [[],[]]
  loop_d(bc,key,dna_1,dna_2,s,k)
  [m,n] = sim.shape
  create_txt_file(s,sim[m-1,n-1],split_path)


# In[ ]:





# In[49]:


def get_dna(data):
    dna_positions = []
    dna_1 = []
    dna_2 = []
    for i in range(len(data)):
        if(re.search('>',data[i][0])):
            dna_positions.append(i)
    for j in range(dna_positions[0]+1,dna_positions[1]):
        dna_1.append(list(data[j][0]))
    for k in range(dna_positions[1]+1,len(data)):
        dna_2.append(list(data[k][0]))
    dna_1 = np.reshape(np.array(dna_1),(np.array(dna_1).size,1))
    dna_2 = np.reshape(np.array(dna_2),(np.array(dna_2).size,1))
    return dna_1,dna_2


# In[3]:


def get_score(dna_1,dna_2):
    scores =  np.zeros((len(dna_2),len(dna_1)))
    for i in range(len(dna_2)):
        for j in range(len(dna_1)):
            if(dna_1[j] == dna_2[i]):
                scores[i,j] = 1
            else:
                scores[i,j] = -1
    return scores


# In[5]:


def get_sim(dna_1,dna_2,scores,d = -2):
    bread_crumbs = {}
    sim = np.zeros((len(dna_2)+1,len(dna_1)+1))
    m,n = sim.shape
    for i in range(n):
        sim[0,i] = d * i
    for j in range(m): 
        sim[j,0] = d * j
    for i in range(m-1):
        for j in range(n-1):
            locations = [];
            diag = sim[i,j]+scores[i,j]
            up = sim[i,j+1] + d
            left = sim[i+1,j] + d
            sim[i+1,j+1] = max(diag,left,up)
            if max(diag,left,up) == diag:
                locations.append([i,j,'d'])
            if max(diag,left,up) == left:
                locations.append([i+1,j,'l'])
            if max(diag,left,up) == up:
                locations.append([i,j+1,'u'])
            bread_crumbs[(i+1,j+1)] = locations
    return sim,bread_crumbs


# In[37]:


def loop_d(bc,key,dna_1,dna_2,s,k):
    if key[0] == 0 or key[1] == 0:
        if (key[0] == 0 and key[1] != 0) or (key[0] != 0 and key[1] == 0):
            get_last_letter(key,dna_1,dna_2,s,k)
            if(key[0] > 1 or key[1] > 1):
                if(key[0] > 1):
                    new_0 = key[0] -1
                    n_key = (new_0,key[1])
                    loop_d(bc,n_key,dna_1,dna_2,s,k)
                elif(key[1] > 1):
                    new_1 = key[1] -1
                    n_key = (key[0],new_1)
                    loop_d(bc,n_key,dna_1,dna_2,s,k)
                    
    else:
        if(len(bc[key]) == 1):
            get_letter(bc[key],dna_1,dna_2,s,k)
            key = (bc[(key)][0][0],bc[key][0][1])
            loop_d(bc,key,dna_1,dna_2,s,k)
            
        elif(len(bc[key])== 2):
            key_1 = (bc[(key)][0][0],bc[key][0][1])
            key_2 = (bc[(key)][1][0],bc[key][1][1])
            new_k = list(s.keys())[-1]+1
            s[new_k] = copy.deepcopy(s[k])
            get_letter([bc[key][0]],dna_1,dna_2,s,k)
            get_letter([bc[key][1]],dna_1,dna_2,s,new_k)
            loop_d(bc,key_1,dna_1,dna_2,s,k)
            loop_d(bc,key_2,dna_1,dna_2,s,new_k)
            
        elif(len(bc[key])== 3):
            key_1 = (bc[(key)][0][0],bc[key][0][1])
            key_2 = (bc[(key)][1][0],bc[key][1][1])
            key_3 = (bc[(key)][2][0],bc[key][2][1])
            new_k_1 = list(s.keys())[-1]+1
            new_k_2 = new_k_1 + 1
            s[new_k_1] = copy.deepcopy(s[k])
            s[new_k_2] = copy.deepcopy(s[k])
            get_letter([bc[key][0]],dna_1,dna_2,s,k)
            get_letter([bc[key][1]],dna_1,dna_2,s,new_k_1)
            get_letter([bc[key][2]],dna_1,dna_2,s,new_k_2)
            loop_d(bc,key_1,dna_1,dna_2,s,k)
            loop_d(bc,key_2,dna_1,dna_2,s,new_k_1)
            loop_d(bc,key_2,dna_1,dna_2,s,new_k_2)


# In[11]:


def get_letter(data,dna_1,dna_2,s,b):
    
    seq_1 = copy.deepcopy(s[b][0])
    seq_2 = copy.deepcopy(s[b][1])
        
    m = data[0][0]
    n = data[0][1]
    if data[0][2] == 'd':
        seq_1.append(dna_1[n])
        seq_2.append(dna_2[m])
    if data[0][2] == 'u':
        seq_1.append(['_'])
        seq_2.append(dna_2[m])
    if data[0][2] == 'l':
        seq_1.append(dna_1[n])
        seq_2.append(['_'])
     
    s[b] = [seq_1,seq_2]


# In[40]:


def get_last_letter(key,dna_1,dna_2,s,b):
    
    seq_1 = copy.deepcopy(s[b][0])
    seq_2 = copy.deepcopy(s[b][1])
    
    m = key[0]
    n = key[1]
    
    if (m == 0 and n != 0):
        seq_1.append(dna_1[n-1])
        seq_2.append(['_'])
    if (m != 0 and n == 0):
        seq_1.append(['_'])
        seq_2.append(dna_2[m-1])
    
    s[b] = [seq_1,seq_2]


# In[31]:


def compare_seq(seq_1,seq_2):
    compare = []
    score = 0;
    for i in range(len(seq_1)):
        if(seq_1[i] == seq_2[i]):
            compare.append('|')
            score = score + 1
        elif(seq_1[i] == '_' or seq_2[i] == '_'):
            compare.append('_')
            score = score - 2
        elif(seq_1[i] != seq_2[i]):
            compare.append(';')
            score = score - 1
    return np.array(compare),score


# In[47]:


def create_txt_file(s,score,split_path):
    file = open(split_path[0]+'/'+str.split(split_path[1],'.')[0]+'_'+"results.txt", "w")
    np.savetxt(file,np.reshape(['Score=',score],[1,2]),fmt = '%s')
    np.savetxt(file,[' '],fmt = '%s')
    for i in range(len(s)):
        x = np.flip(np.array(s[i+1][0]))
        y = np.flip(np.array(s[i+1][1]))
        s1 = ['>Seq_1','_',i+1,' ']
        s2 = ['>Seq_2','_',i+1,' ']
        c,score_s = compare_seq(x,y)
        if(i+1 < 10):
            spaces = ['            ']
        else:
            spaces = ['             ']
        np.savetxt(file,np.reshape(np.append(s1,x),[1,len(x)+len(s1)]),fmt ='%s')
        np.savetxt(file,np.reshape(np.append(spaces,c),[1,len(c)+1]),fmt='%s')
        np.savetxt(file,np.reshape(np.append(s2,y),[1,len(y)+len(s2)]),fmt ='%s')
        np.savetxt(file,np.reshape(['Score=',score_s],[1,2]),fmt = '%s')
        np.savetxt(file,[' '],fmt = '%s')
        np.savetxt(file,[' '],fmt = '%s')
    file.close()


# In[ ]:



if __name__ == "__main__":
    main()

