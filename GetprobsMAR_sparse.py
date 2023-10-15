#!/usr/bin/env python
# coding: utf-8

# In[0]:


import numpy as np
import random as rd


# In[2]:  For sparse probability matrix


M = 40
Probs = np.zeros((M,M))

n_clus = int(M/4)
clus_id = rd.choices(range(n_clus),[1]*n_clus,k=M)
clus = []
for c in range(n_clus):
    clus.append([])
for k in range(M):
    clus[clus_id[k]].append(k)

print(clus)


info = str(clus) + ' \n\n'
for cli in range(len(clus)):
    cl = clus[cli]
    if len(cl) == 0:
        continue
    if len(cl) == 1:
        while True:
            Probs[cl[0],cl[0]] = round(rd.gauss(1/(len(cl)+1),sigma=0.01),5)
            if Probs[cl[0],cl[0]] > 0  and Probs[cl[0],cl[0]] < 1:
                break
        cli_next = (cli + 1)%n_clus
        while(len(clus[cli_next]) == 0):
            cli_next = (cli_next + 1)%n_clus
        Probs[cl[0],rd.choices(clus[cli_next],[1]*len(clus[cli_next]),k=1)[0]] = 1 - Probs[cl[0],cl[0]]
        continue
    outlinker = rd.choices(cl,[1]*len(cl),k=1)[0]
    print(cl, outlinker)
    info += str(cl) + ' ' + str(outlinker) + ' \n'
    for k in cl:
        if k == outlinker:
            n_probs = len(cl)+1
        else:
            n_probs = len(cl)
        oncemore = True
        while(oncemore):
            temp = []
            for ki in range(n_probs-1):
                temp_prob = round(rd.gauss(1/n_probs,sigma=0.01),5)
                if temp_prob < 0 or temp_prob > 1:
                    oncemore = True
                    break
                else:
                    temp.append(temp_prob)
                    oncemore = False
            if not oncemore:
                temp_prob = (1-sum(temp))
                if temp_prob < 0 or temp_prob > 1:
                    oncemore = True
                else:
                    temp.append(temp_prob)
        for ki in range(len(cl)):
            Probs[k,cl[ki]] = temp[ki]
        if k == outlinker:
            cli_next = (cli + 1)%n_clus
            while(len(clus[cli_next]) == 0):
                cli_next = (cli_next + 1)%n_clus
            Probs[k,rd.choices(clus[cli_next],[1]*len(clus[cli_next]),k=1)[0]] = temp[-1]
        

print(Probs) 
info += '\n' + str(Probs)

np.save('Markov_M_'+str(M)+'_probs_sparse',Probs)
with open('Markov_M_'+str(M)+'_probs_sparse_info.txt','w') as f:
    f.write(info)