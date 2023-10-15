#!/usr/bin/env python
# coding: utf-8

# In[0]:


import numpy as np
import random as rd


# In[1]: For dense prabability matrix

M = 40
baseprob = 1/M
Probs = np.zeros((M,M))

for k in range(M):
    oncemore = True
    while(oncemore):
        for l in range(M-1):
            Probs[k,l] = round(rd.gauss(1/M,sigma=0.01),5)
            if Probs[k,l] < 0 or Probs[k,l] > 1:
                oncemore = True
                break
            else:
                oncemore = False
        if not oncemore:
            Probs[k,M-1] = (1-np.sum(Probs[k,:M-1]))
            if Probs[k,M-1] < 0 or Probs[k,l] > 1:
                oncemore = True
print(Probs)
np.save('Markov_M_'+str(M)+'_probs_dense',Probs)