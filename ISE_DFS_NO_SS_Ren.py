#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import numpy.linalg as npla
import scipy.stats as st
import math
# import sympy as smp
import random as rd
import matplotlib.pyplot as plt


# In[2]:


class def_event:
    def __init__(self,E,t):
        self.evnt = E
        self.time = t

class def_serial_episode:
    def __init__(self, evnts_ordr):
        self.evnt = evnts_ordr[:]
#         self.length = len(evnts_ordr)
        self.freq = 0

class def_window:
    def __init__(self, st_time, end_time):
        self.ts = st_time
        self.te = end_time


# In[3]:


class def_episode:
    def __init__(self, evnts, edge_set):#, lastevnt, pred, succ):
        self.freq = 0
        self.evnts = evnts.copy()
        self.edges = edge_set.copy()
#         self.lastevnt = lastevnt
#         self.pred = pred.copy()
#         self.succ = succ.copy()

class def_NOMDFA:
    def __init__(self):
        self.s = (set(),'')
        self.S = list([self.s])
        self.F = list()
        self.D = {}


# In[19]:


def NO_Span_Recursive(ep_p,ep_p_MO_list,Tx,FT):
    global F1
    F = list()
    closed = 1
    for E in F1:
        if E not in ep_p.evnt:
            evnts_temp = ep_p.evnt.copy()
            evnts_temp.append(E)
            ep = def_serial_episode(evnts_temp)
            # F.append(ep) # To find candidates
            ep_MO_list = Compute_MO_Span_List(ep_p_MO_list,E)
            ep.freq = Find_NO_Count(ep_MO_list) #len(ep_MO_list)
            if ep.freq == ep_p.freq:
                closed = 0
            if ep.freq >= FT:
                F += NO_Span_Recursive(ep,ep_MO_list,Tx,FT)
#     if closed == 1:
    F.append(ep_p)
    return(F)


# In[5]:


def Compute_MO_Span_List(ep_p_MO_list,E):
    p_size = len(ep_p_MO_list)
    E_MO_list = MO_list[evnt_id[E]]
    E_size = len(E_MO_list)
    ep_MO_list = list()
    p1=0
    p2=0
    while (p1<p_size and p2<E_size):
        while p2 < E_size and E_MO_list[p2].ts <= ep_p_MO_list[p1].te :
            p2 += 1
        if p2 < E_size and E_MO_list[p2].te-ep_p_MO_list[p1].ts < Tx:
            while(p1 < p_size and ep_p_MO_list[p1].te < E_MO_list[p2].ts):
                p1 += 1
            ep_MO_list.append(def_window(ep_p_MO_list[p1-1].ts,E_MO_list[p2].te))
        else:
            p1 += 1
    return(ep_MO_list)


# In[6]:


def Find_NO_Count(ep_MO_list):
    l = len(ep_MO_list)
    if l == 0:
        return l
    ts = ep_MO_list[0].ts
    te = ep_MO_list[0].te
    fno = 1
    p = 1
    while (p<l):
        if ep_MO_list[p].ts >= ts and ep_MO_list[p].te <= te:
            ts = ep_MO_list[p].ts
            te = ep_MO_list[p].te
        elif ep_MO_list[p].ts > te:
            fno += 1
            ts = ep_MO_list[p].ts
            te = ep_MO_list[p].te
        p += 1
    return fno


# In[7]:


def Find_Closed(Freq_ep):
    hasher = list()
    hashtable = list()
    for ep in Freq_ep:
        fn = ep.freq
        flg = 1
        for k in range(len(hasher)):
            if fn == hasher[k]:
                flg = 0
                hashtable[k].append(ep)
                break
        if flg == 1:
            hasher.append(fn)
            hashtable.append(list([ep]))

    Freq_closed_ep = list()
    for k in range(len(hasher)):
        len_hash_k = len(hashtable[k])
        p = 0
        while p < len_hash_k:
            q = 0
            if q == p:
                q += 1
            while q < len_hash_k:
                len_p = len(hashtable[k][p].evnt)
                len_q = len(hashtable[k][q].evnt)
                if len_p < len_q:
                    sub_ep = 0
                    pi = 0
                    qi = 0
                    while(1):
                        if hashtable[k][p].evnt[pi] == hashtable[k][q].evnt[qi]:
                            pi += 1
                            if pi == len_p:
                                sub_ep = 1
                                break
                        qi += 1
                        if qi == len_q:
                            break
                    if sub_ep == 1:
                        hashtable[k].remove(hashtable[k][p])
                        len_hash_k = len(hashtable[k])
                        q -= 1
                        break
                q += 1
                if q == p:
                    q += 1
            if q >= len_hash_k:
                p += 1

    for k in range(len(hasher)):
        for l in range(len(hashtable[k])):
            Freq_closed_ep.append(hashtable[k][l])

    return Freq_closed_ep


# In[8]:


def Construct_NOMDFA(alpha,prob):
    global alphabet
    alpha_size = len(alpha.evnts)
    pi = {}
    for e in alpha.evnts:
        pi[e] = set()
    for (e1,e2) in alpha.edges:
        pi[e2].add(e1)
    DFA = def_NOMDFA()
    n_states = 1
    for Q in DFA.S:
        DFA.D[str(Q)] = {}
        for ev in alphabet:
            if len(Q[0]) != alpha_size:
                Qs_new = Q[0].copy()
                Ql_new = ev
                if (ev in alpha.evnts):
                    if not (alpha.evnts-Q[0]).intersection(pi[ev]):
                        Qs_new.add(ev)
                if Q[1] == '' or prob[Q[1]][Ql_new] != 0:
                    Q_new = (Qs_new,Ql_new)
                    DFA.D[str(Q)][ev] = Q_new
                    if Q_new not in DFA.S:
                        DFA.S.append(Q_new)
            else:
                DFA.D[str(Q)][ev] = DFA.D[str(DFA.s)][ev]
            
#             print(DFA.D[str(Q)])
        if len(Q[0]) == alpha_size:
            DFA.F.append(Q)
    return DFA


# In[9]:


def GetTransitionMatrix(DFA,prob):
    global final_states
    T0 = list()
    for state in DFA.S:
        T0.append([0]*len(DFA.S))
        for ev in DFA.D[str(state)].keys():
            if state[1] == '':
                T0[-1][DFA.S.index(DFA.D[str(state)][ev])] = 1/len(list(DFA.D[str(state)].keys())) #prob[DFA.S[final_states[0]][1]][ev]
            else:
                T0[-1][DFA.S.index(DFA.D[str(state)][ev])] = prob[state[1]][ev]
    return T0


# In[10]:


def GetStats_Renewal(T0,final_states):

    nstates = len(T0)
    nfinals = len(final_states)
    
    if nfinals == 0:
        return(0,0)
    if nfinals != 1:
        return(-1,-1)
    f = final_states[0]
    
    T1 = np.matrix(T0)
#     print(T1,final_states)
    T1[:,final_states[0]] *= 0
    L = np.eye(nstates)
    
    m1 = npla.solve(L-T1,np.ones(nstates))
    m2 = npla.solve(L-T1,2*(m1-np.ones(nstates)))
#     m2 = npla.solve(L-T1,2*np.matmul(T1,m1).T)

    MX1 = m1[f]**(-1)
    VX1 = m1[f]**(-3)*(m2[f]-m1[f]**2+m1[f])

    return MX1, VX1


# In[20]:


evnt_strm = list()
print('Enter the link of the txt file containing the event stream (for example: ./Data/Datastream_M_40_n_10000_emb_2_N_4_p_0.02_eta_0.5_probs_dense.txt)')
name = input()
with open(name,'r') as f:
    for line in f:
        entry = line.split(',')
        evnt_strm.append(def_event(entry[0],int(entry[1][:-1])))

M = len(evnt_strm)
alph = list()
for m in range(len(evnt_strm)):
    A = evnt_strm[m].evnt
    flg = 1
    for a in alph:
        if a == A:
            flg = 0
            break
    if flg == 1:
        alph.append(A)
LA = len(alph) 

A = alph[0]
evnt_id = {A: 0}
count = 1
for A in alph[1:]:
    evnt_id[A] = count
    count += 1

MO_list = list()
for i in range(LA):
    MO_list.append(list())

Freq_ep = list()
F1 = alph
F1_ep = list()
for A in alph:
    F1_ep.append(def_serial_episode(list([A])))


# FT = 18
# Tx = 6
print('Enter the Expiry time "Tx":')
Tx = float(input())
print('Enter the frequency threshold "FT":')
FT = float(input())

for event in evnt_strm:
    E = event.evnt
    t = event.time
    F1_ep[evnt_id[E]].freq += 1
    MO_list[evnt_id[E]].append(def_window(t,t))

for E in F1:
    if F1_ep[evnt_id[E]].freq >= FT:
        Freq_ep += NO_Span_Recursive(F1_ep[evnt_id[E]],MO_list[evnt_id[E]],Tx,FT)
        
print('Total Candidates: ', len(Freq_ep) + len(F1))


    
Freq_closed_ep = Find_Closed(Freq_ep)

print('Total closed frequent episodes: ', len(Freq_closed_ep))
print('The following are the closed frequent episodes discovered')
for ep in Freq_closed_ep:
    print(ep.evnt ,':', ep.freq)

# print(len(Freq_ep),len(Freq_closed_ep))


# In[225]:


serial_episodes = list()
episodes = list()
for ep in Freq_closed_ep:
    if len(ep.evnt) >=2 and len(ep.evnt) <=4:
        serial_episodes.append(ep)
        
        temp_evnts = set(ep.evnt)
        temp_edges = set()
        for k in range(len(ep.evnt)):
            for l in range(k+1,len(ep.evnt)):
                temp_edges.add((ep.evnt[k],ep.evnt[l]))
        episodes.append(def_episode(temp_evnts,temp_edges))
#         print(ep.evnt ,':', ep.freq)

M = len(alph)
Probs = np.zeros((M,M))

for k in range(1,len(evnt_strm)):
    E_now = evnt_strm[k].evnt
    E_prev = evnt_strm[k-1].evnt
    Probs[evnt_id[E_prev],evnt_id[E_now]] += 1

for k in range(M):
    temp_sum = np.sum(Probs[k,:])
    for l in range(M):
        Probs[k,l] /= temp_sum

Prob_Ep_Events = {}
for k in range(M):
    Prob_Ep_Events[alph[k]] = {}
    for l in range(M):
        Prob_Ep_Events[alph[k]][alph[l]] = Probs[k,l]


# In[226]:


global alphabet
alphabet = alph

mean_estimate = list()
var_estimate = list()

for episode_id in range(len(episodes)):
    DFA = Construct_NOMDFA(episodes[episode_id],Prob_Ep_Events)
    
    final_states = list()
    for fstates in DFA.F:
        final_states.append(DFA.S.index(fstates))
    final_states.sort()
    
    T0 = GetTransitionMatrix(DFA,Prob_Ep_Events)
    MX_Ren, VX_Ren = GetStats_Renewal(T0, final_states)
    # print(MX_Ren, VX_Ren)
    
    mean_estimate.append(MX_Ren)
    var_estimate.append(VX_Ren)

for k in range(len(serial_episodes)):
    print(serial_episodes[k].evnt, ' : ', serial_episodes[k].freq, ' , ', mean_estimate[k])


# In[227]:


n = len(evnt_strm)
SSTh = list()
c = 0.01
for k in range(len(serial_episodes)):
    value = mean_estimate[k]*n + math.sqrt(var_estimate[k]*n)*st.norm.ppf(1-c)
    SSTh.append(value)

count = 0
for k in range(len(serial_episodes)):
    print(serial_episodes[k].evnt, ' : ', serial_episodes[k].freq, ' , ', SSTh[k], mean_estimate[k])
    if serial_episodes[k].freq > SSTh[k]:
        count += 1


# In[228]:


print(count, len(serial_episodes))
for k in range(len(serial_episodes)):
    if serial_episodes[k].freq > SSTh[k]:
        print(serial_episodes[k].evnt, ' : ', serial_episodes[k].freq, ' , ', SSTh[k])


# In[229]:


print(len(serial_episodes),count)
count_FE_2 = 0
count_FE_3 = 0
count_FE_4 = 0
count_SS_2 = 0
count_SS_3 = 0
count_SS_4 = 0
for k in range(len(serial_episodes)):
    if len(serial_episodes[k].evnt) == 2:
        count_FE_2 += 1
#         if serial_episodes[k].freq > SSTh[k]:
#             count_SS_2 += 1
    if len(serial_episodes[k].evnt) == 3:
        count_FE_3 += 1
#         if serial_episodes[k].freq > SSTh[k]:
#             count_SS_3 += 1
    if len(serial_episodes[k].evnt) == 4:
        count_FE_4 += 1
#         if serial_episodes[k].freq > SSTh[k]:
#             count_SS_4 += 1
print("c: ", c)
print('number of Frequent episodes and number of episodes that qualified under Statistical significance threshold')
print(count_FE_2, count_SS_2)
print(count_FE_3, count_SS_3)
print(count_FE_4, count_SS_4)






