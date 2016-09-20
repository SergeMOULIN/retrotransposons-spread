# coding: utf-8 
import numpy as np
import time
from random import *
from munkres import Munkres, print_matrix
munk = Munkres()
from math import exp
from functions import Q, phi, ponderation, dico_densite, DistriDistance, TreeDistance, DistriBuild, TreeBuild, moyenne_saut_reelle
from functions import Ballayage_X, reduction_X
from case_study import Case1, Case2, Case3


def optim_p1(Grid,Case,n1,N2,N3,N4,nb_feuilles):
   '''
   First part of the optimization method. Optimization of the deterioration states distribution. 
   '''
    Delta = np.array(Grid[1], float) - np.array(Grid[0], float) # Size of the interval for each parameter.
    while sum (np.array(Grid[3], float)[1:4] < Delta[1:4]) > 0: # Stop condition based on the size of the intervals. 
        M = np.ones(N3)*float('inf')                            # The N3 smaller distances found yet. 
        N3_bis = (Grid[2][1]+1)*(Grid[2][2]+1)*(Grid[2][3]+1)   
        N3 = min(N3,N3_bis)
        top = np.ones((N3,6))                                   # The N3 best parameter sets.
        Best = np.ones(6)                                       # The best parameter set.
        for i1 in range (Grid[2][1]+1):
            for i2 in range (Grid[2][2]+1):
                for i3 in range (Grid[2][3]+1):
                    Dt = 0                                      # The sum of the obtained distances between simulations and the observed genome.
                    V = np.zeros(6)                        
                    P = [0, i1, i2, i3]
                    for i in range (1,4):
                        if Grid[2][i] == 0:
                            V[i] = Grid[0][i]
                        else:
                            V[i] = (float(P[i]*Delta[i])/Grid[2][i])+Grid[0][i]
                    for j in range (N2):    
                        (s,time)=DistriBuild(V[1],V[2],V[3],nb_feuilles)
                        D = DistriDistance(Case[0],s)     
                        Dt = Dt + D
                    i = np.argmax(M)
                    if Dt < M[i]:
                        M[i] = Dt
                        top[i] = V
        
        Min = float('inf') 
        for i in range (N3):
            V = top[i]
            Dt = 0
            for j in range (N4):
                (s,time)=DistriBuild(V[1],V[2],V[3],nb_feuilles)
                D = DistriDistance(Case[0],s)     
                Dt = Dt + D
            if Dt < Min:
                Best = V
                Min = Dt

        for i in range (1,4):
            Grid[0][i] = max(Best[i] - float(Delta[i])/(2*n1) , float(Grid[0][i])/2)
            Grid[1][i] = Best[i] + float(Delta[i])/(2*n1)
   
        Delta = np.array(Grid[1], float) - np.array(Grid[0], float)
   
    return (Best,Grid)


def optim_p2(Grid,Case,n1,N2,N3,N4,N5,nb_feuilles,v1,v2,v3,dico_genes,prop_cum_taille): 
    '''    
    The second part of the optimization.
    '''
    chromosome = '3L'      # For now it only applies to the 3L chromosome.
    MinX = 0               # The index of the smallest position for X.
    MaxX = len(Case[1])-1    # The index of the largest position for X.
    pond = ponderation(Case)
    dico_dens = dico_densite(dico_genes,Case) 
    Delta = np.array(Grid[1], float) - np.array(Grid[0], float)
    
    while (Grid[3][0] < Delta[0] or Grid[3][4] < Delta[4] or Grid[3][5] < Delta[5]):  
        M = np.ones(N3)*float('inf')
        listeX = Ballayage_X(Grid[0][0], Grid[1][0], Grid[2][0])
        N3_bis = (len(listeX))*(Grid[2][4]+1)*(Grid[2][5]+1)   
        N3 = min(N3,N3_bis)
        top = np.ones((N3,6))
        Best = np.ones(6)
        
        for i0 in listeX:
            for i4 in range (Grid[2][4]+1):
                for i5 in range(Grid[2][5]+1):
                    Dt = 0
                    V = [i0, v1, v2, v3, 0, 0]
                    P = [0, 0, 0, 0, i4, i5]
                    for i in [4,5]:
                        if Grid[2][i] == 0:
                            V[i] = Grid[0][i]
                        else:
                            V[i] = (float(P[i]*Delta[i])/Grid[2][i])+Grid[0][i]
                    for j in range (N2):
                        (s,t,time)=TreeBuild(chromosome,V[0],V[1],V[2],V[3],V[4],V[5],nb_feuilles,dico_genes,prop_cum_taille,Case,dico_dens)
                        D = TreeDistance(Case[0],Case[1],s,t.take([0],axis=1),nb_feuilles,pond)     
                        Dt = Dt + D
                    i = np.argmax(M)
                    if Dt < M[i]:
                        M[i] = Dt
                        top[i] = V

        Min = float('inf')
        for i in range (N3):
            V = top[i]
            Dt = 0
            for j in range (N4):
                (s,t,time)=TreeBuild(chromosome,V[0],V[1],V[2],V[3],V[4],V[5],nb_feuilles,dico_genes,prop_cum_taille,Case,dico_dens)
                D = TreeDistance(Case[0],Case[1],s,t.take([0],axis=1),nb_feuilles,pond)     
                Dt = Dt + D
            if Dt < Min:
                Best = V
                Min = Dt

        Grid[0][0],Grid[1][0] = reduction_X(Grid[0][0],Grid[1][0],Best[0],n1,Grid[2][0],len(Case[1]))
        for i in [4,5]:
            Grid[0][i] = max(Best[i] - float(Delta[i])/(2*n1) , Grid[0][i]/2)
            Grid[1][i] = Best[i] + float(Delta[i])/(2*n1)

        Delta = np.array(Grid[1], float) - np.array(Grid[0], float)

    Score = 0
    moyenne_saut = 0
    T_obs_moyen = 0
    V = Best
    for i in range (N5):
        (s,t,time)=TreeBuild(chromosome,V[0],V[1],V[2],V[3],V[4],V[5],nb_feuilles,dico_genes,prop_cum_taille,Case,dico_dens)
        D = TreeDistance(Case[0],Case[1],s,t.take([0],axis=1),nb_feuilles,pond)
        moyenne_saut = moyenne_saut + moyenne_saut_reelle(t)
        T_obs_moyen = T_obs_moyen + time
        Score = Score + D
    if N5 > 0:
        moyenne_saut = float(moyenne_saut)/N5
        T_obs_moyen = float(T_obs_moyen)/N5
        
    return (Best, Grid, Score, moyenne_saut, T_obs_moyen)

def optim(Grid,Case,n1,N2,N3,N4,N5,nb_feuilles,dico_genes,prop_cum_taille):
    (Best, Grid_2) = optim_p1(Grid,Case,n1,N2,N3,N4,nb_feuilles)
    (Best,Grid_3, Score, moyenne_saut, T_obs_moyen) = optim_p2(Grid_2,Case,n1,N2,N3,N4,N5,nb_feuilles,Best[1],Best[2],Best[3],dico_genes,prop_cum_taille)
    return (Best,Grid_3, Score, moyenne_saut, T_obs_moyen)
