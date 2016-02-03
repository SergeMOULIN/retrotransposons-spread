# coding: utf-8 
import numpy as np
from random import *
from munkres import Munkres, print_matrix
munk = Munkres()
from math import exp
from functions import verif, Q, phi, TreeDistance, TreeBuild 



def optim(Grid,Case,n1,N2,N3,N4,N5,nb_feuilles):

    X_inf = Grid[0][0]
    X_sup = Grid[1][0]
    Delta = np.array(Grid[1], float) - np.array(Grid[0], float)
    
    while sum (np.array(Grid[3], float) < Delta) > 0:  
        M = np.ones(N3)*10**20
        top = np.ones((N3,5))
        Best = np.ones(5)
        for i0 in range (Grid[2][0]+1):
            for i1 in range (Grid[2][1]+1):
                for i2 in range (Grid[2][2]+1):
                    #k = k+ 1
                    #print k
                    for i3 in range (Grid[2][3]+1):
                        for i4 in range (Grid[2][4]+1):
                            Dt = 0
                            V = np.zeros(5)
                            P = [i0, i1, i2, i3, i4]
                            for i in range (5):
                                if Grid[2][i] == 0:
                                    V[i] = Grid[0][i]
                                else:
                                    V[i] = (P[i]*Delta[i]/Grid[2][i])+Grid[0][i]
                            for j in range (N2):    
                                (s,t,time)=TreeBuild(V[0],V[1],V[2],V[3],V[4],nb_feuilles)
                                D = TreeDistance(Case[0],Case[1],s,t.take([0],axis=1),nb_feuilles)     
                                Dt = Dt + D
                                i = np.argmax(M)
                                if Dt < M[i]:
                                    M[i] = Dt
                                    top[i] = V

        Min = 10**20
        for i in range (N3):
            V = top[i]
            Dt = 0
            for j in range (N4):
                (s,t,time)=TreeBuild(V[0],V[1],V[2],V[3],V[4],nb_feuilles)
                D = TreeDistance(Case[0],Case[1],s,t.take([0],axis=1),nb_feuilles)     
                Dt = Dt + D
            if Dt < Min:
                Best = V
                Min = Dt
        
        for i in range (5):
            Grid[0][i] = max(Best[i] - Delta[i]/(2*n1) , Grid[0][i]/2)
            Grid[1][i] = Best[i] + Delta[i]/(2*n1)

        if Grid[0][0] < X_inf:
            Grid[0][0] = X_inf

        if Grid[1][0] > X_sup:
            Grid[1][0] = X_sup

        if verif(Grid[0][0],Grid[1][0],Case[1]) == 0:
            u = abs(Case[1]-Best[0])
            i = np.argmin(u)
            Best[0] = Case[1][i]
            Grid[0][0] = Case[1][i]
            Grid[1][0] = Case[1][i]
            Grid[2][0] = 0 # We stop to find the best X0

        Delta = np.array(Grid[1], float) - np.array(Grid[0], float)

    Score = 0
    V = Best
    for i in range (N5):
        (s,t,time)=TreeBuild(V[0],V[1],V[2],V[3],V[4],nb_feuilles)
        D = TreeDistance(Case[0],Case[1],s,t.take([0],axis=1),nb_feuilles)     
        Score = Score + D
        
    return (Best,Grid, Score)
