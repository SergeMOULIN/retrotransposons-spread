# coding: utf-8 
import numpy as np
from random import *
from munkres import Munkres, print_matrix
munk = Munkres()
from math import exp
from functions import listifiage, Q, phi, TreeDistance, TreeBuild, TreeDraw 
from optim import optim
from case_study import Case1, Case2, Case3
import matplotlib.pyplot as plt
#from sys import argv


def f1(i):
    Case = Case1
    nb_feuilles = len(Case[0])
    X_min = min(Case[1])
    X_max = max(Case[1])
    Grid_1=((X_min, 0.1, 0.01, 0.1, 0.01),  
            (X_max,10.,1.,100.,1.),
            (3,4,3,4,3),
            (10**-3,10**-2,10**-3,10**-1,10**-3))
                      
    Grid = listifiage(Grid_1)
    
    n1 =1.5
    N2 = 100  
    N3 = 100 
    N4 = 1000
    N5 = 20000

    (B,G,S) = optim(Grid, Case, n1, N2 , N3, N4 , N5, nb_feuilles)
    print "Job =", i, "Case =", 1,  "Best =", B, "Grille =", G, "Score =", S   
    return()


def f2(i):
    Case = Case1
    nb_feuilles = len(Case[0])
    X_min = min(Case[1])
    X_max = max(Case[1])
    Grid_1=((X_min, 0.1, 0.01, 0.1, 0.01),  
            (X_max,10.,1.,100.,1.),
            (3,4,3,4,3),
            (10**-3,10**-2,10**-3,10**-1,10**-3))
                      
    Grid = listifiage(Grid_1)
    
    n1 =1.5
    N2 = 100  
    N3 = 100 
    N4 = 1000
    N5 = 20000

    (B,G,S) = optim(Grid, Case, n1, N2 , N3, N4 , N5, nb_feuilles)
    print "Job =", i, "Case =", 1,  "Best =", B, "Grille =", G, "Score =", S   
    return()



def f3(i):
    Case = Case1
    nb_feuilles = len(Case[0])
    X_min = min(Case[1])
    X_max = max(Case[1])
    Grid_1=((X_min, 0.1, 0.01, 0.1, 0.01),  
            (X_max,10.,1.,100.,1.),
            (3,4,3,4,3),
            (10**-3,10**-2,10**-3,10**-1,10**-3))
                      
    Grid = listifiage(Grid_1)
    
    n1 =1.5
    N2 = 100  
    N3 = 100 
    N4 = 1000
    N5 = 20000

    (B,G,S) = optim(Grid, Case, n1, N2 , N3, N4 , N5, nb_feuilles)
    print "Job =", i, "Case =", 1,  "Best =", B, "Grille =", G, "Score =", S   
    return()



def f(i):

    i = i-1
    
    if i in range (0,60):
        f1(i)

    if i in range (60,120):
        f2(i)

    if i in range (120,180):
        f3(i)


#f1(1) # to applie Case 1 juste once, for instance.
#f(argv[-1]) # In case you use a .sge code.


'''
M = 10**20
S1 = []
T1 = []
for k in range (1000000):
    (s,t,n) = TreeBuild(0,B[i][0],0,B[i][1],B[i][2],B[i][3],B[i][4],B[i][5],0,1,Lim_rec,Lim_bra)
    if n == 32:
        D=TreeDistance(Case[0],Case[1],s,t.take([0],axis=1),Pen,Lim_mun) 
        if D < M:
            M = D
            S1 = s
            T1 = t
    
print "S1 = ", S1
print "T1 = ", T1
print "end"
'''


