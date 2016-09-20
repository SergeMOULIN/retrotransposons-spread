# coding: utf-8 
import numpy as np
import time
from random import *
from munkres import Munkres, print_matrix
munk = Munkres()
from math import exp
from functions import listifiage, DistriDistance, TreeDistance, DistriBuild, TreeBuild
from method import optim
from case_study import Case1, Case2, Case3, dico_genes, prop_cum_taille
from sys import argv
#from multiprocessing import Process



def f1(i):
    Case = Case1
    nb_feuilles = len(Case[0])
    X_min = min(Case[1])
    X_max = max(Case[1])
    Grid_1=((0,               0.1,  0.01,   0.1,    0.01,   0.),  
            (len(Case[1])-1,  10.,    1.,   100.,   1.,     0.),
            (3,3,3,3,3,0),
            (10**-3,10**-2,10**-3,10**-1,10**-3,0.))
                      
    Grid = listifiage(Grid_1)
    
    n1 =1.5
    N2 = 0
    N3 = 100
    N4 = 10000
    N5 = 20000
    
    (B,G,S,M,T) = optim(Grid, Case, n1, N2 , N3, N4 , N5, nb_feuilles,dico_genes,prop_cum_taille)
    print "Job =", i, "Case =", 1,  "Best =", B, "Grille =", G, "Score =", S, "moyenne saut =", M, 'T_obs_moyen =', T   
    return()


def f2(i):
    Case = Case2
    nb_feuilles = len(Case[0])
    X_min = min(Case[1])
    X_max = max(Case[1])
    Grid_1=((0,               0.1,  0.01,   0.1,    0.01,   0.),  
            (len(Case[1])-1,  10.,    1.,   100.,   1.,     0.),
            (3,3,3,3,3,0),
            (10**-3,10**-2,10**-3,10**-1,10**-3,0.))

    Grid = listifiage(Grid_1)

    n1 =1.5
    N2 = 0
    N3 = 100
    N4 = 10000
    N5 = 20000
    
    (B,G,S,M,T) = optim(Grid, Case, n1, N2 , N3, N4 , N5, nb_feuilles,dico_genes,prop_cum_taille)
    print "Job =", i, "Case =", 2,  "Best =", B, "Grille =", G, "Score =", S, "moyenne saut =", M, 'T_obs_moyen =', T   
    return()


def f3(i):
    Case = Case3
    nb_feuilles = len(Case[0])
    X_min = min(Case[1])
    X_max = max(Case[1])
    Grid_1=((0,               0.1,  0.01,   0.1,    0.01,   0.),  
            (len(Case[1])-1,  10.,    1.,   100.,   1.,     0.),
            (3,3,3,3,3,0),
            (10**-3,10**-2,10**-3,10**-1,10**-3,0.))
    
    Grid = listifiage(Grid_1)
    
    n1 =1.5
    N2 = 0
    N3 = 100
    N4 = 10000
    N5 = 20000
    
    (B,G,S,M,T) = optim(Grid, Case, n1, N2 , N3, N4 , N5, nb_feuilles,dico_genes,prop_cum_taille)
    print "Job =", i, "Case =", 3,  "Best =", B, "Grille =", G, "Score =", S, "moyenne saut =", M, 'T_obs_moyen =', T   
    return()


def f(i):

    i = i-1
    if i in range (0,20):
        f1(i)

    if i in range (20,40):
        f2(i)

    if i in range (40,60):
        f3(i)

    if i in range (60,80):
        f1(i)

    if i in range (80,100):
        f2(i)

    if i in range (100,120):
        f3(i)

'''
import profile
profile.run('f(1)','f.profile')
import pstats
pstats.Stats('f.profile').sort_stats('time').print_stats()
'''

f(1)
#f(eval(argv[-1]))
