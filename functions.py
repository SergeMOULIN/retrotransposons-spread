# coding: utf-8
import numpy as np
from random import *
from munkres import Munkres, print_matrix
munk = Munkres()
from math import exp
import matplotlib.pyplot as plt


def verif (a,b,case):
    '''
    This function checks if there is at least one element of case study on X_0's interval.
    '''
    res = 0
    for i in range (len(case)):
        if a < case[i] < b:
            res = 1
    return (res)


def listifiage (M):
    '''
    This function turns a tuple into a list.
    '''
    N = list(M)
    for i in range(len(M)):
        N[i] = list(M[i])
    return (N) 


def Q(p,s):
    '''
    This function modifies the parameter of the exponential distribution used
    for time before branch's birth as a function of p and the distance between the "root".
    '''
    y= 1+s*p
    return (y)
    


def phi(S1,B,i):
    S = S1
    '''
    This function changes randomly the distance of the considered branch between the root when a mutation
    happens.
    '''
    R = min(expovariate(1.0/B),1)
    # the proportion of genome that mutate (limited to 1)
    S[i] = R*(1-S1[i])+S[i]-R*S1[i]/3.    #   L' idee ici etant que 1/3 des nucleotides 
                                           #differents de la racine et qui muttent redeviennent comme la racine.
    return (S)




def TreeDistance(S1,T1,S2,T2,nb_feuilles):
    '''
    This function computes the distance between two trees using the munkres
    algorithm.
    '''
    n = m = nb_feuilles
    MC=np.ones((n,m))
    
    for i in range (n):
        for j in range (m):
            MC[i][j]=abs(S1[i]-S2[j])+abs(T1[i]-T2[j])
                
    MC = MC.tolist()
    D = sum([MC[i[0]][i[1]] for i in munk.compute(MC)])

    return (D)




def TreeBuild(X,mu,B,p,L,nb_feuilles):
    '''
    Treebuild is a recursive function creating a tree as a matrix
    
    S is the vector of NW distances to the root. 
    X is the root's position in the genome. 
    mu is the parameter of the exponential distribution used for time before 
    mutation.
    B the parameter of the exponential distribution used for the impact of mutations. 
    L the parameter of the exponential distribution used for the distance where the child is created.
    p the parameter that determine how much mutations led to a decrease in the duplication speed.
    i is the distance to the root measured in number of branches. Example : i = 0 for the root, i = 1 for children of the root, i = 2 for grandchildren ect... 
    n is the number of branches already builded. 
    inten_cum is the vector of cumulated rate parameters. The first n rate parameters are these of mutation speeds and the last n 
    rate parameters are these of duplication speeds. 
    '''
    T = np.zeros((1,3)) # The matrix where the columns are positions, time, and index of the mother branches.
    S=np.zeros(1) # Vector of NW distances to the root.
    inten_cum = np.ones(2) * (1./mu) #  Vector of cumulated rate parameters. 
    inten_cum[1] = (1./mu) + 1       #  Vector of cumulated rate parameters.
    T[0][0]=X # Position.
    T[0][1]=0 # Time of birth.
    T[0][2]=0 # Mother.
    time = 0
    
    while (len (S) < nb_feuilles):
        n = len(S)
        time = time + expovariate(inten_cum[-1]) 
        u = uniform(0,inten_cum[-1])
        i = sum(inten_cum < u)
        
        if i >= n:  # Case of transposition
            i = i - n
            X = T[i][0]
            res = 0
            while res == 0:
                Xd=randint(0,1) # Xd indicates the direction of the move.
                if Xd == 0: 
                    X1=X-expovariate(1/L) 
                else:
                    X1=X+expovariate(1/L)
                if (0 < X1 < 1 ):
                    res = 1
                    X2 = X1

            new_branche = [[X2, time, i]]
            T = np.concatenate((T,new_branche),axis=0) 
            S = np.concatenate((S,[S[i]]),axis=0) 
            inten_cum = inten_cum + (1/mu)
            inten_cum = np.concatenate (([1/mu],inten_cum,[1./Q(p,S[-1])+inten_cum[-1]]),axis=0) 
        else:
            prov = 1./Q(p,S[i])
            S = phi(S,B,i)
            for k in range (n + i , 2*n):
                inten_cum[k] = inten_cum[k] + (1./Q(p,S[i]) - prov)
            if 0.74 < S[i] < 0.76:
                for k in range (i , 2*n):
                    inten_cum[k] = inten_cum[k] - (1/mu)
                    
    BS = [S]
    BT = [T]
    Btime = [time]
    res = 0

    while res == 0:
        n = len(S)
        time = time + expovariate(inten_cum[-1]) 
        u = uniform(0,inten_cum[-1])
        i = sum(inten_cum < u)
        
        if i >= n:  # Case of a transposition
            res = 1

        else:
            prov = 1./Q(p,S[i])
            S = phi(S,B,i)
            for k in range (n + i , 2*n):
                inten_cum[k] = inten_cum[k] + (1./Q(p,S[i]) - prov)
            if 0.74 < S[i] < 0.76:
                for k in range (i , 2*n):
                    inten_cum[k] = inten_cum[k] - (1/mu)
            BS.append(S)
            BT.append(T)
            Btime.append(time)
            

    time = uniform(Btime[0],time)
    i = sum(Btime <= np.ones(len(Btime))*time)-1
    S = BS[i]
    T = BT[i]
    
    return (S,T,time)

 
def TreeDraw(T,time):
    '''
    TreeDraw, is the function that allows to draw graphical representations of the spread of transposable elements.
    '''
    fig = plt.figure()
    fig.patch.set_facecolor('w')
    for i in range(len(T)):
        x = [T[i][0],T[i][0]]
        y = [-T[i][1],-time]
        plt.plot(x,y,color = 'k',markeredgecolor='w',markerfacecolor='w')     # Draw the vertical lines.
    for i in range (1,len(T)):
        mother = int(T[i][2])
        x = [T[mother][0],T[i][0]]
        y = [-T[i][1],-T[i][1]]
        plt.plot(x,y,color = 'k',markeredgecolor='w',markerfacecolor='w')     # Draw the horizontal lines.
    
    delta_x = max(T.take(0,axis=1)) - min(T.take(0,axis=1))
    plt.xlim(min(T.take(0,axis=1))-0.1*delta_x,max(T.take(0,axis=1)+0.1*delta_x))
    delta_y = max(T.take(1,axis=1)) - min(T.take(1,axis=1))
    plt.ylim(-time-0.1*delta_y, min(T.take(1,axis=1)+0.1*delta_y))
    plt.xlabel('Location in the genome',fontsize = 14)
    #plt.xlabel('Location in the chromosome',fontsize = 14)
    plt.ylabel('Time',fontsize = 14)
    plt.show()


