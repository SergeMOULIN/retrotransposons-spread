# coding: utf-8
import numpy as np
from random import *
from munkres import Munkres
munk = Munkres()
from math import exp
from case_study import Case1, Case2, Case3, dico_genes, prop_cum_taille


def in_gene(liste,position):
    '''
    This function checks whether a given nucleotide position is in a gene or not.
    '''
    dedant = False
    for i in range(len(liste)):
        if liste[i][0] < position and position < liste[i][1]:
            dedant = True
    return dedant


def prop_genes(liste):
    '''
    This function returns the total number of nucleotides constituting the genes of a particular chromosome 
    divided by the total number of nucleotides constituting this chromosome.
    In other words, this function returns the proportion of the given chromosome occupied by genes.
    '''
    sum_nuc = 0.0
    for i in range(len(liste)):
        sum_nuc = sum_nuc + liste[i][1] - liste[i][0]
    return sum_nuc


def nb_ET_in_genes(liste,Case):
    '''
    This function returns the number of transposable element located in genes 
    for a given chromosome and a given family.
    '''
    nb_ET = 0
    for i in range(len(Case[1])):
        if in_gene(liste, Case[1][i]):
            nb_ET += 1
    return nb_ET


def densite (liste,Case):
    '''
    This function calculates the T.E. density in the genes divided by T.E. density outside of genes
    for a given chromosome and a given family.
    This value is used to determine how likely to relaunch the position of a T.E. which moves into a gene during simulations.  
    '''
    nb_ET_in = float(nb_ET_in_genes(liste,Case))
    nb_ET_out = float(len(Case[1]) - nb_ET_in)
    prop_in = prop_genes(liste)
    prop_out = 1 - prop_in
    if (prop_in == 0) or (prop_out == 0):
        densite = 1
    else:
        densite = (nb_ET_in/prop_in) / (nb_ET_out/prop_out)
        #print 'in = ', nb_ET_in
        #print 'surface = ', prop_in
    return (densite)


def dico_densite(dico_genes,Case): 
    '''
    This function apply the "density" function to each chromosome.
    This is useful only to prepare an application of the method to several chromosomes (i.e. The entire genome).
    '''
    dico_dens={}
    for k in dico_genes:
        dico_dens[k] = densite(dico_genes[k],Case)
    return dico_dens


def ponderation(Case):
    '''
    This function returns the mean square deviation of distance divided by 
    the mean square deviation of deterioration. 
    The returned value is used in the calculation of the distance between distance between the final state of a simulated tree 
    and the observed chromosome.
    '''
    moyenne_difference_similitude_carree = 0.0
    moyenne_distance_geographique_carree = 0.0
    for i in range(len(Case[0])-1):
        for j in range(i+1,len(Case[0])):
            moyenne_difference_similitude_carree += (Case[0][i] - Case[0][j])**2
            moyenne_distance_geographique_carree += (Case[1][i] - Case[1][j])**2
    if (moyenne_difference_similitude_carree == 0) or (moyenne_distance_geographique_carree == 0):
            ponde = 1
    else:
        ponde = moyenne_difference_similitude_carree / moyenne_distance_geographique_carree
    return ponde


def listifiage (M):
    '''
    Turn a tuple of tuple into a list of list.
    '''
    N = list(M)
    for i in range(len(M)):
        N[i] = list(M[i])
    return (N) 


def Q(p,s):
    '''
    This function modifies the parameter of the exponential distribution used
    for time before branch's birth as a function of p and the similitude with the "root".
    '''
    y= 1+(1-s)*p
    return (y)


def phi(S1,B,i):
    '''
    #This function changes randomly the similarity of the considered branch between the root when a deterioration 
    #happens.
    '''
    S = S1
    R = 1.0 + expovariate(1.0/B) # The number by which the similarity is divided.
    S[i] = S[i]/R   
    return (S)


def reduction_X(Min, Max, top, n1, division,nb_ET):
    '''
    This function recalculates X0 test interval after each step of the optimization method
    (i.e. after each step of "optim_p2").
    '''
    if Max - Min < division:
        Min2 = top
        Max2 = top
    else:
        diff = Max - Min
        ecart = float(diff) / (2*n1)
        Min2 = max(round(top - ecart),0)
        Max2 = min(round(top + ecart),nb_ET - 1)
    return Min2, Max2


def Ballayage_X(Min, Max, division):
    '''
    This function defined in which points X0 have to be tested at each stage of the optimization method
     (i.e. after each step of "optim_p2") depending on the test interval.
    '''
    if Max - Min < division:
        liste2 = range(int(Min), int(Max+1))
    else:
        liste = []
        diff = Max - Min
        ecart = float(diff) / division
        for i in range(division+1):
            
            liste.append(Min + round(i * ecart))
        liste2 =  list(set(liste))   # Remove redundancies.
    return (liste2)



def duplication(chromosome,X,Lx,a,dico_genes,prop_cum_taille,dico_dens):
    '''
    This function determines in which position the T.E. moves when duplicating.
    '''
    u1 = uniform(0,1)
    if u1 < a:    # For a uniform distribution of the daughter copy throughout the genome. This is useful only to prepare an application of the method to the entire genome.
        res = 0  
        while res == 0:
            u2 = uniform(0,1)
            v = np.ones(len(prop_cum_taille[0]))*u2
            i = sum(prop_cum_taille[1] < v) 
            chrom_fille = prop_cum_taille[0][i]
            X1 = uniform(0,1)
            liste = dico_genes[chrom_fille]
            if (not(in_gene(dico_genes[chrom_fille],X1)) or (uniform(0,1) < dico_dens[chrom_fille])):
                res = 1
        return chrom_fille,X1 
    else:         # Case of exponential distribution of the daughter copy near the mother.   
        res = 0  
        while res == 0:
            Xd=randint(0,1) # Xd indicates the direction of displacement.
            if Xd == 0: 
                X1=X-expovariate(1.0/Lx) 
            else:
                X1=X+expovariate(1.0/Lx)
            if (0 < X1 < 1 ) and (not(in_gene(dico_genes[chromosome],X1)) or (uniform(0,1) < dico_dens[chromosome])):
                res = 1
        return chromosome,X1 
             

def moyenne_saut_reelle(T):
    '''
    Calculates the mean traveled distance in simulations (i.e. "J" value in the article).
    '''
    D = 0
    for i in range(1,T.shape[0]):
        mother = T[i,2]
        D += abs(T[i,0]-T[mother,0])
    D = float(D) / (T.shape[0] - 1)
    return D

def DistriDistance(S1,S2):
    S1 = sorted(S1)
    S2 = sorted(S2)
    D = 0
    for i in range (len(S1)):
        D += (S1[i] - S2[i])**2
    return(D)


def TreeDistance(S1,T1,S2,T2,nb_feuilles,ponderation):
    '''
    This function computes the distance between two trees using the munkres
    algorithm made by Yi Cao with weights as functions of the distances
    between branches (in the nucleotids' sequences and in the graph of the
    tree)
    There is also a penality as function of the difference between the number
    of branches for each tree
    '''
    n = m = nb_feuilles
    MC=np.ones((n,m))
    for i in range (n):
        for j in range (m):
            MC[i,j]=ponderation*((S1[i]-S2[j])**2)+(T1[i]-T2[j])**2
    MC = MC.tolist()
    D = sum([MC[couple[0]][couple[1]] for couple in munk.compute(MC)])
    return (D)


def DistriBuild(mu,B,p,nb_feuilles):
    '''
    Simulate a distribution of deterioration state.
    '''
    S=np.ones(1)                     # Vector of similarities to the root.
    inten_cum = np.ones(2) * (1./mu) # The vector of cumulative intensities.
    inten_cum[1] = (1./mu) + 1       
    time = 0
    
    while (len (S) < nb_feuilles):  # The stop condition.
        n = len(S)                   
        time = time + expovariate(inten_cum[-1])  # Time of the next event.  
        u = uniform(0,inten_cum[-1])              
        i = sum(inten_cum < u)                    # Location and type of the next event.
        
        if i >= n:                                # Case of a duplication.
            i = i - n
            S = np.concatenate((S,[S[i]]),axis=0) 
            inten_cum = inten_cum + (1./mu)
            inten_cum = np.concatenate (([1./mu],inten_cum,[1./Q(p,S[-1])+inten_cum[-1]]),axis=0)
        else:
            prov = 1./Q(p,S[i])
            S = phi(S,B,i)
            for k in range (n + i , 2*n):
                inten_cum[k] = inten_cum[k] + (1./Q(p,S[i]) - prov)
            if S[i] < 0.03:
                for k in range (i , 2*n):
                    inten_cum[k] = inten_cum[k] - (1./mu)
                    
    BS = [S]
    Btime = [time]
    res = 0
    
    while res == 0:
        n = len(S)
        time = time + expovariate(inten_cum[-1]) 
        u = uniform(0,inten_cum[-1])
        i = sum(inten_cum < u)
        
        if i >= n:  # Case of a duplication.
            res = 1

        else:
            prov = 1./Q(p,S[i])
            S = phi(S,B,i)
            for k in range (n + i , 2*n):
                inten_cum[k] = inten_cum[k] + (1./Q(p,S[i]) - prov)
            if S[i] < 0.03:
                for k in range (i , 2*n):
                    inten_cum[k] = inten_cum[k] - (1./mu)
            BS.append(S)
            Btime.append(time)
            
    time = float(Btime[0]+time)/2
    i = sum(Btime <= np.ones(len(Btime))*time)-1
    S = BS[i]
    return (S,time)


def TreeBuild(chrom,indiceX,mu,B,p,Lx,a,nb_feuilles,dico_genes,prop_cum_taille,Case,dico_dens):
    '''
    This function simulates the spread of the T.E.
    '''
    Positions = sorted(Case[1])
    T = np.zeros((1,3))  # Matrix whose columns are the position of the T.E., the moment of his birth and the index of the parent branch.
    S = np.ones(1)       # Vector of the similarities to the root.
    C = []               # List of chromosomes on which the copies are positioned. This is useful only to prepare an application of the method to the entire genome).
    inten_cum = np.ones(2) * (1./mu) # The vector of cumulative intensity.
    inten_cum[1] = (1./mu) + 1   
    T[0,0]=Positions[int(indiceX)] # Copies positions. The root is necessarily placed on an existing copy.
    T[0,1]=0                  # The moment of the T.E. birth.
    T[0,2]=0                  # The index of the parent branch.
    C.append(chrom)  
    time = 0
    
    while (len (S) < nb_feuilles):
        n = len(S)
        time = time + expovariate(inten_cum[-1]) 
        u = uniform(0,inten_cum[-1])
        i = sum(inten_cum < u)
        
        if i >= n:  # Case of a duplication.
            i = i - n
            chrom_fille, X1 = duplication(C[i],T[i,0],Lx,a,dico_genes,prop_cum_taille,dico_dens)
            
            new_branche = [[X1, time, i]]
            T = np.concatenate((T,new_branche),axis=0) 
            S = np.concatenate((S,[S[i]]),axis=0) 
            C.append(chrom_fille)
            inten_cum = inten_cum + (1./mu)
            inten_cum = np.concatenate (([1./mu],inten_cum,[1./Q(p,S[-1])+inten_cum[-1]]),axis=0)
        else:
            prov = 1./Q(p,S[i])
            S = phi(S,B,i)
            for k in range (n + i , 2*n):
                inten_cum[k] = inten_cum[k] + (1./Q(p,S[i]) - prov)
            if S[i] < 0.03:
                for k in range (i , 2*n):
                    inten_cum[k] = inten_cum[k] - (1./mu)
                    
    BS = [S]
    BT = [T]
    Btime = [time]
    res = 0

    while res == 0:
        n = len(S)
        time = time + expovariate(inten_cum[-1]) 
        u = uniform(0,inten_cum[-1])
        i = sum(inten_cum < u)
        
        if i >= n:  #Case of a duplication.
            res = 1

        else:
            prov = 1./Q(p,S[i])
            S = phi(S,B,i)
            for k in range (n + i , 2*n):
                inten_cum[k] = inten_cum[k] + (1./Q(p,S[i]) - prov)
            if S[i] < 0.03:
                for k in range (i , 2*n):
                    inten_cum[k] = inten_cum[k] - (1./mu)
            BS.append(S)
            BT.append(T)
            Btime.append(time)
            
    time = uniform(Btime[0],time)
    i = sum(Btime <= np.ones(len(Btime))*time)-1
    S = BS[i]
    T = BT[i]
    
    return (S,T,time)

