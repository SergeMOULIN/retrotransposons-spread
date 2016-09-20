# coding: utf-8 
import numpy as np


dd=eval(open('similIdent20.txt').read()) # This file contains the position and distance to the root of each T.E. 
tailles = {'3L': 24543557}               # Size of the "3L" chromosome.

# The following chapter is useful only to prepare an application of the method to several chromosomes (i.e. The entire genome).
taille_tot = sum(tailles[k] for k in tailles)
liste_chrom = []
liste_cum_taille = []
cum = 0
for k in tailles:
    cum += float(tailles[k])/taille_tot
    liste_chrom.append(k)
    liste_cum_taille.append(cum)    
prop_cum_taille = [liste_chrom, liste_cum_taille]
# End of the chapter.


def make_liste(nom_chromosome):
    '''
    This function produces a list of couple. Each couple represents the position of the beginning and the end of a gene. 
    '''
    nom_input = "genes/dmel-" + nom_chromosome + "-gene-r5.51.fasta"
    dd_genes = open(nom_input).read().split('\n>')
    liste=[]
    for k in dd_genes:
        if "complement(" in k:
            a = eval(k.split("complement(")[1].split('..')[0])
            b = eval(k.split("..")[1].split(')')[0])
        else:
            a = eval(k.split('loc=')[1].split(':')[1].split('..')[0])
            b = eval(k.split('loc=')[1].split('..')[1].split(';')[0])
        a = float(a) / tailles[nom_chromosome]
        b = float(b) / tailles[nom_chromosome]
        liste.append([a,b])
    return (liste)

dico_genes = {}
for chromosome in tailles:
    dico_genes[chromosome] = make_liste(chromosome) 

def sorting (Case):
    '''
    This function sorts T.E. according to their positions on the chromosome.
    '''
    n = Case.shape[1]
    NC = []  # New Case
    for i in range (n):
        NC.append([Case[1][i],Case[0][i]])
    NC = sorted(NC)
    NNC = np.zeros((2,n)) 
    for i in range (n):
        NNC[0][i] = NC[i][1]
        NNC[1][i] = NC[i][0]
    return(NNC)


Case1 = np.zeros((3,16)) 
cpt = -1
for ET in dd.keys() : 
    if 'DM412' == dd[ET]['nom'].upper() and dd[ET]['chrom'] == "3L":
        cpt = cpt + 1L
        Case1[0][cpt] = (dd[ET]['simil'])/100.0
        Case1[1][cpt] = (eval(dd[ET]['position'][0])+eval(dd[ET]['position'][1]))*(1.0/(2*tailles["3L"]))
        

Case2 = np.zeros((3,6)) 
cpt = -1
for ET in dd.keys():
    if 'GYPSY' == dd[ET]['nom'].upper() and dd[ET]['chrom'] == "3L":
        cpt = cpt + 1
        Case2[0][cpt] = (dd[ET]['simil'])/100.0
        Case2[1][cpt] = (eval(dd[ET]['position'][0])+eval(dd[ET]['position'][1]))*(1.0/(2*tailles["3L"]))
        

Case3 = np.zeros((3,32))
cpt = -1
for ET in dd.keys() : 
    if 'ROO' == dd[ET]['nom'].upper() and dd[ET]['chrom'] == "3L":
        cpt = cpt + 1
        Case3[0][cpt] = (dd[ET]['simil'])/100.0
        Case3[1][cpt] = (eval(dd[ET]['position'][0])+eval(dd[ET]['position'][1]))*(1.0/(2*tailles["3L"]))   
        print ET, Case3[0][cpt], Case3[1][cpt]


Case1 = sorting (Case1)
Case2 = sorting (Case2)
Case3 = sorting (Case3)
