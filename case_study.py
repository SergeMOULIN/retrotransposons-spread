# coding: utf-8 
import numpy as np

dd=eval(open('dicoInfosETs.txt').read()) 


Case1 = np.zeros((2,32))
cpt = -1
for ET in dd.keys() : 
    if 'ROO' == dd[ET]['nom'].upper() and dd[ET]['chrom'] == "3L":
        cpt = cpt + 1
        Case1[0][cpt] = (100 - dd[ET]['simil'])/100
        Case1[1][cpt] = (eval(dd[ET]['position'][0])+eval(dd[ET]['position'][1]))*(1.0/(2*32600000))


Case2 = np.zeros((2,16)) 
cpt = -1
for ET in dd.keys() : 
    if 'DM412' == dd[ET]['nom'].upper() and dd[ET]['chrom'] == "3L":
        cpt = cpt + 1
        Case2[0][cpt] = (100 - dd[ET]['simil'])/100
        Case2[1][cpt] = (eval(dd[ET]['position'][0])+eval(dd[ET]['position'][1]))*(1.0/(2*32600000))


Case3 = np.zeros((2,6)) 
cpt = -1
for ET in dd.keys():
    if 'GYPSY' == dd[ET]['nom'].upper() and dd[ET]['chrom'] == "3L":
        cpt = cpt + 1
        Case3[0][cpt] = (100 - dd[ET]['simil'])/100
        Case3[1][cpt] = (eval(dd[ET]['position'][0])+eval(dd[ET]['position'][1]))*(1.0/(2*32600000))

del (dd)


