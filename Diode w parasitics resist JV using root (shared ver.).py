# -*- coding: utf-8 -*-
# ROTINA PARA RESOLVER NUMERICAMENTE A EQ. DE DIODO SCHOTTKY NA PRESENCA DE RESISTENCIAS PARASITAS
# ALEXANDRE DE CASTRO MACIEL, DEP. DE FISICA DA UFPI
# OUTUBRO DE 2015

import matplotlib.pyplot as plt
import numpy as np
import __main__
from scipy import optimize
from matplotlib import font_manager
    
def fun(i,v): # EQUACAO TRANSCENDENTAL NO FORMATO F(i,v)=0
    Rs    = __main__.Rs
    Rp    = __main__.Rp    
    i0    = __main__.i0 
    m     = __main__.m
    kT    = __main__.kT
    return (1.0+Rs/Rp)*i-v/Rp-i0*(np.exp((v-Rs*i)/(m*kT))-1.0)
    
def findroot(vv): # FUNCAO PARA ENCONTRAR RAIZ DE fun(i,v)
    i = 0.0 #GUESS PARA A CORRENTE
    ii_teo  = []
    for v in vv:
        sol = optimize.root(fun,i,args=(v)) #ROTINA PARA ENCONTRAR RAIZES
        i=sol.x # VALOR DE CORRENTE RAIZ DA FUNCAO fun(i,v) PARA v
        ii_teo.append(i)
    return ii_teo

# PARAMETROS FIXOS
T        = 330 #in K
k        = 8.617e-5 #in eV/K
kT       = k*T
# PARAMETROS VARIÁVEIS
Rs       = 4.0e2 # in Ohms
Rp       = 8.5e7 # in Ohms
m        = 1.35
i0       = 1.0e-10 # in Amperes
# INICIALIZANDO VETORES
vv  = []
ii_exp  = []
ii_exp_abs  = []
ii_teo_abs  = []
# 
directory = "diretorio/para/seu/dado/experimental/" # Directory of files. End with /
path = directory+"nome_do_arquivo.csv" #arquivo em duas colunas separadas por virgula.
look = 0
with open(path, 'r') as data:
    for line in data.readlines():
        if(look):
            vv.append(float(line.split(",")[0]))
            ii_exp.append((float(line.split(",")[1])))
            ii_exp_abs.append(abs(float(line.split(",")[1])))
        if(line[0:7] == "Voltage"):
            look = 1
            
ii_teo = findroot(vv)
for i in ii_teo:
    ii_teo_abs.append(abs(i))

ran = 0
vv_p=vv[ran:-1-ran]
ii_exp_p= ii_exp[ran:-1-ran]
ii_exp_abs_p= ii_exp_abs[ran:-1-ran]
ii_teo_p= ii_teo[ran:-1-ran]
ii_teo_abs_p= ii_teo_abs[ran:-1-ran]

font = font_manager.FontProperties(family ='Purisa')
fontx = {'fontname':'Purisa'}

plt.clf()
plt.xkcd()
plt.axvline(0, color='black')
plt.xlabel('Voltagem (V)', **fontx)
plt.ylabel('Corrente (A)', **fontx)
plt.title('Grafico 5', **fontx)
plt.text(0.3, 1e-8, r'$\alpha$ = $1.35$ ',fontsize=18)
plt.text(0.3, 1e-9, r'$R_P$ = $8\times 10^7 \Omega$',fontsize=18)
plt.text(0.3, 1e-10, r'$R_S$ = $400 \Omega$',fontsize=18)
plt.semilogy(vv_p , ii_exp_abs_p , 'g-o', label = 'Dados Experimentais')
plt.semilogy(vv_p , ii_teo_abs_p , 'r-', label = 'Ajuste teorico: Diodo', alpha=.5)
plt.legend(loc='upper left', shadow=True, prop = font)
plt.show()