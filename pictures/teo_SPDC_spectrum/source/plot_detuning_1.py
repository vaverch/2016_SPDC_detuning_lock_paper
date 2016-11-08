
import Tkinter as tk  # gives tk namespace
import numpy as np
import sys,os
import tkFileDialog
import matplotlib.pyplot as plt
import pdb
from scipy.optimize import curve_fit
import math
import numpy as np
sys.path.append('C:\Users\Gerhard.DESKTOP-M06U1PI\Documents\Libary\python')
from Ghad_lib import smooth 

#FUNCTIONS
#plt.rcParams('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],"size":13})
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True
plt.rcParams['text.latex.preamble']=["\\usepackage{siunitx}","\\usepackage{amsmath}","\\usepackage{amsmath}","\\sisetup{math-rm=\mathsf,text-rm=\sffamily}"]

def exp_sym(x, A, x_off, Tau, y_off):
    return A * np.exp(- (abs(x) - x_off) / Tau) + y_off

def si_rate_d_P0_Pp(gamma_s_MHz, gamma_i_MHz, delta_p, Delta):
    P_th_d_P0 = (1 + 4*delta_p**2) * (1 + Delta**2)
    return 1. / P_th_d_P0
#    return  2*math.pi * gamma_s_MHz * gamma_i_MHz / (gamma_s_MHz + gamma_i_MHz) / P_th_d_P0
    
def exp_asym_trap(x, A, x_l, x_r, Tau_left,Tau_right, y_off):
    ind_l,ind_r  = 0, 0
    x_run = x[0]
    while x_l > x_run:
        x_run = x[ind_l]
        ind_l = ind_l + 1        
    while x_r > x_run:
        x_run = x[ind_r]
        ind_r = ind_r + 1
        
    if ind_l > len(x) - 1:
        ind_l = len(x) - 1
    if ind_r > len(x) - 1:
        ind_r = len(x) - 1
        
    y = np.zeros(len(x))
    
    y[0:ind_l]      =        A * np.exp( - abs(x[0:ind_l]      - x_l) / Tau_left) + y_off
    y[ind_l-1:ind_r]      =    A * np.exp( - abs(x[ind_l]      - x_l) / Tau_left) + y_off + 0*x[ind_l-1:ind_r]
    y[ind_r:len(x)] =        A * np.exp( - abs(x[ind_r:len(x)] - x_r) / Tau_right)+ y_off
    
    return y
    
def exp_asym(x, A, x_off, Tau_left,Tau_right, y_off):
    ind = 0
    x_run = x[0]
    while x_off > x_run:
        x_run = x[ind]
        ind = ind + 1        
    
    if ind > len(x) - 1:
        ind = len(x) - 1
    y = np.zeros(len(x))
    
    y[ind:len(x)] = A * np.exp( - abs(x[ind:len(x)] - x_off) / Tau_right) + y_off
    y[0:ind]      = A * np.exp( - abs(x[0:ind]      - x_off) / Tau_left) + y_off
    return y

def test(x, A, B):
    return A + B * x
    
print 'Go for plot!' 
#plt.rcParams.update({'fontname':'Times New Roman','font.size': 11})
title_font = {'fontname':'Arial', 'size':'8', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
axis_font =  {'fontname':'Times New Roman', 'size':'11'}
ticks_font = {'fontname':'Times New Roman', 'size':'11'}
label_font = {'family':'Times New Roman','size':'11'}
#matplotlib.rc('font', **axis_font)
####################################################
######################THE BRIDGE####################
####################################################
Delta_min = -5
Delta_max = 5
Delta = np.arange(Delta_min,Delta_max,0.01)

####################################################
############################# PLOT  ################
####################################################
fig1 = plt.figure(1,figsize=(8./2.55, 5.5/2.55))
fig1.subplots_adjust(left=0.30)
fig1.subplots_adjust(bottom=0.25)
fig1.subplots_adjust(top=0.80)
fig1.subplots_adjust(right=0.95)
ax1 = fig1.add_subplot(111)

gamma_s_MHz = 1
gamma_i_MHz = 1

delta_p = 0
si_rate_norm = si_rate_d_P0_Pp(gamma_s_MHz, gamma_i_MHz, delta_p,Delta)

ax1_handle = ax1.plot(Delta,si_rate_norm,'r', lw=1)

ax1.set_xlabel('\\textrm{Frequency mismatch} $\Delta$',**label_font)
#ax1.set_ylabel('Pair generation rate $r_{si}$',**label_font)
#ax1.set_ylabel('$\\textrm{Pair production rate}$ \n$\\textrm{ }r_\\textrm{si}\\phantom{.}$ $(2 \pi \\frac{\\gamma_\\textrm{s} \\gamma_\\textrm{i}}{\\gamma_\\textrm{s} + \\gamma_\\textrm{i}} \\frac{P_\\textrm{p}}{P_0}){ aa}$',ha='right',va='center',**label_font)
ax1.set_ylabel('$\\textrm{Pair rate } r_\\textrm{si} \\textrm{ (norm.)}$',**label_font)

    
ax1legend = ax1.legend()
ax1.set_ylim([-0.05,1.05])
ax1.locator_params(axis = 'x', nbins = 10)
ax1.locator_params(axis = 'y', nbins = 4)
ax1.grid()

#fig1.canvas.draw()
#xlabels = [item.get_text() for item in ax1.get_xticklabels()]
#ax1.set_xticklabels(xlabels,**ticks_font)
#ylabels = [item.get_text() for item in ax1.get_yticklabels()]
#ax1.set_yticklabels(ylabels,**ticks_font)


fig1.savefig('SPDC_detuning'+'.png')
fig1.savefig('SPDC_detuning'+'.pdf', transparent=True)

plt.show()
print 'Done'