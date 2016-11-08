
import Tkinter as tk  # gives tk namespace
import numpy as np
import sys,os
import tkFileDialog
import matplotlib.pyplot as plt
import pdb
from scipy.optimize import curve_fit
import math
import numpy as np
sys.path.append('C:\Users\gschunk.MPL\Documents\Libary\python')
from Ghad_lib import smooth, lorentz
#FUNCTIONS
#plt.rcParams('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],"size":13})

g_use_tex = True
if g_use_tex == True:
  plt.rcParams['text.usetex'] = True
  plt.rcParams['text.latex.unicode'] = True
  plt.rcParams['text.latex.preamble']=["\\usepackage{siunitx}","\\usepackage{amsmath}","\\usepackage{amsmath}","\\sisetup{math-rm=\mathsf,text-rm=\sffamily}"]

def exp_sym(x, A, x_off, Tau, y_off):
    return A * np.exp(- (abs(x) - x_off) / Tau) + y_off

def Lorentz(gamma_s_MHz, gamma_i_MHz, delta_p, Delta):
    P_th_d_P0 = (1 + 4*delta_p**2) * (1 + Delta**2)
    return 1. / P_th_d_P0
    
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
plt.rcParams.update({'fontname':'Times New Roman','font.size': 11})
title_font = {'fontname':'Arial', 'size':'8', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
axis_font =  {'fontname':'Times New Roman', 'size':'11'}
ticks_font = {'fontname':'Times New Roman', 'size':'11'}
label_font = {'family':'Times New Roman','size':'11'}
legend_font = {'family':'Times New Roman','size':'9'}
#matplotlib.rc('font', **axis_font)
####################################################
######################THE BRIDGE####################
####################################################
Delta_min = -4
Delta_max = 8
DeltaDelta = 0.01
Delta = np.arange(Delta_min,Delta_max,DeltaDelta)

####################################################
############################# PLOT  ################
####################################################
gamma_s_MHz = 1
gamma_i_MHz = 1

delta_p = 0
Delta_idler_list = [10]
for Delta_idler in Delta_idler_list:
#Delta_idler = 3
  fig1 = plt.figure(1,figsize=(5./2.55, 4.3/2.55))
  fig1.subplots_adjust(left=0.30)
  fig1.subplots_adjust(bottom=0.25)
  fig1.subplots_adjust(top=0.90)
  fig1.subplots_adjust(right=0.95)
  ax1 = fig1.add_subplot(111)

  lor_signal = lorentz(Delta,0,gamma_s_MHz,1,0)
  lor_idler = lorentz(Delta,Delta_idler,gamma_s_MHz,1,0)
  lor_idler_signal = lor_signal*lor_idler
  
  #Normalize overlap
  integral = 0
  for intint in range(len(lor_idler_signal)):
    integral = integral + DeltaDelta*lor_idler_signal[intint]
  
  lor_idler_signal = lor_idler_signal / integral
  
  ax1_handle = ax1.plot(Delta,lor_signal,'-r', lw=0.8,label = '$\\textrm{signal}$')
  ax2_handle = ax1.plot(Delta,lor_idler,'-b', lw=0.8,label = '$\\textrm{shifted}$ \n $\\textrm{idler}$')
  ax3_handle = ax1.plot(Delta,lor_idler_signal,':k', lw=1.5,label = '$\\textrm{product}$')
  ax4_handle = ax1.plot(Delta,lor_idler_signal,'-k', lw=1.5,label = '$\\textrm{normal.}$ \n $\\textrm{product}$')
  
#   \n $\\textrm{idler}$'

  if g_use_tex == True:
#    ax1.set_xlabel('\\textrm{Detuning} $(\\nu - \\nu (\\ell_\\textrm{s},\\textrm{q}_\\textrm{s},\\textrm{p}_\\textrm{s}))$',**label_font)
    ax1.set_xlabel('\\textrm{Detuning} $\\frac{\\nu - \\nu (\\ell_\\textrm{s},\\textrm{q}_\\textrm{s},\\textrm{p}_\\textrm{s})}{(\\gamma_\\textrm{s}+\\gamma_\\textrm{i})/2}$',**label_font)
    if Delta_idler == 0:
#      ax1.set_ylabel('$\\textrm{F}_\\textrm{s,i} (\\nu) $',ha='right',va='center',**label_font)
      ax1.set_ylabel('$\\textrm{Amplitude}$',ha='right',va='center',**label_font)
  else:
    ax1.set_xlabel('Frequency detuning (nu - nu_s)',**label_font)
    ax1.set_ylabel('F (nu)',ha='right',va='center',**label_font)
  
  ax1.set_ylim([-0.05,1.05])
  ax1.locator_params(axis = 'x', nbins = 6)
  ax1.locator_params(axis = 'y', nbins = 4)
  ax1.grid()
  
  if g_use_tex == False:
    fig1.canvas.draw()
    xlabels = [item.get_text() for item in ax1.get_xticklabels()]
    ax1.set_xticklabels(xlabels,**ticks_font)
    ylabels = [item.get_text() for item in ax1.get_yticklabels()]
    ax1.set_yticklabels(ylabels,**ticks_font)
  
  if Delta_idler == 10:
    ax1legend = ax1.legend()
    for label in ax1legend.get_texts():
      ax1.legend(loc='upper right',prop=legend_font)
  
  fig1.savefig('SPDC_overlap_Delta_'+ str(Delta_idler)+'_.png')
  fig1.savefig('SPDC_overlap_Delta_'+ str(Delta_idler)+'_.pdf', transparent=True)

  plt.clf()
#plt.show()
print 'Done'