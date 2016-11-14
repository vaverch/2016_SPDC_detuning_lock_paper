# -*- coding: utf-8 -*-
"""
Plot mode spectrum nicely
"""
###############Initialize################
import matplotlib.pyplot as plt
import sys, os
import numpy as np
from scipy.signal import savgol_filter
sys.path.append('C:\\Users\\Gerhard.DESKTOP-M06U1PI\\Documents\\Libary')
os.system('CLS')

print 'Welcome, my Master!', '\n'
###############Define Folders################
pathname = os.getcwd()


exp_folder = 'pumpmode'
plots_folder = exp_folder
exp_path = pathname + '\\' + exp_folder
plots_path = pathname + '\\' + plots_folder

print 'experimental data folder: ',exp_folder,'\n'
print 'plots folder: ',plots_folder,'\n'


if not os.path.exists(exp_path):
	os.makedirs(exp_path)
if not os.path.exists(plots_path):
	os.makedirs(plots_path)

###############READ DATA################
filename = 'pumpsquareunder12015.04.24_17.39_06'
#f_x = open(exp_folder + '\\' + filename + '_x.dat',"w")
#f_y = open(exp_folder + '\\' + filename + '_y.dat',"w")
data = np.genfromtxt(exp_folder + '\\' + filename + '.dat', delimiter='')
x = np.zeros(len(data))
y = np.zeros(len(data))
for ite in range(len(data)):
	x[ite] = data[ite,0]
	y[ite] = data[ite,1]
	
x_mode_MHz = savgol_filter(x, 11, 2)
y_mode_amp = savgol_filter(y, 11, 2)


#############PRE SETTINGS########################
#plt.rcParams.update({'fontname':'Times New Roman','font.size': 12})
title_font = {'fontname':'Arial', 'size':'12', 'color':'black', 'weight':'normal',
'verticalalignment':'bottom'} # Bottom vertical alignment for more space
axis_font = {'fontname':'Times New Roman', 'size':'11'}
ticks_font = {'fontname':'Times New Roman', 'size':'11'}
label_font = {'family':'Times New Roman','size':'11'}
legend_font = {'family':'Times New Roman','size':'8'}


fig1 = plt.figure(2,figsize=(15.2/2.53, 5./2.53))
ax1 = fig1.add_subplot(111)
fig1.subplots_adjust(left=0.15)
fig1.subplots_adjust(bottom=0.25)
fig1.subplots_adjust(top=0.95)
fig1.subplots_adjust(right=0.95)


#############PLOT DATA########################
markersize = 10
linewid = 1

x_mode_s = x_mode_MHz * 50. / 52.447
#x_mode_MHz = x_mode_MHz - 1.
x_mode_s = x_mode_s - min(x_mode_s) - 1.

ax1.plot(x_mode_s,y_mode_amp,'k',markersize=markersize,linewidth = linewid)
#plt.axvline(x=0.402, ymin=0.4, ymax = 0.615, linewidth=2, color='k')
#ax1_handle = ax1.scatter(T_arr, wlidl_arr,c=color_branch_arr,vmin=pltoverlapmin, vmax=pltoverlapmax,s=markersize,label ='q'+str(q_arr[qite])+' p'+str(p_arr[pite]), cmap=cm1,marker=currmarker,lw = 0)
#ax1.set_xlim([-50.,50.])
ax1.set_xlim([0.,100.])
ax1.set_ylim([0.7,1.05])




#############LABELS AND LEGENDS########################
#ax1.grid()


for axis in ['top','bottom','left','right']:
	ax1.spines[axis].set_linewidth(1.)


#ticks:
ax1.locator_params(axis = 'x', nbins = 6)
ax1.locator_params(axis = 'y', nbins = 4)


#ax1.set_title('Temp')
#ax1.set_xlabel('Pump laser detuning $\\nu_p -\\nu(\\ell_p,q_p,p_p)$ (GHz)',**axis_font)
ax1.set_xlabel('Time (s)',**axis_font)
ax1.set_ylabel('Reflected power (norm.)',**axis_font)


locs,labels = plt.yticks()
#leg_handle = plt.legend(['spectrum'],loc='lower right',prop=legend_font)


fig1.canvas.draw()


xlabels = [item.get_text() for item in ax1.get_xticklabels()]
ax1.set_xticklabels(xlabels,**ticks_font)
ylabels = [item.get_text() for item in ax1.get_yticklabels()]
ax1.set_yticklabels(ylabels,**ticks_font)


#############FINALIZE########################
fig1.savefig(plots_path + '//' + filename + '.pdf', transparent=True)
fig1.savefig(plots_path + '//' + filename + '.png', dpi=300, transparent=True)


#plt.show()


#pdb.set_trace()
print '\nFarewell, Master!'
plt.close()
sys.exit()









