import matplotlib.ticker as ticker
from matplotlib.colors import ListedColormap
from matplotlib.ticker import AutoMinorLocator
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib as mpl
import matplotlib
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 2

from scipy import stats
import numpy as np
import pandas as pd
import re
import sys


def plot_sgl(df_combine, labels):

	ncols=2
	fig,ax = plt.subplots(figsize=(8, 8),ncols=ncols,constrained_layout=True)

	for i in range(ncols):
		df = df_combine.loc[df_combine['select']==0]
		df1 = df_combine.loc[df_combine['select']==1]
	
		va2,va1 = [ 'log2(burst frequency_fc)', 'log2(burst size_fc)']
		va2_alias,va1_alias = [ 'log2($bf_{H2}$/$bf_{H1}$))', 'log2($bs_{H2}$/$bs_{H1}$))']
	
		colors = ListedColormap(['gray','#e97817'])
		labels=['Other genes',\
			'Genes w/ log2($mean_{H2}$/$mean_{H1}$) < log2($var_{H2}$/$var_{H1}$)']
		ax[i].scatter(df[va1], df[va2], c='gray',    s=4)
		ax[i].scatter(df1[va1],df1[va2],c='#e97817', s=10, marker='D',linewidths=0.5, edgecolors='k')
	
		ax[i].xaxis.set_minor_locator(AutoMinorLocator())
		ax[i].set_xlabel(va1_alias)
		ax[i].set_ylabel(va2_alias)

		axislim = -6*i +8
		x=np.arange(-axislim,axislim, 0.1)
		y=np.arange(-axislim,axislim, 0.1)
		ax[i].set_xlim(-axislim,axislim)
		ax[i].set_ylim(-axislim,axislim)
		ax[i].plot(x,-y,ls='dotted', c='gray',lw=1)
		ax[i].axhline(y=1, ls='dotted', c='gray',lw=1)
		ax[i].axhline(y=-1, ls='dotted', c='gray',lw=1)
		ax[i].axhline(y=0, ls='dotted', c='gray',lw=1)
		ax[i].axvline(x=0, ls='dotted', c='gray',lw=1)
		ax[i].axvline(x=1, ls='dotted', c='gray',lw=1)
		ax[i].axvline(x=-1, ls='dotted', c='gray',lw=1)

		if i==0:
			ax[i].vlines(-1, -1, 1, lw=2, color='k')
			ax[i].vlines(1, -1, 1,  lw=2, color='k')
			ax[i].hlines(-1, -1, 1, lw=2, color='k')
			ax[i].hlines(1, -1, 1,  lw=2, color='k')
		ax[i].set_aspect('equal', adjustable='box')	

	legend_elements = [Line2D([0], [0], marker='o', color='w', label=labels[0],
			markerfacecolor='gray', markersize=8)]
	legend_elements += [Line2D([0], [0], marker='D', color='w', label=labels[1],
			markerfacecolor='#e97817', markersize=10)]

	plt.legend(handles=legend_elements, bbox_to_anchor=(0.1, 1.05), frameon=False)
	plt.savefig('kp_distribution.pdf')
	plt.savefig('kp_distribution.png')
	plt.show()

def plot2(df):

	fig,ax = plt.subplots(ncols=2,sharex=True,sharey=True,constrained_layout=True,figsize=(6, 3))

	d = {'Group1':['log2(mean_fc)', 'log2(burst frequency_fc)'], 'Group2':['log2(var_fc)', 'log2(burst frequency_fc)']}
	for i, (key, value) in enumerate(d.items()):

		va2,va1 = value
		df_va1 = df[va1]
		df_va2 = df[va2]

		norm = mpl.colors.Normalize(vmin=0, vmax=10)
		if i > 2: p1 = ax[i].scatter(df_va1, df_va2, c=abs(df['log2(burst frequency_fc)']), \
					cmap='winter_r', s=6, norm=norm)
		else: ax[i].scatter(df_va1, df_va2, c='tab:blue',s=6)
		ax[i].xaxis.set_minor_locator(AutoMinorLocator())
		ax[i].set_xlim(-3,3)
		ax[i].set_ylim(-3,3)
		ax[i].set_xlabel(va1)
		ax[i].set_ylabel(va2)

		x=np.arange(-3, 3, 0.1)
		y=np.arange(-3, 3, 0.1)
		ax[i].plot(x,y,ls='dotted', c='gray',lw=1)
		ax[i].axhline(y=1, ls='dotted', c='gray',lw=1)
		ax[i].axhline(y=-1, ls='dotted', c='gray',lw=1)
		ax[i].axhline(y=0, ls='dotted', c='gray',lw=1)
		ax[i].axvline(x=0, ls='dotted', c='gray',lw=1)
		ax[i].axvline(x=1, ls='dotted', c='gray',lw=1)
		ax[i].axvline(x=-1, ls='dotted', c='gray',lw=1)
		ax[i].plot(x,y,ls='dotted', c='gray', lw=1)

#	cbar = plt.colorbar(p1, ax=ax[1])
#	cbar.ax.set_ylabel('log2(burst frequency_fc)', rotation=270)
	plt.savefig('kp_distribution.png')
	plt.show()

def main():

	kinetics_input_list = sys.argv[1].split(',')
	labels = sys.argv[2].split(',')
	output = sys.argv[3]

	df_combine = pd.DataFrame()	
	for i, kinetics_input in enumerate(kinetics_input_list):

		kinetics = pd.read_csv(kinetics_input, header=0,index_col=0)
		kinetics[['gene','allele']] = kinetics.index.to_series().str.split('-', n=1, expand=True)
		kinetics['burst size'] = kinetics['ksyn']/kinetics['koff']
		kinetics['burst frequency'] = kinetics['kon'] * kinetics['koff'] / (kinetics['kon'] + kinetics['koff'])
		kinetics['fraction'] = kinetics['kon'] / (kinetics['kon'] + kinetics['koff'])
	
		H1_kinetics = kinetics.loc[kinetics['allele']=='H1_allele']
		H2_kinetics = kinetics.loc[kinetics['allele']=='H2_allele']
		df = pd.merge(H1_kinetics, H2_kinetics, on='gene', suffixes=('_m','_p'))
		fc_cols=[]
		for kp in ['kon','koff','ksyn','burst size', 'burst frequency','mean','var']:
			df['%s_fc'%kp] = df['%s_p'%kp] / df['%s_m'%kp]
			df['log2(%s_fc)'%kp] = df['%s_fc'%kp].transform(np.log2)	
	

		#select0 = (abs(df['log2(burst frequency_fc)']) > np.log2(2)) | (abs(df['log2(burst size_fc)']) > np.log2(2)) 
		#select = select0
		#df_select = df.loc[select]
		select0 = True
		df_select  = df	

		select00 = ( 2*abs(df['log2(mean_fc)'])) < abs(df['log2(var_fc)'])
		select1 = select0 & select00	
	
		df_select['select'] = 0
		df_select.loc[select1,'select'] = 1
		df_select['sample'] = labels[i]
		if len(df_combine) < 1: df_combine = df_select
		else: df_combine = pd.concat([df_combine, df_select])

	output_cols = ['sample','gene','log2(burst size_fc)', 'log2(burst frequency_fc)', 'log2(mean_fc)', 'log2(var_fc)','select']
	df_combine[output_cols].drop_duplicates().to_csv('%s_burstfc.csv'%output,index=False)
	plot_sgl(df_combine,labels)

if __name__ == "__main__":
	main()









