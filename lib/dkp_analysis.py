import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 2

from scipy import stats
import pingouin as pg
import pandas as pd
import seaborn as sns
import numpy as np
import pickle
import sys
import re

def bar_vis_allele(kpe_gene, gene, output):


	kp_list = ['burst_size', 'burst_frequency']
	kp_label = {'burst_size':'Burst size', 'burst_frequency':'Burst frequency'}

	fig,ax = plt.subplots(constrained_layout=True,figsize=(len(kp_list)*1.8,2),\
		ncols=len(kp_list),sharex=True)

	palette ={"H1_allele": "#ff7a7a", "H2_allele": "#2876df"}
	for i, ikp in enumerate(kp_list):
		sns.barplot(data=kpe_gene, x='sample',y=ikp, hue='allele',ax=ax[i], palette=palette)

		x_coords = [p.get_x() + 0.5 * p.get_width() for p in ax[i].patches]
		y_coords = [p.get_height() for p in ax[i].patches]

		upper_error = kpe_gene['%s_upper'%ikp].tolist()
		lower_error = kpe_gene['%s_low'%ikp].tolist()
		asymmetric_error = [lower_error, upper_error]

		ax[i].errorbar(x=x_coords, y=y_coords, yerr=asymmetric_error, fmt="none", c="k", linewidth=1, capsize=2)
		ax[i].set_title(kp_label[ikp],loc='right')
		ax[i].set_xlabel('')
		ax[i].legend()
		ax[i].set_ylim([0,1.6*max(kpe_gene[ikp].values)])
	ax[0].set_ylabel('Values(a.u.)')
	ax[0].get_legend().remove()
	ax[1].set_ylabel('')
	handles, labels = ax[1].get_legend_handles_labels()
	sns.move_legend(ax[1], "lower center", bbox_to_anchor=(.5, 1.1), ncol=3, title=None, frameon=False)
	plt.tight_layout()
	plt.savefig('%s_bar_kp.pdf'%(output), dpi=300)
	plt.show()

def main():

	kpe_list  = sys.argv[1].split(',')
	sample_list   = sys.argv[2].split(',')
	gene       = sys.argv[3]
	output     = sys.argv[4]

	kpe_gen = []
	kpe_all = pd.DataFrame()
	for i, kpe_temp_in in enumerate(kpe_list):
		kpe_temp = pd.read_csv(kpe_temp_in, header=0, index_col=0).dropna()
		kpe_temp[['gene','allele']] = kpe_temp.index.to_series().str.split('-',expand=True)
		kpe_temp['sample'] = sample_list[i]
 
		kpe_var_in = kpe_temp_in.strip('_pb_qc.est') + '_pb.var'
		kpe_var = pd.read_csv(kpe_var_in, header=0, index_col=0)
		print(kpe_temp)
		print(kpe_var)
		kpe_temp_var = pd.merge(kpe_temp[['gene','allele','sample']], kpe_var, left_index=True, right_index=True, suffixes=('','_var'))
	
		kpe_temp_var['burst_size'] = kpe_temp_var['ksyn']/kpe_temp_var['koff']
		kpe_temp_var['burst_frequency'] = kpe_temp_var['kon'] * kpe_temp_var['koff'] /(kpe_temp_var['kon'] + kpe_temp_var['koff'])

		kpe_temp_var['burst_size_upper'] = kpe_temp_var['ksyn_upper']/kpe_temp_var['koff_low'] - kpe_temp_var['burst_size']
		kpe_temp_var['burst_size_low'] = kpe_temp_var['burst_size'] - kpe_temp_var['ksyn_low']/kpe_temp_var['koff_upper']

		kpe_temp_var['burst_frequency_upper'] =  kpe_temp_var['kon_upper'] * kpe_temp_var['koff_upper']/(kpe_temp_var['kon_low'] + kpe_temp_var['koff_low']) - kpe_temp_var['burst_frequency']
		kpe_temp_var['burst_frequency_low'] = kpe_temp_var['burst_frequency'] - kpe_temp_var['kon_low'] * kpe_temp_var['koff_low']/(kpe_temp_var['kon_upper'] + kpe_temp_var['koff_upper']) 

		kpe_gene = kpe_temp_var.loc[kpe_temp_var['gene']==gene]
		print(kpe_gene)
		kpe_all = kpe_all.append(kpe_gene,ignore_index=False)
	
	bar_vis_allele(kpe_all,gene, output)

if __name__ == "__main__":
	main()




