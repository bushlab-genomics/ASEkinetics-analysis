import matplotlib.ticker as ticker
from matplotlib.colors import ListedColormap
from matplotlib.ticker import AutoMinorLocator
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib as mpl
import matplotlib
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 16
plt.rcParams['axes.linewidth'] = 2

from scipy import stats
import numpy as np
import pandas as pd
import re
import sys

def pear_corr(df, output):

	fig,ax = plt.subplots(ncols=3,nrows=3,sharex=True,sharey='row',\
			constrained_layout=True,figsize=(8, 8))
	xvar = 'PvsM_es'
	fig_dict = {
		'Group1':['log(mean_fc)', 'log(var_fc)','NA'], \
		'Group2':['log(kon_fc)', 'log(burst frequency_fc)','NA'], \
		'Group3':['log(koff_fc)', 'log(ksyn_fc)','log(burst size_fc)']}

	fig_dict_alias = {
		'log(mean_fc)':'log2($mean_{H2}$/$mean_{H1}$)', \
		'log(var_fc)':'log2($var_{H2}$/$var_{H1}$)', \
		'log(kon_fc)':'log2($k_{+H2}$/$k_{+H1}$)', \
		'log(burst frequency_fc)':'log2($bs_{H2}$/$bs_{H1}$)', \
		'log(koff_fc)':'log2($k_{-H2}$/$k_{-H1}$)', \
		'log(ksyn_fc)': 'log2($r_{H2}$/$r_{H1}$)', \
		'log(burst size_fc)': 'log2($bs_{H2}$/$bs_{H1}$)'}

	fig_dict_title = {
		'log(mean_fc)':'Mean', \
		'log(var_fc)':'Variance', \
		'log(kon_fc)':'k+', \
		'log(burst frequency_fc)':'Burst frequency', \
		'log(koff_fc)':'k-', \
		'log(ksyn_fc)': 'r', \
		'log(burst size_fc)': 'Burst size'}

	for row, (key, val_list) in enumerate(fig_dict.items()):
		for col, val in enumerate(val_list):	
			
			if val == 'NA':
				ax[row,col].axis('off')
				continue
			df_va1 = df[xvar].values
			df_va2 = df[val].values

			p_stat,   p_pval = stats.pearsonr(df_va1, df_va2)
			sp_stat, sp_pval = stats.spearmanr(df_va1, df_va2)

			if key in [ 'Group1', 'Group2']: color_code= '#d62728'
			else: color_code = 'k'
			ax[row,col].scatter(df_va1, df_va2, c=color_code,s=2)	
			ax[row,col].xaxis.set_minor_locator(AutoMinorLocator())
			ax[row,col].set_ylabel(fig_dict_alias[val])	

			if sp_pval < 0.01: sp_pval_str = str('{0:.2e}'.format(sp_pval))
			else:sp_pval_str = str('{0:.3f}'.format(sp_pval))
			sp_stat_str = str('{0:.3f}'.format(sp_stat))

			ax[row,col].text(0.95, 0.95, \
				'correlation coefficient=%s\np-value=%s'%(sp_stat_str, sp_pval_str),\
				horizontalalignment='right',verticalalignment='top',\
				transform = ax[row,col].transAxes)

			ax[row,col].set_title(fig_dict_title[val])
			x = np.arange(-2, 2, 0.1)
			y = p_stat * x
			ax[row,col].set_xlim(-1.5,1.5)
			ax[row,col].set_ylim(-5,5)
			ax[row,col].plot(x,y,ls='-', lw=2, c='k')
			ax[row,col].axhline(y=0, ls='--',lw=1, c='k')
			ax[2,col].set_xlabel('eQTL effect size(H2/H1)')
	plt.savefig(output + '.pdf')
	plt.savefig(output + '.png')
	plt.show()

def main():

	eqtl_input = sys.argv[1]
	vcf_input  = sys.argv[2]
	kinetics_input = sys.argv[3]
	output = sys.argv[4]
	
	eqtl = pd.read_csv(eqtl_input, header=0, sep='\t')
	vcf = pd.read_csv(vcf_input, header=None, sep='\t', comment='#' ,usecols=list(range(6)) + [9],\
			names=['chr','pos','ID','ref','alt','qual','genotype'],compression='gzip', dtype={'chr':str})
	kinetics = pd.read_csv(kinetics_input, header=0,index_col=0)
	kinetics[['gene','allele']] = kinetics.index.to_series().str.split('-', n=1, expand=True)
	kinetics['burst size'] = kinetics['ksyn']/kinetics['koff']
	kinetics['burst frequency'] = kinetics['kon'] * kinetics['koff'] / (kinetics['kon'] + kinetics['koff'])
	kinetics['fraction'] = kinetics['kon'] / (kinetics['kon'] + kinetics['koff'])

	eqtl['ID'] = eqtl[['chromosome','position','ref','alt']].astype(str).agg(':'.join, axis=1)
	
	vcf['ID'] = vcf[['chr','pos','ref','alt']].astype(str).agg(':'.join, axis=1)
	eqtl_vcf = pd.merge(eqtl, vcf[['ID','genotype']], on='ID')
	eqtl_vcf['PvsM_es'] = eqtl_vcf.apply(lambda x: x['beta'] if x['genotype']=='0|1' else \
			 (-x['beta'] if x['genotype']=='1|0' else np.nan), axis = 1)

	H1_kinetics = kinetics.loc[kinetics['allele']=='H1_allele']
	H2_kinetics = kinetics.loc[kinetics['allele']=='H2_allele']
	pair_kinetics = pd.merge(H1_kinetics, H2_kinetics, on='gene', suffixes=('_m','_p'))
	fc_cols=[]
	for kp in ['kon','koff','ksyn','burst size', 'burst frequency','mean','var']:
		pair_kinetics['%s_fc'%kp] = pair_kinetics['%s_p'%kp] / pair_kinetics['%s_m'%kp]
		pair_kinetics['log(%s_fc)'%kp] = pair_kinetics['%s_fc'%kp].transform(np.log10)	
		fc_cols +=  [ '%s_fc'%kp, 'log(%s_fc)'%kp ]
	fc_cols += ['gene']
	pair_kinetics_genotype = pd.merge(pair_kinetics[fc_cols], \
			eqtl_vcf[['molecular_trait_id','PvsM_es','ID','rsid','pvalue','genotype']], \
			left_on='gene', right_on='molecular_trait_id')
	
	pair_kinetics_genotype.to_csv(output, index=False)

def main1():

	df_in = sys.argv[1]
	output = sys.argv[2]
	pair_kinetics_genotype = pd.read_csv(df_in, header=0 )
	pair_kinetics_filter = pair_kinetics_genotype.sort_values('pvalue').groupby('gene').head(3)
	
	pear_corr(pair_kinetics_filter, output)	

if __name__ == "__main__":
#	main()
	main1()
