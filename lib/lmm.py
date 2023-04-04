import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib as mpl
import matplotlib
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2
from scipy import stats
from scipy.sparse import coo_matrix
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
import statsmodels.formula.api as smf
from sklearn import linear_model
from sklearn.cluster import KMeans
import numpy as np
import pandas as pd
import pingouin as pg
import csv
import re
import sys


def pref(x):

	if   x['asb_allele'] == 'H1':
		return  1,0
	elif x['asb_allele'] == 'H2':
		return  0,1

def qref(x):

	if   x['genotype'] == '0|1':
		return  x['Corrected.AR'], 1 - x['Corrected.AR']
	elif x['genotype'] == '1|0':
		return  1 - x['Corrected.AR'], x['Corrected.AR']

def lmm(asb_w_kinetics, variable, kp_list, out):
	
	res_dict = {}
	asb_w_kinetics.reset_index(inplace=True)

	res = pd.DataFrame(columns=['kp','fct','pval','effect size'])
	kp_alias = {'kon':'k+','koff':'k-','ksyn':'r',\
		'burst_size':'Burst size','burst_frequency':'Burst frequency'}
	for kp in kp_list:	
		variable_line = ' + '.join(variable)
		lm = smf.mixedlm(" %s ~ %s "%(kp, variable_line), asb_w_kinetics, groups=asb_w_kinetics['gene'])
		lmf =  lm.fit(method=["powell"])
		lmf.summary().tables[1].to_csv('%s_%s_summary.csv'%(out,kp), quoting=csv.QUOTE_ALL)
		res_dict[kp] = lmf.summary().tables[1][:-1].astype(float)
		print(res_dict[kp])

	return res_dict

def lmm_vis(res_dict, out):

	res = pd.DataFrame(columns=['kp','fct','pval','effect size'])
	kp_alias = {'kon':'k+','koff':'k-','ksyn':'r',\
		'burst_size':'Burst size','burst_frequency':'Burst frequency'}

	kp_list = ['kon','koff','ksyn','burst_size','burst_frequency']
	
	sig_tf_all = []
	for i, kp in enumerate(kp_list):
		res_kp = res_dict[kp]
		res_kp['P>|z|'] = res_kp['P>|z|'].astype(float)
		sig_tf = res_kp.loc[(res_kp['P>|z|'] < 0.05/(len(kp_list)))].index
		sig_tf_all = set(sig_tf_all).union(set(sig_tf))

		p_vals = res_kp['P>|z|'].astype(float)
		rejected, pvalue_corrected = pg.multicomp(p_vals,alpha=0.05,method='fdr_bh')
		res_dict[kp]['pval-corrected'] = pvalue_corrected

	sig_tf_all.remove('Intercept')
	sig_tf_all = res_dict[kp].loc[sig_tf_all].sort_values(by='pval-corrected').index

	fig, ax = plt.subplots(ncols=len(res_dict),figsize=(8,3))
	for i, kp in enumerate(kp_list):
		res_kp = res_dict[kp]
		for irow, itf in enumerate(sig_tf_all):
			x_coord = (res_kp.loc[itf,'[0.025'],res_kp.loc[itf,'0.975]'])
			color = 'tab:red' if res_kp.loc[itf, 'pval-corrected'] < 0.05 else(\
				'tab:orange' if res_kp.loc[itf, 'P>|z|'] < 0.05  else 'k' )

			ax[i].scatter(res_kp.loc[itf, 'Coef.'], irow, marker='s',c=color)	
			ax[i].plot(x_coord, (irow, irow), c=color)
			if i > 0:  ax[i].set_yticks([])
			ax[i].axvline(x=0,color='k',linewidth=1)
			ax[i].axhline(y=irow,color='gray',linewidth=1, ls='--',alpha=0.4)

		max_val = np.max(abs(res_kp.loc[sig_tf_all, ['[0.025','0.975]']].values))
		ax[i].set_xlim([-1*max_val,max_val])
		ax[i].set_title(kp_alias[kp])
	ax[0].set_yticks(np.arange(len(sig_tf_all)))
	ax[int(np.floor(len(res_dict)/2))].set_xlabel('effect size')
	ax[0].set_yticklabels(sig_tf_all, style='italic')
	fig.subplots_adjust(wspace=0.1, hspace=0.2)
	plt.tight_layout()
	plt.xticks()
	plt.savefig(out+'.pdf',dpi=300,bbox_inches ='tight')
	plt.savefig(out+'.png',dpi=300,bbox_inches ='tight')
	plt.show()

def main():

	asb_input = sys.argv[1]
	vcf_input = sys.argv[2]
	kinetics_input = sys.argv[3]
	output = sys.argv[4]
	
	asb = pd.read_csv(asb_input, header=0)
	asb = asb.loc[asb['peak.isASB']==True]
	vcf = pd.read_csv(vcf_input, header=None, sep='\t',usecols=list(range(6)) + [9,11,12,13],\
		names=['chr','pos','ID','ref','alt','qual','genotype','reg_starts','reg_ends','gene'])

	asb_w_genotype = pd.merge(asb,vcf,on='ID', suffixes=('_asb',''))
	asb_w_genotype[['H1_occu','H2_occu']] = asb_w_genotype.apply(qref, axis=1, result_type="expand")
	asb_w_genotype_reform=pd.melt(asb_w_genotype,id_vars=['peak_name','TF','gene'],\
					value_vars=['H1_occu','H2_occu'],var_name='allele',value_name='occu')
	asb_w_genotype_reform['allele']  = asb_w_genotype_reform['allele'].str.replace('_occu','_allele')
	asb_w_genotype_reform.drop_duplicates(inplace=True)

	kinetics = pd.read_csv(kinetics_input, header=0,index_col=0)
	kinetics.dropna(inplace=True)
	kinetics[['gene','allele']] = kinetics.index.to_series().str.split('-', n=1, expand=True)
	kinetics['burst_size'] = kinetics['ksyn']/kinetics['koff']
	kinetics['burst_frequency'] = kinetics['kon'] * kinetics['koff'] / (kinetics['kon'] + kinetics['koff'])

	kp_list = ['kon','koff','ksyn','burst_size','burst_frequency']
	variable = asb_w_genotype_reform['TF'].unique()

	asb_w_genotype_reform['allele'] = asb_w_genotype_reform[['gene','allele']].agg('-'.join, axis=1)
	mat = asb_w_genotype_reform[['allele','TF','occu']].groupby(['allele','TF'])['occu'].sum()
	mat = mat.to_frame(name='sum_occu').reset_index()

	mat['occu'] = mat['sum_occu']
	mat = mat.pivot(index='allele', columns='TF',values='occu').fillna(0)
	mat_kinetics = pd.merge(mat, kinetics,left_index=True, right_index=True)
	variable_nozero = variable[mat_kinetics[variable].sum(axis=0)>1]
	print(len(variable_nozero))
	print(mat_kinetics['gene'].nunique())
	res_dict = lmm(mat_kinetics, variable_nozero, kp_list, output)

	lmm_vis(res_dict,output)	
	
if __name__ == "__main__":
	main()


