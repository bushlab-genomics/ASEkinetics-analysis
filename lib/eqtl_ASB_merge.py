from matplotlib.colors import ListedColormap
from adjustText import adjust_text
from scipy.sparse import coo_matrix
from sklearn import linear_model
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib as mpl
import pingouin as pg
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 2
import seaborn as sns
sns.set_style("white")
sns.despine(top=False,right=False,left=False)

from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import pingouin as pg
import dask.dataframe as dd
import pandas as pd
import numpy as np
import re
import sys

def ase_kp2fc(kinetics):

	kinetics[['gene','allele']] = kinetics.index.to_series().str.split('-', n=1, expand=True)
	kinetics['burst size'] = kinetics['ksyn']/kinetics['koff']
	kinetics['burst frequency'] = kinetics['kon'] * kinetics['koff'] / (kinetics['kon'] + kinetics['koff'])
	kinetics_grp = kinetics.groupby('gene')

	record_df = pd.DataFrame(columns = ['gene','mean_fc', 'var_fc', 'bs_fc', 'bf_fc'])
	for igrp, grp in kinetics_grp:
		select = (grp['allele'] == 'maternal_allele')
		grp_m = grp.loc[select]
		grp_p = grp.loc[~select]
		mean_fc = grp_p['mean'].values / grp_m['mean'].values
		mean_fc = np.log2(mean_fc)
		var_fc = grp_p['var'].values / grp_m['var'].values
		var_fc = np.log2(var_fc)

		bs_fc = grp_p['burst size'].values / grp_m['burst size'].values
		bs_fc = np.log2(bs_fc)
		bf_fc = grp_p['burst frequency'].values / grp_m['burst frequency'].values
		bf_fc = np.log2(bf_fc)

		record_df.loc[len(record_df),:]=[igrp, mean_fc[0], var_fc[0], bs_fc[0], bf_fc[0]]

	return record_df

def main():

	# NA12878 input
	asb_input = sys.argv[1]
	eqtl_input = sys.argv[2]
	vcf_input = sys.argv[3]
	kinetics_input = sys.argv[4]

	LD_ref_input = sys.argv[5]
	
	vcf = pd.read_csv(vcf_input, sep='\t', comment="#", header=None, usecols=list(range(5)) + [9], 
		names=['CHROM','POS','ID','REF','ALT','genotype'],dtype={'CHROM':'str'},compression='gzip')

	asb  = pd.read_csv(asb_input, header=0)
	eqtl = pd.read_csv(eqtl_input, header=0, sep='\t', dtype={'chromosome':'str', 'position':'str'} )
	eqtl['variant'] = eqtl['variant'].apply(lambda x: x.replace('_',':'))
	eqtl['variant'] = eqtl['variant'].str.strip('chr')

	eqtl_vcf = pd.merge(eqtl, vcf[['ID','genotype']], left_on='variant', right_on='ID')
	eqtl_vcf.rename(columns={'genotype':'genotype_eqtl','ID':'ID_eqtl'},inplace=True)

	asb.rename(columns={'genotype':'genotype_asb','ID':'ID_asb'},inplace=True)
	df = pd.merge(asb, eqtl_vcf, left_on='gene', right_on='gene_id')

	kinetics = pd.read_csv(kinetics_input, header=0,index_col=0)	
	kpe_fc = ase_kp2fc(kinetics)
	
	LD_ref = pd.read_csv(LD_ref_input, header=None, names=['asb_loci','eqtl_loci','Rsq','Dp'])

	df_kinectics = pd.merge(df, kpe_fc, on='gene')
	df_kp_highLD = pd.merge(df_kinectics, LD_ref, right_on=['asb_loci','eqtl_loci'], \
				left_on=['ID_asb', 'ID_eqtl'])

	output_cols =  ['ID_asb','AR','peak_name','peak.isASB','TF','gene','genotype_asb']
	output_cols += ['rsid','ID_eqtl','pvalue','beta','se','genotype_eqtl','Rsq','Dp'] 
	output_cols += ['mean_fc', 'var_fc', 'bs_fc', 'bf_fc']
	df_kp_highLD[output_cols].to_csv('NA12878_isASB_EP-TSS2k_eqtl_highLD_kp.list',index=False)

def main1():

	df_kp_highLD_input = sys.argv[1]

	# NA18502 input
	vcf_input = sys.argv[2]
	kinetics_input = sys.argv[3]

	vcf = pd.read_csv(vcf_input, sep='\t', comment="#", header=None, usecols=list(range(5)) + [9],\
		names=['CHROM','POS','ID','REF','ALT','genotype'],compression='gzip')
	kinetics = pd.read_csv(kinetics_input, header=0,index_col=0)	
	kpe_fc = ase_kp2fc(kinetics)

	df_kp_highLD = pd.read_csv(df_kp_highLD_input, header=0)
	df_kp_highLD = pd.merge(df_kp_highLD, vcf[['ID','genotype']], left_on='ID_asb', right_on='ID')
	df_kp_highLD.rename(columns={'genotype':'genotype_asb_yri'},inplace=True)


	df_kp_highLD = pd.merge(df_kp_highLD, vcf[['ID','genotype']], left_on='ID_eqtl', right_on='ID')
	df_kp_highLD.rename(columns={'genotype':'genotype_eqtl_yri'},inplace=True)

	df_kp_highLD.to_csv('GM_isASB_EP-TSS2k_eqtl_highLD_kp.list',index=False)

def scatter():

	df_kp_highLD_input = sys.argv[1]
	df_kp_highLD = pd.read_csv(df_kp_highLD_input, header=0)

	df_kp_highLD['beta_pm'] = df_kp_highLD.apply(lambda x: -1 * x.beta if x.genotype_eqtl=='0|1' \
							else x.beta, axis=1)
	record = pd.DataFrame(columns=['TF','beta_v_mean','beta_v_var','beta_v_bs','beta_v_bf'])
	for igrp, grp in df_kp_highLD.groupby('TF'):
		grp_uniq = grp[['beta_pm','mean_fc','var_fc','bs_fc','bf_fc','gene']].drop_duplicates()

		beta_v_mean_pc,pval = stats.spearmanr(grp_uniq['beta_pm'],grp_uniq['mean_fc'])
		beta_v_var_pc,pval = stats.spearmanr(grp_uniq['beta_pm'],grp_uniq['var_fc'])
		beta_v_bs_pc,pval = stats.spearmanr(grp_uniq['beta_pm'],grp_uniq['bs_fc'])
		beta_v_bf_pc,pval = stats.spearmanr(grp_uniq['beta_pm'],grp_uniq['bf_fc'])
	
		record.loc[igrp,:] = [igrp, beta_v_mean_pc, beta_v_var_pc, beta_v_bs_pc, beta_v_bf_pc]

	record.dropna(inplace=True)
	fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)
	ax.scatter(record['beta_v_bs'], record['beta_v_bf'],s=4, c='k')
	
	texts = [ ax.text(row['beta_v_bs'], row['beta_v_bf'],row['TF']) for i, row in record.iterrows() ]
	adjust_text(texts, arrowprops=dict(arrowstyle="->", color='r', lw=1))

	ax.set_xlim(-1,1)
	ax.set_ylim(-1,1)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.get_yaxis().set_tick_params(which='both', direction='out')
	ax.get_xaxis().set_tick_params(which='both', direction='out')
	ax.axvline(x=0,c='k')
	ax.axhline(y=0,c='k')
	ax.set_xlabel('corr(beta, burst size)')
	ax.set_ylabel('corr(beta, burst freq)')
	plt.savefig('beta_v_burst.pdf')
	plt.show()

def plot():

	df_kp_highLD_input = sys.argv[1]
	tf_order_in = sys.argv[2]
	df_kp_highLD = pd.read_csv(df_kp_highLD_input, header=0)
	tf_order     = np.loadtxt(tf_order_in, dtype=str)

	df_kp_highLD = df_kp_highLD.loc[df_kp_highLD['pvalue']<10e-7]

	param = 'Dp'
	df_kp_highLD_rank = df_kp_highLD.groupby('TF')[param].mean()
	df_kp_highLD_rank = df_kp_highLD_rank.sort_values().index
	
	fig,ax = plt.subplots(constrained_layout=True)
	ax = sns.boxplot(y="TF", x=param, data=df_kp_highLD,order=tf_order, width=1,fliersize=0)
	ax.grid(True, axis='y',linestyle='--')
	ax.set(xlim=(-0.1,1.1))

	fig.set_figwidth(5)
	fig.set_figheight(8)
	plt.savefig('%s_distribution.pdf'%(param))
	plt.show()

if __name__ == "__main__":
#	plot()
	scatter()
#	main1()
#	main()


