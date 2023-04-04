import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib as mpl
import matplotlib
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2

import seaborn as sns
sns.set_style("white")
sns.despine(top=False,right=False,left=False)

from VCT import *
from scipy import stats
from scipy.sparse import coo_matrix
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
import statsmodels.formula.api as smf
import umap
from sklearn import linear_model
from sklearn.cluster import KMeans
import numpy as np
import pandas as pd
import pingouin as pg
import re
import sys

def qref(x):
	
	if   x['genotype'] == '0|1':
		return  x['Corrected.AR'], 1 - x['Corrected.AR']
	elif x['genotype'] == '1|0':
		return  1 - x['Corrected.AR'], x['Corrected.AR']

def pref(x):

	if   x['asb_allele'] == 'H1':
		return  1,0
	elif x['asb_allele'] == 'H2':
		return  0,1

def single_tf_di(tf, asb_w_kinetics_all, output):
	
	asb_w_kinetics = asb_w_kinetics_all[asb_w_kinetics_all['TF']==tf]
	kp=['kon','koff','ksyn','burst_size','burst_frequency']
	kp_alias = {'kon':'k+','koff':'k-','ksyn':'r',\
		'burst_size':'Burst size','burst_frequency':'Burst frequency'}

	asb_w_kinetics.to_csv(output+'_'+tf+'.csv',index=False)
	nfig=len(kp)
	fig,ax = plt.subplots(ncols=nfig, sharex=True,constrained_layout=True,figsize=(4*nfig, 3))
	for i, ikp in enumerate(kp):
		x  = [0,1]
		bool_w_occu = (asb_w_kinetics['occu'] > 0)
		subdf = asb_w_kinetics[['peak_name','TF','occu',ikp,'allele','gene']].drop_duplicates()
		subdf_w_occu = subdf.loc[bool_w_occu]
		subdf_wo_occu = subdf.loc[~bool_w_occu]
		
		subdf_reform = pd.merge(subdf_w_occu, subdf_wo_occu, on=['gene','peak_name'],suffixes=('','_0'))
		subdf_reform.to_csv('%s_%s.csv'%(output,tf), index=False)

		rvs1 = subdf_reform[ikp].values
		rvs0 = subdf_reform[ikp+'_0'].values

		if len(rvs1) < 6 or len(rvs0) < 6: continue	
		wx_res = pg.wilcoxon(x=rvs1,y=rvs0)

		pval = wx_res.loc['Wilcoxon','p-val']
		pval = round(pval,4)

		es = wx_res.loc['Wilcoxon','CLES'] - 0.5
		es = round(es,4)

		ax[i].scatter(np.zeros(len(rvs0)),rvs0,s=4, c='k')
		ax[i].scatter(np.ones(len(rvs1)), rvs1,s=4, c='k')
		for ii in range(len(rvs0)):
			x = [0,1]
			y = [rvs0[ii], rvs1[ii]]
			ax[i].plot(x, y, color="black", lw=1, alpha=0.5)

		ax[i].set_title(kp_alias[ikp])
		ax[i].set_xticks([0,1])
		ax[i].set_xlim(-1,2)
		ax[i].set_ylim(min(np.hstack([rvs0,rvs1]))-0.1, max(np.hstack([rvs1,rvs0]))*1.3)
		ax[i].text(0.95, 0.95,'U=%s\npval=%s'%(es, pval), \
			horizontalalignment='right',verticalalignment='top',transform = ax[i].transAxes)
	ax[2].set_xlabel('TF binding status')
	fig.set_figwidth(10)
	plt.savefig(output + '_' + tf + '.pdf')

def lmm(asb_w_kinetics, variable, kp_list, out):

	asb_w_kinetics.reset_index(inplace=True)

	trans = umap.UMAP().fit(asb_w_kinetics[variable])
	kmeans = KMeans(n_clusters=2, random_state=0).fit(trans.embedding_)
	asb_w_kinetics['group'] = kmeans.labels_

	res = pd.DataFrame(columns=['kp','fct','pval','effect size'])
	kp_alias = {'kon':'k+','koff':'k-','ksyn':'r',\
		'burst_size':'Burst size','burst_frequency':'Burst frequency'}

	variable_nozero = variable[asb_w_kinetics[variable].sum(axis=0) > 0]
	asb_w_kinetics.to_csv('temp.pair')

	for kp in kp_list:	
		variable_line = ' + '.join(variable_nozero)
		lm = smf.mixedlm(" %s ~ %s "%(kp, variable_line), asb_w_kinetics, groups=asb_w_kinetics['group'])
		lmf =  lm.fit(method=["lbfgs"])
		lmf.summary().tables[1].to_csv('%s_summary.csv'%kp)

def MWU_fct_analysis(asb_w_kinetics, variable, kp_list, out):

	res = pd.DataFrame(columns=['kp','fct','pval','effect size'])
	kp_alias = {'kon':'k+','koff':'k-','ksyn':'r',\
		'burst_size':'Burst size','burst_frequency':'Burst frequency'}

	i = 0
	for kp in kp_list:
		for fct in variable:
			bool_w_fct = (asb_w_kinetics['TF'] == fct)
			subdf = asb_w_kinetics.loc[bool_w_fct,['peak_name','TF','occu',kp,'allele','gene']].drop_duplicates()
			bool_w_occu = (subdf['occu'] > 0)
			subdf_w_occu = subdf.loc[bool_w_occu]
			subdf_wo_occu = subdf.loc[~bool_w_occu]

			subdf_reform = pd.merge(subdf_w_occu, subdf_wo_occu, on=['gene','peak_name'],suffixes=('','_0'))
			if fct == 'YY1': subdf_reform.to_csv('YY1')
			rvs1 = subdf_reform[kp]
			rvs0 = subdf_reform[kp+'_0']

			if len(rvs1) < 10 or len(rvs0) < 10: continue
			wx_res = pg.wilcoxon(x=rvs1,y=rvs0)

			res.loc[i, 'kp' ] = kp_alias[kp]
			res.loc[i, 'fct'] = fct
			res.loc[i, 'pval'] = wx_res.loc['Wilcoxon','p-val']
			res.loc[i, 'effect size'] = wx_res.loc['Wilcoxon','CLES']
			res.loc[i, 'effect size'] -= 0.5
			i += 1

	for ikp, ikp_alias in kp_alias.items():
		kp_bool = (res['kp']==ikp_alias)
		p_vals = res.loc[kp_bool,'pval'].values.astype(float)
		rejected, pvalue_corrected = pg.multicomp(p_vals,alpha=0.1,method='fdr_bh')
		res.loc[kp_bool,'pval-corrected'] = pvalue_corrected

	if res['fct'].nunique() > 1:
		res_filter_fct = res.loc[res['pval-corrected']<0.1,'fct'].unique()
	else: res_filter_fct = res['fct'].unique()

	order = res.sort_values(by='pval-corrected')['fct'].unique()
	g = sns.FacetGrid(res, col="kp",despine=False)
	g.map_dataframe(sns.barplot,x='effect size',y='fct',order=order,color='k', orient='h', dodge=False, errwidth=0)
	g.set_xlabels('')

	for ax in g.axes.flat:
		ax.grid(True, axis='x',linestyle='--')
		ax.grid(True, axis='y',linestyle='--')

	g.set(xlim=(-0.4, 0.4))
	g.map(plt.axvline,x=0,color='k',linewidth=1)
	g.set_xticklabels(fontdict={'fontstyle':'italic'})

	g.set_titles('{col_name}')
	g.fig.text(0.4, 0.01, 'Effect size')
	g.fig.set_figwidth(10)
	g.fig.set_figheight(np.ceil(len(res)/24))
	plt.tight_layout()
	plt.savefig(out+'_mwu.pdf',dpi=300)

	return res, res_filter_fct

def lr_analysis(asb_w_kinetics,variable,kp_list,out):

	df_h2 = pd.DataFrame(columns=['kp','mean','var','std'])
	for i, ikp in enumerate(kp_list):
		y = asb_w_kinetics[ikp]
		#y = asb_w_kinetics[ikp].astype(float).values.reshape((-1,1))
		X = asb_w_kinetics[variable].astype(float)
		K = np.dot(X, X.T)/len(variable)
		dummy_cov = np.ones(len(X)).reshape((-1,1))
		#lmm = lmm_cov.LMM(Y=y, X=dummy_cov, G=None, K=K, regressX=True, inplace=False, forcefullrank=False)
		h2, h2_mean, h2_var, grid, post_h2 =est_h2(y.values, K, covariates=None, nGridH2=10000, plot=True, verbose=True)
		df_h2.loc[i,:] = [ikp, h2_mean, h2_var, np.sqrt(h2_var)]
	df_h2.to_csv(out+'_h2.csv',index=False)

def allele_ratio(asb_w_genotype_reform,out):
	
	ar = asb_w_genotype_reform.groupby(['TF','allele'])['occu'].sum()
	ar_percent = ar / asb_w_genotype_reform.groupby(['TF','allele'])['occu'].count()

	ar_df = ar.to_frame().reset_index()
	ar_df_pivot = ar_df.pivot(index='TF', columns='allele', values='occu')
	ar_df_pivot.reset_index(inplace=True)

	ar_df_pivot.sort_values(by='H1_allele', inplace=True)	
	labels=ar_df_pivot['TF']
	m_ar = ar_df_pivot['H1_allele']
	p_ar = ar_df_pivot['H2_allele']

	ar_percent_df = ar_percent.to_frame().reset_index()
	ar_percent_df_pivot = ar_percent_df.pivot(index='TF', columns='allele', values='occu')
	ar_percent_df_pivot.reset_index(inplace=True)

	ar_percent_df_pivot.set_index('TF',inplace=True)	
	m_ar_percent = ar_percent_df_pivot.loc[labels,'H1_allele']
	p_ar_percent = ar_percent_df_pivot.loc[labels,'H2_allele']
		
	fig, ax = plt.subplots(figsize=(8,10),constrained_layout=True, ncols=2)
	ax[0].barh(labels, m_ar, label='H1')
	ax[0].barh(labels, p_ar, left=m_ar,label='H2')

	ax[1].barh(labels, m_ar_percent, label='H1')
	ax[1].barh(labels, p_ar_percent, left=m_ar_percent,label='H2')
	ax[1].legend(frameon=False, bbox_to_anchor=[1.01,0.5])

	plt.savefig(out+'_ASB_dist.pdf',dpi=300,bbox_inches ='tight')
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
	asb_w_genotype[['H1_occu','H2_occu']] = asb_w_genotype.apply(pref, axis=1, result_type="expand")
#	asb_w_genotype[['H1_occu','H2_occu']] = asb_w_genotype.apply(qref, axis=1, result_type="expand")

	asb_w_genotype_reform=pd.melt(asb_w_genotype,id_vars=['peak_name','TF','gene'],\
		value_vars=['H1_occu','H2_occu'],var_name='allele',value_name='occu')
	asb_w_genotype_reform['allele']  = asb_w_genotype_reform['allele'].str.replace('_occu','_allele')

	asb_w_genotype_reform.drop_duplicates(inplace=True)

	kinetics = pd.read_csv(kinetics_input, header=0,index_col=0)
	kinetics.dropna(inplace=True)
	kinetics[['gene','allele']] = kinetics.index.to_series().str.split('-', n=1, expand=True)
	kinetics['burst_size'] = kinetics['ksyn']/kinetics['koff']
	kinetics['burst_frequency'] = kinetics['kon'] * kinetics['koff'] / (kinetics['kon'] + kinetics['koff'])
	asb_w_genotype_kp = pd.merge(asb_w_genotype_reform, kinetics,on=['gene','allele'])

	kp_list = ['kon','koff','ksyn','burst_size','burst_frequency']
	variable = asb_w_genotype_kp['TF'].unique()

#	mwu_res, filtered_fct = MWU_fct_analysis(asb_w_genotype_kp,variable,kp_list,output)
#	mwu_res.sort_values(by='pval-corrected').to_csv(output+'_mwu_sig_gene.csv', index=False)
#	for ifct in filtered_fct:
#		single_tf_di(ifct, asb_w_genotype_kp, output)
	
	asb_w_genotype_reform['allele'] = asb_w_genotype_reform[['gene','allele']].agg('-'.join, axis=1)
	mat = asb_w_genotype_reform[['allele','TF','occu']].groupby(['allele','TF'])['occu'].sum()
	mat = mat.to_frame(name='sum_occu').reset_index()
	mat['occu'] = mat['sum_occu']
	mat['occu'] = mat['sum_occu'].apply(lambda x: 1 if x > 0 else 0)

	mat = mat.pivot(index='allele', columns='TF',values='occu').fillna(0)
	mat_kinetics = pd.merge(mat, kinetics,left_index=True, right_index=True)
	lr_analysis(mat_kinetics,variable,kp_list,output)
	mat_kinetics.to_csv(output+'.csv')
#	lmm(mat_kinetics, variable, kp_list, output)

if __name__ == "__main__":
	main()



