import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib as mpl
import matplotlib
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2

from scipy import stats
import umap 
from sklearn import linear_model
from sklearn.cluster import KMeans

import numpy as np
import pandas as pd
import re
import sys


def pref(x):

	if   x['asb_allele'] == 'H1':
		return  1,0
	elif x['asb_allele'] == 'H2':
		return  0,1

def lmm(asb_w_kinetics, variable, kp_list, out):

	asb_w_kinetics.reset_index(inplace=True)
	trans = umap.UMAP().fit(asb_w_kinetics[variable])
	kmeans = KMeans(n_clusters=3, random_state=0).fit(trans.embedding_)

#	fig, ax = plt.subplots(ncols=3)	
#	ax[0].scatter(trans.embedding_[:, 0], trans.embedding_[:, 1], s=3,c=kmeans.labels_,cmap='Spectral')	
#	ax[1].scatter(asb_w_kinetics['log(burst_size_fc)'],asb_w_kinetics['log(burst_frequency_fc)'], s=3,c=kmeans.labels_,cmap='Spectral')	
#	ax[2].scatter(asb_w_kinetics['log(kon_fc)'],asb_w_kinetics['log(koff_fc)'], s=3,c=kmeans.labels_,cmap='Spectral')	
	cls0 = kmeans.labels_ == 0
	cls1 = kmeans.labels_ == 1
	sorted_variables1 = asb_w_kinetics.loc[cls0, variable].sum(axis=0).sort_values(ascending=False).index[0:50]
	sorted_variables2 = asb_w_kinetics.loc[cls1, variable].sum(axis=0).sort_values(ascending=False).index[0:50]

	fig, ax = plt.subplots(ncols=2)
	ax[0].imshow(asb_w_kinetics.loc[cls0,sorted_variables1].values)
	ax[1].imshow(asb_w_kinetics.loc[cls1,sorted_variables2].values)
	
#	fig, ax = plt.subplots(ncols=2)	
#	ax[0].scatter(trans.embedding_[:, 0], trans.embedding_[:, 1], c=kmeans.labels_,cmap='Spectral')	
#	ax[1].scatter(trans1.embedding_[:, 0], trans1.embedding_[:, 1], c=kmeans.labels_,cmap='Spectral')	
	plt.show()
	

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

	res = pd.DataFrame()
	df_h2 = pd.DataFrame(columns=['kp','h2'])
	#df_h2 = pd.DataFrame(columns=['kp','h2'])

	for i, ikp in enumerate(kp_list):
		y = asb_w_kinetics[ikp].astype(float).values.reshape((-1,1))
		X = asb_w_kinetics[variable].astype(float)
		K = np.dot(X, X.T)/len(variable)
		dummy_cov = np.ones(len(X)).reshape((-1,1))
		lmm = lmm_cov.LMM(Y=y, X=dummy_cov, G=None, K=K, regressX=True, inplace=False, forcefullrank=False)
		h2 = lmm.findH2()
		df_h2.loc[i,:] = [ikp, h2['h2'][0]]
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

	kp_list = [ 'kon','koff','ksyn']
	variable = asb_w_genotype_kp['TF'].unique()
	
	asb_w_genotype_reform['allele'] = asb_w_genotype_reform[['gene','allele']].agg('-'.join, axis=1)
	mat = asb_w_genotype_reform[['allele','TF','occu']].groupby(['allele','TF'])['occu'].sum()
	mat = mat.to_frame(name='sum_occu').reset_index()
	mat['occu'] = mat['sum_occu'].apply(lambda x: 1 if x > 0 else 0)
	mat = mat.pivot(index='allele', columns='TF',values='occu').fillna(0)
	mat_kinetics = pd.merge(mat, kinetics,left_index=True, right_index=True)

	lmm(mat_kinetics, variable, kp_list, output)
	return 0

	H1_kinetics = mat_kinetics.loc[mat_kinetics['allele']=='H1_allele']
	H2_kinetics = mat_kinetics.loc[mat_kinetics['allele']=='H2_allele']
	pair_kinetics = pd.merge(H1_kinetics, H2_kinetics, on='gene', suffixes=('_m','_p'))
	fc_cols = []
	for kp in ['kon','koff','ksyn','burst_size', 'burst_frequency']:
		pair_kinetics['%s_fc'%kp] = pair_kinetics['%s_p'%kp] / pair_kinetics['%s_m'%kp]
		pair_kinetics['log(%s_fc)'%kp] = pair_kinetics['%s_fc'%kp].transform(np.log10)
		fc_cols +=  [ '%s_fc'%kp, 'log(%s_fc)'%kp ]
		fc_cols += ['gene']

	diff_cols = []
	for var in variable:
		pair_kinetics['%s_diff'%var] = pair_kinetics['%s_p'%var] - pair_kinetics['%s_m'%var]
		diff_cols += ['%s_diff'%var]
	
	print(pair_kinetics)
	lmm(pair_kinetics, diff_cols, fc_cols, output)

if __name__ == "__main__":
	main()



