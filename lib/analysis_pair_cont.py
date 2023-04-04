import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.cm import get_cmap
import matplotlib as mpl
from pygam import LinearGAM,s
import matplotlib
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 2
matplotlib.use('Agg')

import seaborn as sns
sns.set_style("white")
sns.despine(top=False,right=False,left=False)

from scipy import stats
from scipy.sparse import coo_matrix
from statsmodels.stats.multitest import fdrcorrection
from sklearn import linear_model
import numpy as np
import pandas as pd
import pingouin as pg
import re
import sys

def matrix(indf):

	TF_list = indf['TF'].unique().tolist()
	allele_list=indf['allele'].unique().tolist()
	for idx, allele in enumerate(allele_list):
		select = indf['allele'] == allele
		indf.loc[select, 'allele_idx'] = int(idx)

	for jdx, TF in enumerate(TF_list):
		select = indf['TF'] == TF
		indf.loc[select, 'TF_idx'] = int(jdx)

	ncol=len(TF_list)
	nrow=len(allele_list)
	col = indf['TF_idx']
	row = indf['allele_idx']
	data = indf['AR']
	coo = coo_matrix((data, (row, col)), shape=(nrow, ncol))
	tab = pd.DataFrame(coo.toarray(),index=allele_list,columns=TF_list)

	return tab

def pref(x):

	if not x['peak.isASB']: return  1,1

	if   (x['Corrected.AR'] > 0.5) & (x['genotype'] == '1|0'):
		return  0,1
	elif (x['Corrected.AR'] > 0.5) & (x['genotype'] == '0|1'):
		return  1,0
	elif (x['Corrected.AR'] < 0.5) & (x['genotype'] == '1|0'):
		return  1,0
	elif (x['Corrected.AR'] < 0.5) & (x['genotype'] == '0|1'):
		return  0,1

def AR_pref(x):
	if   (x['Corrected.AR'] > 0.5) & (x['genotype'] == '1|0'):
		return  1 - x['Corrected.AR'], x['Corrected.AR']
	elif (x['Corrected.AR'] > 0.5) & (x['genotype'] == '0|1'):
		return  x['Corrected.AR'], 1 - x['Corrected.AR']
	elif (x['Corrected.AR'] < 0.5) & (x['genotype'] == '1|0'):
		return  1 - x['Corrected.AR'], x['Corrected.AR']
	elif (x['Corrected.AR'] < 0.5) & (x['genotype'] == '0|1'):
		return  x['Corrected.AR'], 1 - x['Corrected.AR']

def single_tf_di(tf, asb_w_kinetics_all, output):
	
	asb_w_kinetics = asb_w_kinetics_all[asb_w_kinetics_all['TF']==tf]
	kp=['kon','koff','ksyn','burst size','burst frequency']

	asb_w_kinetics.to_csv(output+'_'+tf+'.csv',index=False)
	nfig=len(kp)
	fig,ax = plt.subplots(ncols=nfig, sharex=True,constrained_layout=True,figsize=(3*nfig, 3))
	for i, ikp in enumerate(kp):
		x  = [0,1]
		bool_w_AR = (asb_w_kinetics['AR'] > 0)
		subdf = asb_w_kinetics[['peak_name','TF','AR',ikp,'allele','gene']].drop_duplicates()
		subdf_w_AR = subdf.loc[bool_w_AR]
		subdf_wo_AR = subdf.loc[~bool_w_AR]
		
		subdf_reform = pd.merge(subdf_w_AR, subdf_wo_AR, on=['gene','peak_name'],suffixes=('','_0'))
		rvs1 = subdf_reform[ikp]
		rvs0 = subdf_reform[ikp+'_0']

		if len(rvs1) < 6 or len(rvs0) < 6: continue	
		wx_res = pg.wilcoxon(x=rvs1,y=rvs0,alternative='two-sided')

		pval = wx_res.loc['Wilcoxon','p-val']
		pval = round(pval,5)

		es = wx_res.loc['Wilcoxon','CLES']
		es -= 0.5
		es = round(es,5)

		ax[i].scatter(np.zeros(len(rvs0)),rvs0,s=1)
		ax[i].scatter(np.ones(len(rvs1)), rvs1,s=1)

		for ii in range(len(rvs0)):
			x = [0,1]
			y = [rvs0[ii], rvs1[ii]]
			ax[i].plot(x, y, color="black", alpha=0.1)

		ax[i].set_title(ikp)
		ax[i].set_xlabel('TF binding status')
		ax[i].set_xticks([0,1])
		ax[i].set_xlim(-1,2)
		ax[i].text(0.5, 0.8,'U=%s\npval=%s'%(es, pval), \
			horizontalalignment='center',verticalalignment='center',transform = ax[i].transAxes)
	fig.set_figwidth(10)
	plt.savefig(output + '_' + tf + '.pdf')


def lr_analysis(asb_w_kinetics,variable,kp_list,out):

	res = pd.DataFrame()

	for i, ikp in enumerate(kp_list):
		y = asb_w_kinetics[ikp].astype(float)
		X = asb_w_kinetics[variable].astype(float)

		gam = LinearGAM().fit(X, y)
		gam.summary()
		print(X.shape)
		print(len(gam.coef_))
		print(gam.coef_)

def MWU_fct_analysis(asb_w_kinetics, variable, kp_list, out):

	res = pd.DataFrame(columns=['kp','fct','pval','effect size'])
	
	i = 0
	for kp in kp_list:
		for fct in variable:
			bool_w_fct = (asb_w_kinetics['TF'] == fct)
			bool_w_AR = (asb_w_kinetics['AR'] > 0)
			subdf = asb_w_kinetics.loc[bool_w_fct,['peak_name','TF','AR',kp,'allele','gene']].drop_duplicates()
			subdf_w_AR = subdf.loc[bool_w_AR]
			subdf_wo_AR = subdf.loc[~bool_w_AR]

			subdf_reform = pd.merge(subdf_w_AR, subdf_wo_AR, on=['gene','peak_name'],suffixes=('','_0'))
			rvs1 = subdf_reform[kp]
			rvs0 = subdf_reform[kp+'_0']

			if len(rvs1) < 10 or len(rvs0) < 10: continue
			wx_res = pg.wilcoxon(x=rvs1,y=rvs0,alternative='two-sided')

			res.loc[i, 'kp' ] = kp
			res.loc[i, 'fct'] = fct
			res.loc[i, 'pval'] = wx_res.loc['Wilcoxon','p-val']
			res.loc[i, 'effect size'] = wx_res.loc['Wilcoxon','RBC']
			i += 1
	for kp in kp_list:
		kp_bool = (res['kp']==kp)
		p_vals = res.loc[kp_bool,'pval'].values.astype(float)
		rejected, pvalue_corrected = pg.multicomp(p_vals,alpha=0.1,method='fdr_bh')
		res.loc[kp_bool,'pval-corrected'] = pvalue_corrected

	res_filter_fct = res.loc[res['pval']<0.05,'fct'].unique()
	res_filter = res[ res['fct'].isin(res_filter_fct)]

	select = res_filter['kp']=='ksyn'
	order = res_filter.loc[select].sort_values(by='effect size')['fct'].unique()
	
	np.savetxt(out+'_reg.list',order,fmt='%s')
	g = sns.FacetGrid(res_filter, col="kp",despine=False)
	g.map_dataframe(sns.barplot,x='effect size',y='fct',order=order, color='mediumpurple', orient='h', dodge=False, errwidth=0)
	g.fig.subplots_adjust(wspace=0, hspace=0)
	for ax in g.axes.flat:
		ax.grid(True, axis='x',linestyle='--')
		ax.grid(True, axis='y',linestyle='--')

	xlimit = res_filter['effect size'].apply(lambda x: abs(x)).max()

	g.set(xlim=(-0.5,0.5))
	g.map(plt.axvline,x=0,color='k',linewidth=1)
	g.set_titles('{col_name}')
	g.fig.set_figwidth(10)
	g.fig.set_figheight(np.ceil(len(res_filter)/24))
	plt.savefig(out+'_mwu.pdf',dpi=300,bbox_inches ='tight')
	return res_filter

def allele_ratio(asb_w_genotype_reform,out):
	
	ar = asb_w_genotype_reform.groupby(['TF','allele'])['AR'].sum()
	ar_percent = ar / asb_w_genotype_reform.groupby(['TF','allele'])['AR'].count()

	ar_df = ar.to_frame().reset_index()
	ar_df_pivot = ar_df.pivot(index='TF', columns='allele', values='AR')
	ar_df_pivot.reset_index(inplace=True)

	ar_df_pivot.sort_values(by='maternal_allele', inplace=True)	
	labels=ar_df_pivot['TF']
	m_ar = ar_df_pivot['maternal_allele']
	p_ar = ar_df_pivot['paternal_allele']

	ar_percent_df = ar_percent.to_frame().reset_index()
	ar_percent_df_pivot = ar_percent_df.pivot(index='TF', columns='allele', values='AR')
	ar_percent_df_pivot.reset_index(inplace=True)

	ar_percent_df_pivot.set_index('TF',inplace=True)	
	m_ar_percent = ar_percent_df_pivot.loc[labels,'maternal_allele']
	p_ar_percent = ar_percent_df_pivot.loc[labels,'paternal_allele']
		
	fig, ax = plt.subplots(figsize=(8,10),constrained_layout=True, ncols=2)
	ax[0].barh(labels, m_ar, label='maternal')
	ax[0].barh(labels, p_ar, left=m_ar,label='paternal')

	ax[1].barh(labels, m_ar_percent, label='maternal')
	ax[1].barh(labels, p_ar_percent, left=m_ar_percent,label='paternal')
	ax[1].legend(frameon=False, bbox_to_anchor=[1.01,0.5])

	plt.savefig(out+'_ASB_dist.pdf',dpi=300,bbox_inches ='tight')
	plt.show()

def main():

	asb_input = sys.argv[1]
	vcf_input = sys.argv[2]
	kinetics_input = sys.argv[3]
	output = sys.argv[4]
	
	asb = pd.read_csv(asb_input, header=0,index_col=0)
	asb = asb.loc[asb['peak.isASB']==True]
	vcf = pd.read_csv(vcf_input, header=None, sep='\t',usecols=list(range(6)) + [9,11,12,13],\
		names=['chr','pos','ID','ref','alt','qual','genotype','reg_starts','reg_ends','gene'])
	kinetics = pd.read_csv(kinetics_input, header=0,index_col=0)
	kinetics[['gene','allele']] = kinetics.index.to_series().str.split('-', n=1, expand=True)
	kinetics['burst size'] = kinetics['ksyn']/kinetics['koff']
	kinetics['burst frequency'] = kinetics['kon'] * kinetics['koff'] / (kinetics['kon'] + kinetics['koff'])
	kinetics['fraction'] = kinetics['kon'] / (kinetics['kon'] + kinetics['koff'])

	asb_w_genotype = pd.merge(asb,vcf,on='ID')
	asb_w_genotype[['maternal_AR','paternal_AR']] = asb_w_genotype.apply(AR_pref, axis=1, result_type="expand")
	
	asb_w_genotype_reform=pd.melt(asb_w_genotype,id_vars=['peak_name','TF','gene'],\
		value_vars=['maternal_AR','paternal_AR'],var_name='allele',value_name='AR')
	asb_w_genotype_reform['allele']  = asb_w_genotype_reform['allele'].str.replace('_AR','_allele')
	asb_w_genotype_reform.drop_duplicates(inplace=True)
	asb_w_genotype_kp = pd.merge(asb_w_genotype_reform, kinetics,on=['gene','allele'])


	kp_list = ['kon','koff','ksyn','burst size','burst frequency']
	variable = asb_w_genotype_kp['TF'].unique()
	
	asb_w_genotype_reform['allele'] = asb_w_genotype_reform[['gene','allele']].agg('-'.join, axis=1)
	mat = matrix(asb_w_genotype_reform[['allele','TF','AR']])
	print(mat)
	mat_kinetics = pd.merge(mat, kinetics,left_index=True, right_index=True)

	lr_analysis(mat_kinetics,variable,kp_list,output)

if __name__ == "__main__":
	main()



