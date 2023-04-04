import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.cm import get_cmap
import matplotlib as mpl
import matplotlib
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 2
matplotlib.use('Agg')

import seaborn as sns
sns.set_style("white")
sns.despine(top=False,right=False,left=False)

from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
from sklearn import linear_model
import numpy as np
import pandas as pd
import pingouin as pg
import re
import sys

def MWU_fct_analysis(qbic_w_kp, variable, kp_list, out):

	res = pd.DataFrame(columns=['kp','fct','pval','effect size'])
	i = 0

	for kp in kp_list:
		for fct in variable:
			bool_w_fct = (qbic_w_kp['TF_gene'] == fct)
			bool_w_pref = (qbic_w_kp['binding_count'] > 0)

			subdf = qbic_w_kp.loc[bool_w_fct,['ID_pos','TF_gene','binding_count',kp,'allele','gene']].drop_duplicates()
			subdf_w_pref = subdf.loc[bool_w_pref]
			subdf_wo_pref = subdf.loc[~bool_w_pref]

			subdf_reform = pd.merge(subdf_w_pref, subdf_wo_pref, on=['ID_pos','gene'],suffixes=('','_0'))
			rvs1 = subdf_reform[kp]
			rvs0 = subdf_reform[kp+'_0']
			
			if len(rvs1) < 10 or len(rvs0) < 10: continue
			wx_res = pg.wilcoxon(rvs1,rvs0,alternative='two-sided')
			res.loc[i, 'kp' ] = kp
			res.loc[i, 'fct'] = fct
			res.loc[i, 'pval'] = wx_res.loc['Wilcoxon','p-val']
			res.loc[i, 'effect size'] = wx_res.loc['Wilcoxon','RBC']
			i += 1

	for kp in kp_list:
		kp_bool = (res['kp']==kp)
		rejected, pvalue_corrected = fdrcorrection(res.loc[kp_bool,'pval'])
		res.loc[kp_bool,'pval-corrected'] = pvalue_corrected

	res_filter = res
	order = variable	
#	order = res_filter.sort_values(by='pval-corrected')['fct'].unique()

	np.savetxt(out+'_reg.list',order,fmt='%s')
	g = sns.FacetGrid(res_filter, col="kp",despine=False)
	g.fig.subplots_adjust(wspace=0, hspace=0)
	g.map_dataframe(sns.barplot,x='effect size',y='fct',order=order, color='dodgerblue', orient='h', dodge=False, errwidth=0)
	for ax in g.axes.flat:
		ax.grid(True, axis='x',linestyle='--')
		ax.grid(True, axis='y',linestyle='--')

	g.set(xlim=(-1,1))
	g.map(plt.axvline,x=0,color='k',linewidth=1)
	g.set_titles('{col_name}')
	g.fig.set_figwidth(10)
	g.fig.set_figheight(np.ceil(len(res_filter)/24))
	plt.savefig(out+'_mwu.pdf',dpi=300,bbox_inches ='tight')
	return res_filter

def single_tf_di(tf, qbic_subset, output):

	qbic_tf = qbic_subset[qbic_subset['TF_gene']==tf]
	kp=['kon','koff','ksyn','burst size','burst frequency']
	qbic_tf.to_csv(output+'_'+tf+'.csv',index=False)

	nfig=len(kp)
	fig,ax = plt.subplots(ncols=nfig, sharex=True,constrained_layout=True,figsize=(3*nfig, 3))
	for i, ikp in enumerate(kp):
		bool_w_fct = (qbic_tf['TF_gene'] == tf)
		bool_w_pref = (qbic_tf['binding_count'] > 0)

		subdf = qbic_tf.loc[bool_w_fct,['ID_pos','TF_gene','binding_count',ikp,'allele','gene']].drop_duplicates()
		subdf_w_pref = subdf.loc[bool_w_pref]
		subdf_wo_pref = subdf.loc[~bool_w_pref]

		subdf_reform = pd.merge(subdf_w_pref, subdf_wo_pref, on=['gene','ID_pos'],suffixes=('','_0'))
		rvs1 = subdf_reform[ikp]
		rvs0 = subdf_reform[ikp+'_0']
			
		if len(rvs1) < 10 or len(rvs0) < 10: continue
		wx_res = pg.wilcoxon(rvs1,rvs0,alternative='two-sided')
		pval = wx_res.loc['Wilcoxon','p-val']
		pval = round(pval,5)
		es = wx_res.loc['Wilcoxon','RBC']
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

def scatter(tf, qbic_subset, output):

	qbic_tf = qbic_subset.loc[qbic_subset['TF_gene']==tf]
	if len(qbic_tf) <2 : return 0

	kp_list = ['kon','koff','ksyn','burst size','burst frequency']
	cols = len(kp_list)    
	fig,axs = plt.subplots(ncols=cols,constrained_layout=True, figsize=(3*cols, 3))

	for i, ikp in enumerate(kp_list):
		x = qbic_subset[['gene','allele',ikp,'binding_zscore']].drop_duplicates()['binding_zscore'].values
		y = qbic_subset[['gene','allele',ikp,'binding_zscore']].drop_duplicates()[ikp].values
	
		x = x.reshape((-1, 1))
		reg = linear_model.LinearRegression().fit(x, y)
		m = reg.coef_[0]
		m = round(m, 5)
		c = reg.intercept_
		r_sq = reg.score(x, y)
		r_sq = round(r_sq,5)

		axs[i].scatter(x,y,c='k',s=12, marker="^")
		axs[i].set_xlabel('zscore')
		axs[i].set_xlim([-1,20])
		axs[i].set_ylabel(ikp)
		axs[i].set_title('effect size=%s\nR2=%s'%(m,r_sq),loc='right')

		# fitting curve
		x1 = np.linspace(1,15,100)
		y1 = m*x1 + c
		axs[i].plot(x1,y1, '--',c='dodgerblue',label='fitting curve')
	
	fig.set_figwidth(10)		
	plt.savefig(output + '_allele_es' + '.pdf')

def main():

	qbic_input = sys.argv[1]
	kinetics_input = sys.argv[2]
	tf_list_in = sys.argv[3]
	tf_list = np.loadtxt(tf_list_in, dtype=str)
	output = sys.argv[4]

	qbic = pd.read_csv(qbic_input,header=0)	
	qbic[['maternal_allele','paternal_allele']] = qbic['genotype'].str.split('|',expand=True)

#	snp_list_in = sys.argv[5]
#	snp_list = np.loadtxt(snp_list_in, dtype=str)
#	qbic = qbic.loc[qbic['ID_pos'].isin(snp_list)]

	kinetics = pd.read_csv(kinetics_input, header=0,index_col=0)
	kinetics.dropna(inplace=True)
	kinetics[['gene','allele']] = kinetics.index.to_series().str.split('-', n=1, expand=True)

	kp_cols = ['kon','koff','ksyn']
	kinetics['burst size'] = kinetics['ksyn']/kinetics['koff']
	kinetics['burst frequency'] = kinetics['kon'] * kinetics['koff'] / (kinetics['kon'] + kinetics['koff'])
	kinetics['fraction'] = kinetics['kon'] / (kinetics['kon'] + kinetics['koff'])

	pivot_cols = ['maternal_allele','paternal_allele']
	id_vars = set(qbic.columns.tolist()) - set(pivot_cols)
	qbic_reform = pd.melt(qbic, id_vars = id_vars,value_vars = pivot_cols,var_name = 'allele', value_name = 'allele_count')

	qbic_reform.drop_duplicates(inplace=True, ignore_index=True)
	qbic_reform = qbic_reform[['gene','allele', 'binding_status', 'allele_count', 'ID_pos','TF_gene']]
	qbic_reform_all = pd.merge(qbic_reform, kinetics, on=['gene','allele'])
	qbic_reform_all.drop_duplicates(inplace=True, ignore_index=True)
	qbic_reform_w_kp = qbic_reform_all.dropna()

	qbic_subset = qbic_reform_w_kp[qbic_reform_w_kp['binding_status'].isin(['unbound>bound','bound>unbound'])]
	select = qbic_subset['binding_status'] == 'unbound>bound'
	select1 = qbic_subset['binding_status'] == 'bound>unbound'

	qbic_subset.loc[select,'binding_count']  = qbic_subset.loc[select, 'allele_count'].astype(int)
	qbic_subset.loc[select1,'binding_count'] = 1 - qbic_subset.loc[select1, 'allele_count'].astype(int)

	kp_list = ['kon','koff','ksyn','burst size','burst frequency']
	variable = tf_list
	mwu_res = MWU_fct_analysis(qbic_subset, variable, kp_list, output)
	mwu_res.to_csv(output+'_mwu_sig_gene.csv', index=False)
	
	for tf in mwu_res['fct'].unique():
		single_tf_di(tf, qbic_subset, output)

if __name__ == "__main__":
	main()

