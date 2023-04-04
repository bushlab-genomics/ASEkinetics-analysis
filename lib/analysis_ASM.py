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
#matplotlib.use('Agg')

import seaborn as sns
sns.set_style("white")
sns.despine(top=False,right=False,left=False)

from scipy import stats
from scipy.sparse import coo_matrix
from statsmodels.stats.multitest import fdrcorrection
import numpy as np
import pandas as pd
import pingouin as pg
import re
import sys

def single_tf_di(tf, asb_w_kinetics, output):

	kp=['kon','koff','ksyn','burst size','burst frequency']

	asb_w_kinetics.to_csv(output+'_'+tf+'.csv',index=False)
	nfig=len(kp)
	fig,ax = plt.subplots(ncols=nfig, sharex=True,constrained_layout=True,figsize=(3*nfig, 3))
	for i, ikp in enumerate(kp):
		x  = [0,1]
		bool_w_occu = (asb_w_kinetics['occu'] > 0)
		subdf = asb_w_kinetics[['CG_position','occu',ikp,'allele','gene']].drop_duplicates()
		subdf_w_occu = subdf.loc[bool_w_occu]
		subdf_wo_occu = subdf.loc[~bool_w_occu]

		subdf_reform = pd.merge(subdf_w_occu, subdf_wo_occu, on=['gene','CG_position'],suffixes=('','_0'))

		rvs1 = subdf_reform[ikp]
		rvs0 = subdf_reform[ikp+'_0']

		if len(rvs1) < 6 or len(rvs0) < 6: continue
		wx_res = pg.wilcoxon(x=rvs1,y=rvs0)

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


def MWU_fct_analysis(asb_w_kinetics, kp_list, out):

	res = pd.DataFrame(columns=['kp','fct','pval','effect size'])

	i = 0
	for kp in kp_list:
		bool_w_occu = (asb_w_kinetics['occu'] > 0)
		subdf = asb_w_kinetics[['CG_position','occu',kp,'allele','gene']].drop_duplicates()
		subdf_w_occu = subdf.loc[bool_w_occu]
		subdf_wo_occu = subdf.loc[~bool_w_occu]
		subdf_reform = pd.merge(subdf_w_occu, subdf_wo_occu, on=['gene','CG_position'],suffixes=('','_0'))
		rvs1 = subdf_reform[kp]
		rvs0 = subdf_reform[kp+'_0']

		if len(rvs1) < 6 or len(rvs0) < 6: continue
		wx_res = pg.wilcoxon(x=rvs1,y=rvs0)

		res.loc[i, 'kp' ] = kp
		res.loc[i, 'fct'] = 'met'
		res.loc[i, 'pval'] = wx_res.loc['Wilcoxon','p-val']
		res.loc[i, 'effect size'] = wx_res.loc['Wilcoxon','RBC']
		i += 1
	for kp in kp_list:
		kp_bool = (res['kp']==kp)
		p_vals = res.loc[kp_bool,'pval'].values.astype(float)
		rejected, pvalue_corrected = pg.multicomp(p_vals,alpha=0.1,method='fdr_bh')
		res.loc[kp_bool,'pval-corrected'] = pvalue_corrected

	res_filter_fct = res.loc[res['pval-corrected']<0.1,'fct'].unique()
	res_filter_fct = res['fct'].unique()
	res_filter = res[ res['fct'].isin(res_filter_fct)]
	order = res_filter.sort_values(by='pval-corrected')['fct'].unique()

	np.savetxt(out+'_reg.list',order,fmt='%s')
	g = sns.FacetGrid(res_filter, col="kp",despine=False)
	g.map_dataframe(sns.barplot,x='effect size',y='fct',order=order, color='mediumpurple', orient='h', dodge=False, errwidth=0)
	g.fig.subplots_adjust(wspace=0, hspace=0)
	for ax in g.axes.flat:
		ax.grid(True, axis='x',linestyle='--')
		ax.grid(True, axis='y',linestyle='--')

	xlimit = res_filter['effect size'].apply(lambda x: abs(x)).max()

	g.set(xlim=(-1,1))
	g.map(plt.axvline,x=0,color='k',linewidth=1)
	g.set_titles('{col_name}')
	g.fig.set_figwidth(10)
	g.fig.set_figheight(np.ceil(len(res_filter)/24))
	plt.savefig(out+'_mwu.pdf',dpi=300,bbox_inches ='tight')
	return res_filter

def pref(x):

	if   (x['genotype'] == '0|1'):
		return  x['Reference_allele_fraction_methylated'],x['Variant_allele_fraction_methylated']
	elif (x['genotype'] == '1|0'):
		return  x['Variant_allele_fraction_methylated'],x['Reference_allele_fraction_methylated']

def main():

	ASM_input = sys.argv[1]
	vcf_input = sys.argv[2]
	kinetics_input = sys.argv[3]
	output = sys.argv[4]
	
	ASM = pd.read_csv(ASM_input, header=0, sep='\t')
	ASM['chr'] = ASM['Chromosome'].str.strip('chr')
	ASM['loc'] = ASM[['chr','SNP_postion_ends']].astype(str).agg(':'.join, axis=1)	

	vcf = pd.read_csv(vcf_input, header=None, sep='\t',usecols=list(range(6)) + [9,11,12,13],\
		names=['chr','pos','ID','ref','alt','qual','genotype','reg_starts','reg_ends','gene'])
	vcf['loc'] = vcf[['chr','pos']].astype(str).agg(':'.join, axis=1)	
 
	kinetics = pd.read_csv(kinetics_input, header=0,index_col=0)
	kinetics[['gene','allele']] = kinetics.index.to_series().str.split('-', n=1, expand=True)
	kinetics['burst size'] = kinetics['ksyn']/kinetics['koff']
	kinetics['burst frequency'] = kinetics['kon'] * kinetics['koff'] / (kinetics['kon'] + kinetics['koff'])
	kinetics['fraction'] = kinetics['kon'] / (kinetics['kon'] + kinetics['koff'])

	ASM_w_genotype = pd.merge(ASM,vcf,on='loc')
	ASM_w_genotype[['maternal_met','paternal_met']] = ASM_w_genotype.apply(pref, axis=1, result_type="expand")

	ASM_w_genotype_reform=pd.melt(ASM_w_genotype,id_vars=['gene','CG_position'],\
		value_vars=['maternal_met','paternal_met'],var_name='allele',value_name='met')
	ASM_w_genotype_reform['allele']  = ASM_w_genotype_reform['allele'].str.replace('_met','_allele')
	ASM_w_genotype_reform[['met_percent','CG_num']] = ASM_w_genotype_reform['met'].str.split('(',expand=True)
	ASM_w_genotype_reform['met_percent'] = ASM_w_genotype_reform['met_percent'].astype(float)
	ASM_w_genotype_reform['occu'] = ASM_w_genotype_reform.apply(lambda row: 1 if row['met_percent'] > 0.5 else 0, axis=1)
	
	ASM_w_genotype_reform.drop_duplicates(inplace=True)
	ASM_w_genotype_kp = pd.merge(ASM_w_genotype_reform, kinetics,on=['gene','allele'])

	kp_list = ['kon','koff','ksyn','burst size','burst frequency']
	mwu_res = MWU_fct_analysis(ASM_w_genotype_kp,kp_list,output)
	mwu_res.sort_values(by='pval-corrected').to_csv(output+'_mwu_sig_gene.csv', index=False)

	single_tf_di('met', ASM_w_genotype_kp, output)

if __name__ == "__main__":
	main()



