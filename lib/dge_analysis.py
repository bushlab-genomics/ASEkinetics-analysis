import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2

from scipy import stats
import pingouin as pg
import pandas as pd
import seaborn as sns
import numpy as np
import pickle
import sys
import re

def violin_analysis_allele(ase_gene, cols, gene, output):
	
	fig,ax = plt.subplots(constrained_layout=True,figsize=(6,3),ncols=2,sharex=True,sharey=True)
	for i, row in ase_gene.iterrows():
		ase_val = row[cols].astype(int).values
		if len(ase_val)<1 or (max(ase_val) < 2): continue

		if row['allele'] == 'H1_allele':
			color= '#2876df'
			xpos = 0
		elif row['allele'] == 'H2_allele':
			color='#ff7a7a'
			xpos = 1
		if row['sample'] == 'NA12878':
			ax[0].set_title('NA12878',loc='right')
			violin_plot1 = ax[0].violinplot(ase_val, [xpos], showmeans=True)
			for pc in violin_plot1['bodies']:
				pc.set_facecolor(color)
				pc.set_edgecolor(color)

		elif row['sample'] == 'NA18502':
			ax[1].set_title('NA18502',loc='right')
			violin_plot2 = ax[1].violinplot(ase_val, [xpos], showmeans=True)
			for pc in violin_plot2['bodies']:
				pc.set_facecolor(color)
				pc.set_edgecolor(color)

	ax[0].set_ylabel('UMI')
	ax[0].set_xticks([0,1])
	ax[1].set_xticks([0,1])
	ax[0].set_xticklabels(['H1-Haplotype','H2-Haplotype'])
	ax[1].set_xticklabels(['H1-Haplotype','H2-Haplotype'])

	plt.savefig('%s_%s_hist.pdf'%(output, gene))
	plt.show()

def hist_analysis_allele(ase_gene, cols, gene, output):
	
	fig,ax = plt.subplots(constrained_layout=True,figsize=(8,4),ncols=2,sharex=True,sharey=True)
	for i, row in ase_gene.iterrows():
		ase_val = row[cols].astype(int).values
		if len(ase_val)<1 or (max(ase_val) < 2): continue
		gene_density = stats.gaussian_kde(ase_val)
		x_density = np.arange(0,max(ase_val),0.8)

		ase_mean = np.mean(ase_val)
		ase_std = np.std(ase_val)

		if row['allele'] == 'H1_allele':
			color= '#7a7aff'
			ypos= 0.3
		elif row['allele'] == 'H2_allele':
			color='#ff7a7a'
			ypos= 0.35
		if row['sample'] == 'NA12878':
			ax[0].set_title('NA12878',loc='right')
			ax[0].plot(x_density, gene_density(x_density), color=color, ls='-', lw=2)
			ax[0].axvline(ase_mean, color=color, ls='--', lw=1)
			ax[0].plot((ase_mean-ase_std, ase_mean+ase_std),(ypos,ypos),'r|-',color=color,lw=2)
			ax[0].scatter(ase_mean, ypos, s=8,color=color)

		elif row['sample'] == 'NA18502':
			ax[1].set_title('NA18502',loc='right')
			ax[1].plot(x_density, gene_density(x_density), color=color, ls='-', lw=2,alpha=0.6)
			ax[1].vlines(ase_mean, ymin=0, ymax=gene_density(ase_mean), color=color, ls='--', lw=1)
			ax[1].plot((ase_mean-ase_std, ase_mean+ase_std),(ypos,ypos),'r|-',color=color,lw=2)
			ax[1].scatter(ase_mean, ypos, s=8, color=color)

	ax[0].set_ylabel('Probability density')
	ax[0].set_xlabel('UMI')
	ax[1].set_xlabel('UMI')
	ax[1].set_ylim(-0.05,0.4)

	plt.savefig('%s_%s_hist.pdf'%(output, gene))
	plt.show()	

def main():

	ase_list  = sys.argv[1].split(',')
	sample_list   = sys.argv[2].split(',')
	gene       = sys.argv[3]
	output     = sys.argv[4]

	ase_gen = []
	ase_all = pd.DataFrame()
	for i, ase_temp_in in enumerate(ase_list):
		ase_temp = pd.read_csv(ase_temp_in, header=0, index_col=0).dropna()
		cols = ase_temp.columns
		ase_temp[['gene','allele']] = ase_temp.index.to_series().str.split('-',expand=True)
		stat_file = re.sub('.ase.reform', '.qc.stats', ase_temp_in)
		ase_stat = pd.read_csv(stat_file, header=0)
		gene_qc = ase_stat.loc[ase_stat['qc']==True,'gene'].unique()
		ase_temp.loc[ase_temp['gene'].isin(gene_qc),'sample'] = sample_list[i]

		ase_gene = ase_temp.loc[ase_temp['gene']==gene]
		ase_all = ase_all.append(ase_gene,ignore_index=False)

#	tempfile = open("%s.obj",'wb')
#	pickle.dump(ase_all, tempfile)
#	pickle.dump(cols, tempfile)
#	pickle.dump(gene, tempfile)
#	pickle.dump(output, tempfile)
#
#	tempfile = open("%s.obj",'rb')
#	ase_all = pickle.load(tempfile)
#	cols = pickle.load(tempfile)
#	gene = pickle.load(tempfile)
#	output = pickle.load(tempfile)
		
	violin_analysis_allele(ase_all, cols, gene, output)

if __name__ == "__main__":
	main()


