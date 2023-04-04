import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib as mpl
import matplotlib
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 2

import seaborn as sns
import numpy as np
import pandas as pd
import sys

def main():

	infile = sys.argv[1]
	alpha  = float(sys.argv[2])
	output = sys.argv[3]

	df = pd.read_csv(infile, sep=':|,', header=None, engine='python', \
		names=['kp','pval','effect size','logfc_bs','logfc_bf','n_test'])
	df['logpval'] = df['pval'].transform(np.log2)
	df['logpval'] *= -1
	df['pass' ] =  df['pval'].apply(lambda x: 1 if x < alpha else 0)

	df_res = df.groupby(['logfc_bs','logfc_bf','kp'])['pass'].mean().reset_index()

	plot_cols = ['bf','bs','eqtl']
	title_cols = ['Burst frequency','Burst size','eQTL effect size']

	fig, axs = plt.subplots(ncols=len(plot_cols),figsize=(2.5*len(plot_cols),2.5),
			sharex=True, sharey=True)
	cbar_ax = fig.add_axes([1, .3, .03, .4])	
	for i, col in enumerate(plot_cols):
		df_sub = df_res[df_res['kp']==col].pivot('logfc_bf','logfc_bs','pass')
		g = sns.heatmap(df_sub, cmap='Blues', ax=axs[i], vmin=0, vmax=1, 
				cbar=i == 0, cbar_ax=None if i else cbar_ax) 
		axs[i].invert_yaxis()
		axs[i].set_title(title_cols[i],fontsize=8)
		axs[i].set_ylabel('')
		axs[i].set_xlabel(r'log2($bs$/$bs_{0}$)')
		axs[i].set_aspect('equal')
	axs[0].set_ylabel(r'log2($bf$/$bf_{0}$)')
	cbar_ax.set_title("Power")
	plt.subplots_adjust(wspace=0.06, hspace=0.06)
	plt.savefig('%s_power_%s.pdf'%(output,alpha), dpi=1000, bbox_inches = "tight")
	
if __name__ == "__main__":
	main()


