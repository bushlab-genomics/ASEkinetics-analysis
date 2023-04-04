import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib as mpl
import matplotlib
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 2
import seaborn as sns
import numpy as np
import pandas as pd
import sys

def main():
	infile = sys.argv[1]
	alpha  = float(sys.argv[2])
	output = sys.argv[3]
	df = pd.read_csv(infile, sep=':|,', header=None, engine='python', index_col=0,\
		names=['kp','pval','effect size','n_sample','n_cell','ntest'])
	df.loc[:,'pass'] = 0
	select = (df['pval']<alpha)
	df.loc[select, 'pass'] = 1

	df_res = df.groupby(['n_sample','n_cell','kp'])['pass'].mean().reset_index()
	plot_cols = ['bf','bs','eqtl']
	title_cols = ['Burst frequency','Burst size','eQTL effect size']
	
	fig, axs = plt.subplots(ncols=len(plot_cols),figsize=(2.5*len(plot_cols),2),
				sharex=True, sharey=True)
	cbar_ax = fig.add_axes([1, .3, .03, .4])
	for i, col in enumerate(plot_cols):
		df_sub = df_res[df_res['kp']==col].pivot('n_cell','n_sample','pass')
		g = sns.heatmap(df_sub, cmap='Blues', ax=axs[i], vmin=0, vmax=1,\
			cbar=i == 0, cbar_ax=None if i else cbar_ax , annot=True, annot_kws={"fontsize":8})
		#for _, spine in g.spines.items():
		#	spine.set_visible(True)
		axs[i].invert_yaxis()
		axs[i].set_title(title_cols[i])
		axs[i].set_ylabel('# Cells')
		axs[i].set_xlabel('# Subjects')
		axs[i].set_aspect('equal')
#	axs[0].set_ylabel('# cells')
	cbar_ax.set_title("Power")

#	plt.subplots_adjust(wspace=0.01, hspace=0.01)
	plt.savefig('%s_%s.png'%(output,alpha), dpi=300, bbox_inches = "tight")
	plt.show()

if __name__ == "__main__":
	main()

