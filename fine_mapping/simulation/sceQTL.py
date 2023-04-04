import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib as mpl
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 2
import seaborn as sns
import numpy as np
import pandas as pd
import sys

def main():
	gt_in   = sys.argv[1]
	expression_in = sys.argv[2]
	output        = sys.argv[3]

	gt = pd.read_csv(gt_in,index_col=0, header=0)
	gt['alt_count'] = gt.apply(lambda row: 1 if (row['snp1']=='1|0') or (row['snp1']=='0|1') else \
                                        ( 2 if row['snp1']=='1|1' else 0), axis=1)

	expression = pd.read_csv(expression_in,index_col=0, header=0)
	cell_index = expression.columns

	expression[['sub_idx','allele']] = expression.index.to_series().str.split('-',n=1,expand=True)
	exp_grp = expression.groupby('sub_idx')[cell_index].sum().reset_index()
	exp_grp['sub_idx'] = exp_grp['sub_idx'].astype(int)

	exp_grp_gt = pd.merge(exp_grp, gt, left_on='sub_idx', right_index=True)
	exp_gt_reform = pd.melt(exp_grp_gt, id_vars='alt_count', value_vars=cell_index, \
				var_name='cell', value_name='exp')

	exp_gt_mean = exp_gt_reform.groupby('alt_count')['exp'].mean()
	exp_gt_var = exp_gt_reform.groupby('alt_count')['exp'].var()

	my_pal = {2: "#ff7a7a", 0:'#2876df',1:'#ae78a0'}

	fig, ax = plt.subplots(figsize=(2,2))
	g = sns.violinplot(x='alt_count',y='exp', data=exp_gt_reform, ax=ax,palette=my_pal)
#	sns.swarmplot(x='alt_count',y='exp', data=exp_gt_reform, ax=ax, \
#			color='white', edgecolor='black', linewidth=1, size=2)
	ax.set_title('Simulated scRNA-seq profile')
	ax.set_ylabel('UMI')
	ax.set_xlabel('Genotype')

	x=[0,1,2]
	ax.plot(x, exp_gt_mean, 'k.:',lw=2, label='mean', alpha=0.8, markersize=8)
	ax.plot(x, exp_gt_var,  'k.--', lw=2, label='variance', alpha=0.8, markersize=8)
	ax.legend(frameon=False, bbox_to_anchor=(1, 1.01))
	plt.savefig('%s.png'%output, dpi=300, bbox_inches = "tight")
	plt.show()
if __name__ == "__main__":
	main()

