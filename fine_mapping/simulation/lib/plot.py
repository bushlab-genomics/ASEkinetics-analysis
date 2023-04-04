import matplotlib as mpl
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 2
import pandas as pd
import numpy as np
import sys
	
def main():

	res_input = sys.argv[1]
	res_df = pd.read_csv(res_input, header=None, sep=' ',names=['maf', 'sample', 'slope',  'pval'])

	res_df['logpval'] = res_df['pval'].transform(np.log10)
	res_df['logpval'] *= -1

	sns.color_palette("mako", as_cmap=True)
	res_reform = res_df.pivot(index='sample', columns='maf', values='logpval')	
	ax = sns.heatmap(res_reform,  cbar_kws={'label': 'log10(pval)'}, \
			cmap='mako')#, vmin=0, vmax=20)
	ax.invert_yaxis()
	plt.savefig(res_input + '.pdf')
	plt.show()

if __name__ == "__main__":

	main()

