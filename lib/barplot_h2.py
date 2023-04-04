import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2
import pandas as pd
import numpy as np
import sys

def main():

	h2_input1 = sys.argv[1]
	df1_h2 = pd.read_csv(h2_input1, header=0)
	h2_input2 = sys.argv[2]
	df2_h2 = pd.read_csv(h2_input2, header=0)

	discrete_x = df2_h2.index	
	x1 = discrete_x - 0.2
	x2 = discrete_x + 0.2

	fig, ax = plt.subplots(constrained_layout=True,figsize=(7, 4))
	ax.bar(x1, df1_h2['mean'], width=0.3, fill=False, hatch='///')
	ax.errorbar(x1, df1_h2['mean'], yerr=df1_h2['std'], fmt='s',capsize=4,c='k',lw=2)
	ax.bar(x2, df2_h2['mean'], width=0.3, fill=False, hatch='.')
	ax.errorbar(x2, df2_h2['mean'], yerr=df2_h2['std'], fmt='s',capsize=4,c='k',lw=2)

	circ1 = mpatches.Patch(facecolor='w', edgecolor='k', hatch=r'\\\\',label='core promoter region')
	circ2 = mpatches.Patch(facecolor='w', edgecolor='k', hatch=r'.',label='Distal enhancer profile')
	ax.legend(handles = [circ1,circ2],frameon=False, bbox_to_anchor=(0.6, 1.01))	
	ax.set_xticks(discrete_x)
	ax.set_xticklabels(['k+','k-','r','burst size','burst frequency'])
	ax.set_yticks(np.arange(0,1.1,0.1))
	ax.set_ylabel('Proportion of variance in transcription kinetics')
	plt.tight_layout()
	plt.savefig('h2_comparsion' + '.pdf')
	plt.show()

if __name__ == "__main__":
	main()

