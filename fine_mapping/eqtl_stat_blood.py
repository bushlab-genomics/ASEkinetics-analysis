import matplotlib.ticker as ticker
from matplotlib.colors import ListedColormap
from matplotlib.ticker import AutoMinorLocator
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib as mpl
import matplotlib
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2

from scipy import stats
import numpy as np
import pandas as pd
import re
import sys

def main():

	eqtl_input = sys.argv[1]
	output = sys.argv[2]
	
	eqtl = pd.read_csv(eqtl_input, header=0, sep='\t')
	eqtl['ID'] = eqtl[['chromosome','position','ref','alt']].astype(str).agg(':'.join, axis=1)
	
	eqtl_filter =  eqtl.loc[eqtl['pvalue'] < 1e-6]
	eqtl_summary = eqtl_filter.groupby('molecular_trait_id')['variant'].nunique().reset_index()

	eqtl_summary.to_csv('test.csv')

def main1():

	res_in = sys.argv[1]
	res = pd.read_csv(res_in, index_col=0, header=0)
	print(res['variant'].mean())

	fig, ax= plt.subplots()
	ax.hist(res['variant'], bins=200)
	ax.set_xlim(0,2000)
	ax.set_xlabel('# eQTL > 1e-6')	
	ax.set_ylabel('Count')
	plt.show()

if __name__ == "__main__":
	main1()



