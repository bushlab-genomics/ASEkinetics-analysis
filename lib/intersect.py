from matplotlib.colors import ListedColormap
from scipy.sparse import coo_matrix
from sklearn import linear_model
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
import pingouin as pg
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 2
import seaborn as sns
sns.set_style("white")
sns.despine(top=False,right=False,left=False)

import pybedtools
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import pingouin as pg
import pandas as pd
import numpy as np
import re
import sys

def main():

	asb_input = sys.argv[1]
	qbic_input = sys.argv[2]
	
	asb  = pd.read_csv(asb_input, header=0,index_col=0)
	qbic = pd.read_csv(qbic_input, header=0)

	df = pd.merge(qbic, asb, left_on = ['gene','allele'], 
			right_on = ['gene','allele'])

	df['consensus'] = (df['occu'] == df['binding_count'])
	print(df[['gene','allele','occu','binding_count']])
	print(df['consensus'].sum()/df['consensus'].count())


if __name__ == "__main__":
	main()



