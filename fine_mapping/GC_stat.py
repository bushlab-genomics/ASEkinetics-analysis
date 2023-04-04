import pandas as pd
import re
import sys

def main():

	GC_input = sys.argv[1]
	vcf_input  = sys.argv[2]
	output = sys.argv[3]
	
	GC = pd.read_csv(GC_input, header=None, sep='\t',usecols=list(range(5)), names=['chr','pos','ID','ref','alt'])
	vcf = pd.read_csv(vcf_input, header=None, sep='\t', comment='#' ,usecols=list(range(6)) + [9],\
			names=['chr','pos','ID','ref','alt','qual','genotype'],compression='gzip', dtype={'chr':str})

	GC['loc'] = GC[['chr','pos','ref','alt']].astype(str).agg(':'.join, axis=1)
	
	vcf['loc'] = vcf[['chr','pos','ref','alt']].astype(str).agg(':'.join, axis=1)
	GC_vcf = pd.merge(GC, vcf[['loc','genotype']], on='loc')

	GC_vcf.to_csv('%s_GC.tsv'%output, index=False)

if __name__ == "__main__":
	main()



