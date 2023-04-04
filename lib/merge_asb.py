from scipy import stats
import pingouin as pg
import pandas as pd
from pybedtools import BedTool
import numpy as np
import re
import sys

def pref(x):

        if   (x['Corrected.AR'] > 0.5) & (x['genotype'] == '1|0'):
                return  'H2'
        elif (x['Corrected.AR'] > 0.5) & (x['genotype'] == '0|1'):
                return  'H1'
        elif   (x['Corrected.AR'] < 0.5) & (x['genotype'] == '1|0'):
                return  'H1'
        elif (x['Corrected.AR'] < 0.5) & (x['genotype'] == '0|1'):
                return  'H2'
# keep the ASB consistent with phase infor
def filterASB(asb_input, genotype_input, output):

	asb = pd.read_csv(asb_input, header=0,index_col=0)
	genotype = pd.read_csv(genotype_input, sep='\t', comment="#", header=None, usecols=list(range(5)) + [9],
		names=['CHROM','POS','ID','REF','ALT','genotype'],dtype={'CHROM':'str'}, compression='gzip')
	
	genotype['ID'] = genotype[['CHROM','POS','REF','ALT']].astype(str).agg(':'.join,axis=1)
	asb_genotype = pd.merge(asb, genotype[['ID','genotype']], on='ID')
	asb_genotype['asb_allele'] = asb_genotype.apply(lambda row: pref(row), axis=1)	
	asb_genotype_varisASB = asb_genotype.loc[asb_genotype['Allvar.isASB']==True]	
	asb_genotype_peakisASB = asb_genotype_varisASB.groupby('peak_name')['asb_allele'].nunique()
	peakisASB = asb_genotype_peakisASB.loc[asb_genotype_peakisASB <=1].index

	asb_genotype['peak.isASB']=False
	asb_genotype.loc[asb_genotype['peak_name'].isin(peakisASB),'peak.isASB'] = True	
	asb_genotype.to_csv(output+'_isASB.final',index=False)

#keep ASB consistent within the peak
def mergePeak(asb_input,het_input,bed_input,subject,output):

	asb = pd.read_csv(asb_input, header=0,index_col=0)
	het = BedTool(het_input)
	het_cols = ['chr','start','end','ID','score','strand']
	peaks = BedTool(bed_input)
	peaks_cols = ['p_chr','p_start','p_end']

	het_peaks = het.intersect(peaks, wa=True, wb=True)
	het_peaks_df = pd.read_table(het_peaks.fn, header=None, \
			usecols=range(9),names=het_cols + peaks_cols)

	het_peaks_df['peak_name'] = het_peaks_df[peaks_cols].astype(str).agg(':'.join,axis=1)

	asb_peaks = pd.merge(het_peaks_df, asb, left_on='ID', right_on='%s.ID'%subject)
	asb_peaks_grp = asb_peaks.groupby('peak_name')['%s.isASB'%subject].mean()
	#select peak that only contains SNPs showing ASB
	ASBpeaks = asb_peaks_grp.loc[asb_peaks_grp ==1].index

	asb_peaks['Allvar.isASB']=False
	asb_peaks.loc[asb_peaks['peak_name'].isin(ASBpeaks), 'Allvar.isASB'] = True
	asb_peaks[asb.columns.tolist() + ['peak_name','Allvar.isASB']].to_csv(output+'.merge',index=False)
	
if __name__ == "__main__":

	args = sys.argv
	globals()[args[1]](*args[2:])


