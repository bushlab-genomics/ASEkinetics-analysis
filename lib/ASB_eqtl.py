import pandas as pd
import numpy as np
import re
import sys
	
def main():

	asb_input  = sys.argv[1]
	eqtl_input = sys.argv[2]
	label      = sys.argv[4]

	asb_raw = pd.read_csv(asb_input, header=0)
	asb = asb_raw[asb_raw['peak.isASB']==True]
	asb['chr'] = asb['CHROM'].str.strip('chr')
	eqtl = pd.read_csv(eqtl_input, header=0)

	asb['loc'] = asb[['chr','POS']].astype(str).agg(':'.join, axis=1)
	eqtl[['chr','pos','ref','alt']] = eqtl['ID'].str.split(':', expand=True)
	eqtl['loc'] = eqtl[['chr','pos']].astype(str).agg(':'.join, axis=1)

	asb_eqtl = pd.merge(asb,eqtl, on='loc',suffixes=('','_eqtl'))
	select0 = (abs(asb_eqtl['log(burst frequency_fc)']) > np.log10(2)) | (abs(asb_eqtl['log(burst size_fc)']) > np.log10(2)) 
	select00 = (abs(asb_eqtl['log(mean_fc)']) > np.log10(2)) | (abs(asb_eqtl['log(var_fc)']) > np.log10(2))

	select1 = ( asb_eqtl['log(burst frequency_fc)']  < 1 ) & ( asb_eqtl['PvsM_es'] < 0 ) 
	select2 = ( asb_eqtl['log(burst frequency_fc)']  > 1 ) & ( asb_eqtl['PvsM_es'] > 0 )
	select3 = ( asb_eqtl['log(burst size_fc)']  < 1 ) & ( asb_eqtl['PvsM_es'] < 0 ) 
	select4 = ( asb_eqtl['log(burst size_fc)']  > 1 ) & ( asb_eqtl['PvsM_es'] > 0 )

	select = select0 & select00 & (select1 | select2 | select3 | select4)  & ( asb_eqtl['pvalue'] < 1e-8)
	asb_eqtl_select = asb_eqtl.loc[select].drop_duplicates()
	asb_eqtl_select.drop_duplicates().to_csv('%s_asb_eqtl.csv'%label,index=False)

	GC_input = sys.argv[3]
	GC = pd.read_csv(GC_input, header=0)

	asb_GC = pd.merge(asb,GC,left_on='ID', right_on='loc',suffixes=('','_GC'))
	asb_GC.drop_duplicates().to_csv('%s_asb_GC.csv'%label,index=False)

	qbic_input = sys.argv[5]
	qbic = pd.read_csv(qbic_input, header=0)	
	asb_eqtl_qbic = pd.merge(asb_eqtl, qbic[['TF_gene','binding_status','ID_pos']],
				left_on = 'ID', right_on='ID_pos')
	asb_eqtl_qbic_temp= \
		asb_eqtl_qbic[['genotype','binding_status','asb_allele']].astype(str).agg('_'.join, axis=1)

	select = (asb_eqtl_qbic_temp == '0|1_unbound>bound_H2') | \
		(asb_eqtl_qbic_temp == '0|1_bound>unbound_H1') | \
		(asb_eqtl_qbic_temp == '1|0_unbound>bound_H1') | \
		(asb_eqtl_qbic_temp == '1|0_bound>unbound_H2') 	
	asb_eqtl_qbic = asb_eqtl_qbic.loc[select].drop_duplicates()
	asb_eqtl_qbic.to_csv('%s_asb_eqtl_qbic.csv'%label,index=False)


	asb_GC_qbic = pd.merge(asb_GC, qbic[['TF_gene','binding_status','ID_pos']],
				left_on = 'ID' , right_on='ID_pos')
	asb_GC_qbic_temp = \
		asb_GC_qbic[['genotype','binding_status','asb_allele']].astype(str).agg('_'.join, axis=1)

	select = (asb_GC_qbic_temp == '0|1_unbound>bound_H2') | \
		(asb_GC_qbic_temp == '0|1_bound>unbound_H1') | \
		(asb_GC_qbic_temp == '1|0_unbound>bound_H1') | \
		(asb_GC_qbic_temp == '1|0_bound>unbound_H2') 
	asb_GC_qbic = asb_GC_qbic.loc[select].drop_duplicates()
	asb_GC_qbic.to_csv('%s_asb_GC_qbic.csv'%label,index=False)


	asb_qbic = pd.merge(asb,qbic[['TF_gene','binding_status','ID_pos']],\
				left_on='ID', right_on='ID_pos')
	asb_qbic_temp = \
		asb_qbic[['genotype','binding_status','asb_allele']].astype(str).agg('_'.join, axis=1)
	
	select = (asb_qbic_temp == '0|1_unbound>bound_H2') | \
		(asb_qbic_temp == '0|1_bound>unbound_H1') | \
		(asb_qbic_temp == '1|0_unbound>bound_H1') | \
		(asb_qbic_temp == '1|0_bound>unbound_H2')
	asb_qbic = asb_qbic.loc[select].drop_duplicates()
	asb_qbic.to_csv('%s_asb_qbic.csv'%label,index=False)

	yri_gt_input = sys.argv[6]
	yri_gt = pd.read_csv(yri_gt_input, sep='\t', usecols=list(range(7))+[9], \
		names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL','FILTER','genotype'])
	shared_eqtl_qbic = pd.merge(asb_eqtl_qbic, yri_gt[['ID','genotype']], on=['ID'], suffixes=('','_yri'))
	shared_GC_qbic   = pd.merge(asb_GC_qbic,   yri_gt[['ID','genotype']], on=['ID'], suffixes=('','_yri'))
	shared_qbic      = pd.merge(asb_qbic,      yri_gt[['ID','genotype']], on=['ID'], suffixes=('','_yri'))

	shared_eqtl_qbic.to_csv('shared_asb_eqtl_qbic.csv',index=False)
	shared_GC_qbic.to_csv('shared_asb_GC_qbic.csv',index=False)
	shared_qbic.to_csv('shared_asb_qbic.csv',index=False)

if __name__ == "__main__":

	main()


