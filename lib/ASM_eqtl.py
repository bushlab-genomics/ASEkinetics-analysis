import pandas as pd
import numpy as np
import re
import sys
	
def main():

	asm_input  = sys.argv[1]
	eqtl_input = sys.argv[2]
	label      = sys.argv[4]

	asm = pd.read_csv(asm_input, header=0, sep='\t')
	asm['chr'] = asm['Chromosome'].str.strip('chr')
	eqtl = pd.read_csv(eqtl_input, header=0)

	asm['loc'] = asm[['chr','SNP_postion_ends']].astype(str).agg(':'.join, axis=1)
	eqtl[['chr','pos','ref','alt']] = eqtl['ID'].str.split(':', expand=True)
	eqtl['loc'] = eqtl[['chr','pos']].astype(str).agg(':'.join, axis=1)

	asm_eqtl = pd.merge(asm, eqtl, on='loc',suffixes=('','_eqtl'))

	select0 = (abs(asm_eqtl['log(burst frequency_fc)']) > np.log10(2)) | (abs(asm_eqtl['log(burst size_fc)']) > np.log10(2)) 
	select00 = (abs(asm_eqtl['log(mean_fc)']) > np.log10(2)) | (abs(asm_eqtl['log(var_fc)']) > np.log10(2))
	select00 = True

	select1 = ( asm_eqtl['log(burst frequency_fc)']  < 1 ) & ( asm_eqtl['PvsM_es'] < 0 ) 
	select2 = ( asm_eqtl['log(burst frequency_fc)']  > 1 ) & ( asm_eqtl['PvsM_es'] > 0 )
	select3 = ( asm_eqtl['log(burst size_fc)']  < 1 ) & ( asm_eqtl['PvsM_es'] < 0 ) 
	select4 = ( asm_eqtl['log(burst size_fc)']  > 1 ) & ( asm_eqtl['PvsM_es'] > 0 )

	select = select0 & select00 & (select1 | select2 | select3 | select4) # & ( asm_eqtl['pvalue'] < 1e-8)
	asm_eqtl_select = asm_eqtl.loc[select].drop_duplicates()
	asm_eqtl_select.to_csv('%s_eqtl_asm.csv'%label,index=False)

	GC_input = sys.argv[3]
	GC = pd.read_csv(GC_input, header=0)

	asm_eqtl_GC_select = pd.merge(asm_eqtl_select,GC[['ID','genotype']],left_on='rsid', right_on='ID')
	asm_eqtl_GC_select.to_csv('%s_eqtl_asm_GC.csv'%label,index=False)

if __name__ == "__main__":

	main()

