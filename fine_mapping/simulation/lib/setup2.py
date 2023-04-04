from pb_est import *
from simulate_fun import *
from evaluation import *
import sys

def compare(kp_df):

	kp_cols = ['kon','koff','ksyn','bs','bf']
	bool_select = (kp_df['val'] == '0')
	res = pd.DataFrame(columns=['kp','pval','effect size'])
	for i, kp in enumerate(kp_cols):
		h1_val = kp_df.loc[bool_select, kp]
		h2_val = kp_df.loc[~bool_select, kp]
		statistic, pval = sp.stats.ttest_ind(h1_val,h2_val)
		res.loc[i, 'kp' ] = kp
		res.loc[i, 'effect size' ] = statistic
		res.loc[i, 'pval' ] = pval

	return res

def evaluate(org_kp, org_data, method='pb'):

	n = len(org_data)
	sim_data = pb_simulation(org_kp,n)

	obj  = mRNAkinetics(sim_data)
	obj.MaximumLikelihood()
	sim_kp = obj.get_estimate()

	chisq, p = GoF(sim_data, org_data)
	ks_p = kstest(sim_data, org_data)
	simlr = simLikelihoodRatioTest(sim_kp, sim_data, org_kp, org_data)	
	return p, simlr, ks_p

def main():
	maf = sys.argv[1]
	sample = int(sys.argv[2])
	n_cells = int(sys.argv[3])

	# input kinetics
	kinetics_input = sys.argv[4]	
	kinetics = np.loadtxt(kinetics_input)
	sample_bs, sample_bf, sample_r = kinetics

	out_file = sys.argv[5]

	# comparable kinetics
	sample_bs2 = sample_bs * 2 
	sample_bf2 = sample_bf 
	sample_r2  = sample_r
	H1_bkp = [sample_bs, sample_bf, sample_r]
	H2_bkp = [sample_bs2,sample_bf2,sample_r2]

	genotype = simulate_genotype(maf, sample)
	sim_obj = simulate_eQTL()
	sim_obj.simulate_kp(H1_bkp,H2_bkp)

	expression_df = pd.DataFrame()
	for i, row in genotype.iterrows():
		temp_exp = sim_obj.simulate_expression(n_cells, row['snp1'])
		temp_df = pd.DataFrame(data=temp_exp, columns=list(range(n_cells)), 
					index=['%s-H1'%str(i),'%s-H2'%str(i)])

		if len(expression_df) < 1: expression_df = temp_df
		else: expression_df = pd.concat([expression_df, temp_df])

	genotype_count = genotype.apply(lambda row: 1 if (row['snp1']=='1|0') or (row['snp1']=='0|1') else \
					( 2 if row['snp1']=='1|1' else 0), axis=1)
	genotype_allele_df = genotype.copy()
	genotype_allele_df[['H1','H2']] = genotype_allele_df['snp1'].str.split('|', n=1, expand=True)
	genotype_allele_df.reset_index(inplace=True)
	genotype_allele_df = pd.melt(genotype_allele_df, value_vars=['H1','H2'], var_name='allele',\
					value_name='val', id_vars='index')
	genotype_allele_df['uniq_idx'] = genotype_allele_df[['index','allele']].astype(str).agg('-'.join,axis=1)

	bulk_exp = expression_df.apply(lambda row: np.mean(row), axis=1)
	bulk_exp_df = bulk_exp.to_frame('exp')
	bulk_exp_df[['gene_idx','allele']] = bulk_exp_df.index.to_series().str.split('-', n=1, expand=True)
	
	bulk_exp_final = bulk_exp_df.groupby('gene_idx')['exp'].mean()

	genotype_count.sort_index(inplace=True)
	bulk_exp_final.sort_index(inplace=True)

	# test eqtl association
	test_obj = eqtlAssociation(genotype_count, bulk_exp_final)
	slope, pval = test_obj.sgl_pair_lr(genotype_count.values, bulk_exp_final.values.reshape(-1))
		
	kp_df = pd.DataFrame(columns=['kon','koff','ksyn','bs','bf'])
	# estimate kineicts
	for i, row in expression_df.iterrows():
		obj1  = mRNAkinetics(row.values)
		obj1.MaximumLikelihood()
		kpe = obj1.get_estimate()
		kon, koff, ksyn = kpe
		bf = (kon * koff)/(kon+koff)
		bs = ksyn/koff
		record = [ kon, koff, ksyn, bs, bf ]
		kp_df.loc[i,:] = record
	kp_df = pd.merge(kp_df, genotype_allele_df[['uniq_idx','val']],left_index=True,right_on='uniq_idx')
	res_df = compare(kp_df)

	res_df.loc[len(res_df)]  = ['eqtl',pval, slope]
	res_df['n_sample'] = sample
	res_df['n_cell']   = n_cells
	res_df.to_csv(out_file,index=False,header=False)

if __name__ == "__main__":
	main()






