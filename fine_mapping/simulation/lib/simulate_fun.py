import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import numpy as np
import random
import scipy as sp
import sys

class eqtlAssociation:

	def __init__(self, genotype, expression):
		self.genotype = genotype
		self.expression = expression

	def ProcessInput(self):
		if self.genotype.shape[0] != self.expression.shape[0]:
			print('non-equivalent size of input matrices')
			return 0

	def sgl_pair_lr(self, expression_vec, genotype_vec):
		print(expression_vec)
		print(genotype_vec)
		slope, intercept, r, pval, se = stats.linregress(expression_vec,genotype_vec)
		return slope, pval

	def show_lr_results(self,thres=1e-5):
		select = self.lr_res['pval'] < thres
		print(self.lr_res.loc[select])

	def AssocationTest(self):
		self.ProcessInput()
		exp_vec = self.expression.values.reshape(-1)
		#calculate effect size from linear regression
		lr_res = self.genotype.apply(lambda row: self.sgl_pair_lr(exp_vec,row), axis=1)
		lr_res_df = pd.DataFrame.from_records(lr_res.values,columns=['lr','pval'], index=lr_res.index)
		self.lr_res = lr_res_df
		self.show_lr_results(1)

def simulate_genotype(maf, n):

	maf = float(maf)
	n   = int(n)

	# number of aa
	p = 1 - maf
	q = maf
	n_aa = int(n*q*q)
	n_Aa = int(n*2*q*p)
	n_AA = n - n_aa - n_Aa
	sample_idx = range(n)

	# initial genotype profile
	genotype_profile = pd.DataFrame(index=sample_idx, columns=['snp1'])	
	genotype_profile.loc[0:n_aa-1,         'snp1'] = '1|1'
	genotype_profile.loc[n_aa:n_aa+n_Aa-1, 'snp1'] = np.random.choice(['0|1','1|0'], size=n_Aa)
	genotype_profile.loc[n_aa+n_Aa:n,      'snp1'] = '0|0'

	return genotype_profile.sample(frac=1)

def rand_sample_hist(vals):

	values,indices = np.histogram(vals,bins=100)
	values=values.astype(np.float32)
	weights=values/np.sum(values)

	new_random = np.random.choice(indices[1:],1,p=weights)
	return new_random

def sample_PB(m_kon, m_koff, m_syn, n=1000):

	m_array_p = np.ones(n)
	m_array_x = np.ones(n)

	m_array_p = np.random.beta(m_kon,m_koff,n)
	m_array_x = np.random.poisson(m_syn*m_array_p)
	return m_array_x
	
class simulate_eQTL:

#	def __init__(self,kinetics_input):
#		kinetics = pd.read_csv(kinetics_input, header=0,index_col=0)
#		kinetics['bs'] = kinetics['ksyn']/kinetics['koff']
#		kinetics['bf'] = kinetics['kon'] * kinetics['koff'] / (kinetics['kon'] + kinetics['koff'])
#		log_kp = list(map(lambda x: 'log10' + x, kinetics.columns))
#
#		for kp in ['kon','koff', 'ksyn', 'bs', 'bf']:
#			kinetics['log10'+kp] = kinetics[kp].transform(np.log10)
#		self.kinetics = kinetics

	def simulate_kp(self, H1_bkp=None, H2_bkp=None):

#		if not H1_bkp or H2_bkp:
#			sample_log10bs = rand_sample_hist(self.kinetics['log10bs'].values)[0] 
#			sample_log10bf = rand_sample_hist(self.kinetics['log10bf'].values)[0] 
#			sample_log10r = rand_sample_hist(self.kinetics['log10ksyn'].values)[0] 
#
#			sample_r    = np.power(sample_log10r, 10)
#			sample_bs   = np.power(sample_log10bs, 10)
#			sample_bf   = np.power(sample_log10bf, 10)
		
		sample_bs,sample_bf,sample_r = H1_bkp
		sample_bs2,sample_bf2,sample_r2 = H2_bkp

		sample_koff = sample_r/sample_bs
		sample_kon  = 1.0/(1.0/sample_bf - 1.0/sample_koff)		

		H1_kp  = [sample_kon,sample_koff,sample_r,sample_bs,sample_bf]

		sample_koff2 = sample_r2/sample_bs2
		sample_kon2  = 1.0/(1.0/sample_bf2 - 1.0/sample_koff2)		

		H2_kp = [sample_kon2,sample_koff2,sample_r2,sample_bs2,sample_bf2]
	
		self.H1_kp = H1_kp
		self.H2_kp = H2_kp

		#self.plot(H1_kp, H2_kp)

		print("0:", self.H1_kp)
		print("1:", self.H2_kp)

	def plot(self,kp1, kp2):

		sample_kon,sample_koff,sample_r,sample_bs,sample_bf = kp1
		sample_kon2,sample_koff2,sample_r2,sample_bs2,sample_bf2 = kp2

		kpe_exp1 = sample_PB(sample_kon,sample_koff,sample_r, 2000)
		gene_density1 = stats.gaussian_kde(kpe_exp1)	
		x_density1 = np.arange(-0.1,max(kpe_exp1),0.5)

		kpe_exp2 = sample_PB(sample_kon2,sample_koff2,sample_r2, 2000)
		gene_density2 = stats.gaussian_kde(kpe_exp2)	
		x_density2 = np.arange(-0.1,max(kpe_exp2),0.5)

		fig,ax = plt.subplots(constrained_layout=True,figsize=(3,3))
		ax.plot(x_density1, gene_density1(x_density1), color='#1f77b4', ls='-', lw=2)
		ax.hist(kpe_exp1, color='#1f77b4',density=True, alpha=0.5)
		ax.plot(x_density2, gene_density2(x_density2), color='#ff7f0e', ls='-', lw=2)
		ax.hist(kpe_exp2, color='#ff7f0e', density=True, alpha=0.5)
		plt.savefig('ase_distribution.png')

	def simulate_expression(self, m, genotype):

		sample_kon,sample_koff,sample_r,sample_bs,sample_bf = self.H1_kp
		sample_kon2,sample_koff2,sample_r2,sample_bs2,sample_bf2 = self.H2_kp

		if genotype == '1|1':
			H1_exp = sample_PB(sample_kon2,sample_koff2,sample_r2, m)
			H2_exp = sample_PB(sample_kon2,sample_koff2,sample_r2, m)
			mean_exp = np.mean(H1_exp) + np.mean(H2_exp)

		elif genotype == '0|1':
			H1_exp = sample_PB(sample_kon,sample_koff,sample_r, m)
			H2_exp = sample_PB(sample_kon2,sample_koff2,sample_r2, m)
			mean_exp = np.mean(H1_exp) + np.mean(H2_exp)
			print(genotype,"H1_mean:", np.mean(H1_exp),"H1_var:",np.var(H1_exp),
				"H2_mean:",np.mean(H2_exp), "H2_var:",np.var(H2_exp))

		elif genotype == '1|0':
			H1_exp = sample_PB(sample_kon,sample_koff,sample_r, m)
			H2_exp = sample_PB(sample_kon2,sample_koff2,sample_r2, m)
			mean_exp = np.mean(H1_exp) + np.mean(H2_exp)
			print(genotype, "H1_mean:", np.mean(H2_exp),"H1_var:",np.var(H2_exp),
				"H2_mean:",np.mean(H1_exp), "H2_var:",np.var(H1_exp))

		elif genotype == '0|0':
			H1_exp = sample_PB(sample_kon,sample_koff,sample_r, m)
			H2_exp = sample_PB(sample_kon,sample_koff,sample_r, m)
			mean_exp = np.mean(H1_exp) + np.mean(H2_exp)
		return np.vstack([H1_exp, H2_exp])

if __name__ == '__main__':
	main()



