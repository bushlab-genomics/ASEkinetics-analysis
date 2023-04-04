#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(LDlinkR)

file=args[1]
df=read.table(file, header=FALSE, sep=',',col.names=c('asb_rsid','eqtl_rsid'))
print(df)
for(i in 1:nrow(df)){
	rsid1=df[i,"asb_rsid"]
	rsid2=df[i,"eqtl_rsid"]
	print(rsid1)
	print(rsid2)
	res = LDpair(var1 = "rs1158528361", var2 = "rs11187580", pop = "CEU", token = '92e8f3bc25c6')
	print(res)
}
