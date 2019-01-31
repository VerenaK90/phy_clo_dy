source Phylogenetic_inference.R

load example data: Example.data.RData

The example data is a list containing information of mutated loci in a pair of a primary and recurrent tumor from a male person. The readcounts are provided in the list elements “read.count.matrix.prim” and “read.count.matrix.rec”, coverage and B-allele frequencies are given in “Coverage.B.alleles.prim” and “Coverage.B.alleles.rec” and the chromosomal coordinates of the SNVs/indels in “Locus.prim” and “Locus.rec”

Run the following code to fit up to three subclones per tumor sample (runs several hours):

fitting.result <- model.fits(read.count.matrix.prim = example.data$read.count.matrix.prim, read.count.matrix.rec = example.data$read.count.matrix.rec,
                        coverage.ratio.baf.prim = example.data$Coverage.B.alleles.prim, coverage.ratio.baf.rec = example.data$Coverage.B.alleles.rec,
                        locus = example.data$Locus.prim, sex = "male", directory = "~/")


You obtain a list of fitting results. Trees of increasing complexity were fit to the datasets. Each element of the list corresponds to a different candidate tree. The candidate tree is stored in the “tree” element of each sublist. Each row of the “Parameters.prim” and “Parameters.rec” elements contain the subclone sizes inferred from multiple fits (one fit per row). The last column corresponds to normal tissue. The likelihood of each fit is stored in the “Likelihood” element. Finally, the “which.prim” and “which.rec” elements encrypt the columns of the tree matrix belonging to the primary and relapse subclones respectively. 

E.g., if you want to have a look at the best result of candidate tree 15, run the following

tree <- fitting.result[[15]]$tree (prints the binary tree matrix of the candidate tree)
colnames(tree) <- rep(“”, ncol(tree))
# name the tree tips according to the tumor sample they belong to:
colnames(tree)[fitting.result[[15]]$which.prim] <- rep("Primary", length(fitting.result[[15]]$which.prim))
colnames(tree)[fitting.result[[15]]$which.rec] <- rep("Relapse", length(fitting.result[[15]]$which.rec))
tree[,c(fitting.result[[15]]$which.prim, fitting.result[[15]]$which.rec)]

# obtain the fit with the highest likelihood for this tree: 

best.tree <- which.min(fitting.result[[15]]$Likelihood)
# parameters of primary tumor:
fitting.result[[15]]$Parameters.prim[best.tree,]
# parameters of recurrent tumor: 
fitting.result[[15]]$Parameters.rec[best.tree,]

To sort individual mutations on this specific tree, we have to run the expectation function once again:

sort.mutations <- Expectation(tree = fitting.result[[15]]$tree, read.count.matrix.prim = example.data$read.count.matrix.prim, read.count.matrix.rec = example.data$read.count.matrix.rec, parms.prim = fitting.result[[15]]$Parameters.prim[1,], parms.rec = fitting.result[[15]]$Parameters.rec[1,], coverage.ratio.baf.prim = example.data$Coverage.B.alleles.prim, coverage.ratio.baf.rec = example.data$Coverage.B.alleles.rec, locus = example.data$Locus.prim, sex = "male", which.prim = fitting.result[[15]]$which.prim, which.rec = fitting.result[[15]]$which.rec, likelihood.space = F)

The value returned by Expectation() contains several matrices indicating the fraction of alleles mutated per subclone (“rel.mutated.alleles.prim/rec”), the copy number state (“copy.number.prim/rec”) and the readcounts at maximum likelihood (“reacounts.prim/rec”). The columns are in the order of the parameter estimates.


Select a tree:

To infer the most likely tree, we employ a modified Bayesian Information Criterion and require that at least 50% of the clonal mutations map unambiguously to the trunc of the tree. This can be tuned with the parameter min.non.ambiguous.clonals. Run the following code to store for each different candidate tree the index of the best fit along with the log-likelihood and the modified BIC value. In addition we store the mean squared error of the variant allele frequencies (measured vs. predicted) of clonal and all mutations along with the percentage of clonal mutations that are ambiguously mapped on the tree

best.sol.per.cat <- Best.solutions(fitting.result)

To select a solution, different criteria can be applied. We required a good fit of clonal mutations along with a good BIC. Note that there is not a “best” solution, but several solutions can fit the data well. 


best.sol.per.cat. <- best.sol.per.cat[best.sol.per.cat$`mean sq error clonal` <= quantile(p=0.1,best.sol.per.cat$`mean sq error clonal`, na.rm = T) ,]
cand.models <- as.numeric(rownames(best.sol.per.cat.)[which(best.sol.per.cat.[,"delta.mBIC"] <= min(best.sol.per.cat.[,"delta.mBIC"][is.finite(best.sol.per.cat.[,"delta.mBIC"])]) + 10 & best.sol.per.cat.[,"delta.mBIC"]!=0 &
                       is.finite(best.sol.per.cat.[,"delta.mBIC"]))])
print(cand.models)


Detect clonal mutations:

First, in order to sort the mutations on the selected tree we have to run the expectation step again
which.solution <- best.sol.per.cat[cand.models,"Index"]
sort.mutations <- Expectation(tree = fitting.result[[cand.models]]$tree, read.count.matrix.prim = example.data$read.count.matrix.prim, read.count.matrix.rec = example.data$read.count.matrix.rec, parms.prim = fitting.result[[cand.models]]$Parameters.prim[which.solution,], parms.rec = fitting.result[[cand.models]]$Parameters.rec[which.solution,], coverage.ratio.baf.prim = example.data$Coverage.B.alleles.prim, coverage.ratio.baf.rec = example.data$Coverage.B.alleles.rec, locus = example.data$Locus.prim, sex = "male", which.prim = fitting.result[[cand.models]]$which.prim, which.rec = fitting.result[[cand.models]]$which.rec, likelihood.space = F)
parms.prim <- fitting.result[[cand.models]]$Parameters.prim[which.solution,]
parms.rec <- fitting.result[[cand.models]]$Parameters.rec[which.solution,]


## Then we can detect clonal mutations. To do this, we look at the fitting.result $rel.mutated.alleles.prim and require a mutated allele in all subclones:
clonal.mutations.prim <- rownames(sort.mutations$rel.mutated.alleles.prim)[rowSums((sort.mutations$rel.mutated.alleles.prim[,seq(2,ncol(sort.mutations$rel.mutated.alleles.prim)-2,2),drop=F]!=0)*1)==(length(parms.prim)-1)]
clonal.mutations.rec <- rownames(sort.mutations$rel.mutated.alleles.rec)[rowSums((sort.mutations$rel.mutated.alleles.rec[,seq(2,ncol(sort.mutations$rel.mutated.alleles.rec)-2,2),drop=F]!=0)*1)==(length(parms.rec)-1)]
clonal.mutations <- clonal.mutations.prim[clonal.mutations.prim %in% clonal.mutations.rec]
head(clonal.mutations)

Analyze copy number changes:

For example, we might want to know how the copy numbers at chromosome 10 look like:

chr10.loci<- rownames(example.data$Locus.prim)[example.data$Locus.prim[,1]==10]
head(sort.mutations$copy.number.prim[intersect(chr10.loci, rownames(sort.mutations$copy.number.prim)),])
head(sort.mutations$copy.number.rec[intersect(chr10.loci, rownames(sort.mutations$copy.number.rec)),])


Finally, if we now want to test an additional tree, we can extend the tree by additional subclones. This is achieved by providing the best tree and extending it. All positions in the tree are tested to add this subclone. E.g., we want to extend tree # 96:


tree.to.extend <- fitting.result[[96]]$tree
which.prim <- fitting.result[[96]]$which.prim
which.rec <- fitting.result[[96]]$which.rec

extended.tree <- best.model(read.count.matrix.prim = example.data$read.count.matrix.prim, read.count.matrix.rec = example.data$read.count.matrix.rec, coverage.ratio.baf.prim = example.data$Coverage.B.alleles.prim, coverage.ratio.baf.rec = example.data$Coverage.B.alleles.rec, locus = example.data$Locus.prim,
                       sex = "male", directory = "~/", tree.best.fit = tree.to.extend, which.prim.best.fit = which.prim, which.rec.best.fit = which.rec, add.prim.clone = T, add.rec.clone = F)