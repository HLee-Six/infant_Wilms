# 2022.10.15
# wilms organoids
# phylogeny construction
# PD37590 (normal colon case)

# for an explanation of the steps, see word document entitled 'Methods'

library('GenomicRanges')
library('BSgenome.Hsapiens.UCSC.hg19')
library("dplyr")
library("VGAM")
library('ape')
library('treemut')

setwd('/Users/hl11/Documents/Research/Projects/wilms_organoids/Drafts/202212_submission_version/github_upload/')

# read in the two pileups.
aa10 <- read.csv('PD37590_mq10bq10.tsv', sep='\t', header=T, stringsAsFactors = F)
aa30 <- read.csv('PD37590_mq30bq30.tsv', sep='\t', header=T, stringsAsFactors = F)
aa10$mutkey <- paste0(aa10$Chrom, '_', aa10$Pos, '_', aa10$Ref, '_', aa10$Alt)
aa30$mutkey <- paste0(aa30$Chrom, '_', aa30$Pos, '_', aa30$Ref, '_', aa30$Alt)

# read in info about the samples
# I am going to build a phylogeny just for the crypts from the right colon.
# there are no normal samples at all for this patient.
crypts <- read.csv('colon_paper_supplementary_table2.csv', header=T, stringsAsFactors = F)
rcrypts <- crypts[crypts$patient=='PD37590' & crypts$site=='Right','crypt']
tcrypts <- crypts[crypts$patient=='PD37590' & crypts$site!='Right','crypt']

allsamps <- gsub('_VAF', '', grep("VAF", colnames(aa30), value=T))
samps <- rcrypts[rcrypts %in% allsamps]
normals <- tcrypts[tcrypts %in% allsamps]

# get useful vectors of column names to work on the pileup.
mtrs <- paste0(samps, '_MTR')
deps <- paste0(samps, '_DEP')
vafs <- paste0(samps, '_VAF')
ofs <- paste0(samps, '_OFS')

## separate mutations into two groups based on their vaf in the bulk normal samples.
# mutations that are 'clonal-in-bulk-normal' include germline + early embryonic mutations.
# mutations that are 'subclonal-in-bulk-normal' include early embryonic + somatic.

# calculate the vaf of every mutation across all normal samples.
filtmat <- aa30[,c('mutkey', paste0(rep(normals, each=2), c('_MTR', '_DEP')))]
filtmat$norm_mtr <- rowSums(filtmat[grep('MTR', colnames(filtmat))])
filtmat$norm_dep <- rowSums(filtmat[grep('DEP', colnames(filtmat))])
filtmat$norm_vaf <- filtmat$norm_mtr/filtmat$norm_dep
# there are too many samples for filtmat to be a useful thing.
filtmat <- filtmat[,grep('PD37590',colnames(filtmat), invert = T)]

# only keep mutations that are well covered in the bulk normals. If they are poorly covered, they could be dodgy sites, or we can't confidently exclude them from being germline.
# set a cutoff for good depth for the normal samples. This is done manually based on the histogram for this sample.
hist(filtmat$norm_dep, 100)
hist(filtmat$norm_dep, 100, xlim=c(0,500))  
abline(v=200, col='red')
abline(v=400, col='red')
norm_cutoff_min <- 200
norm_cutoff_max <- 400

# choose a vaf cutoff to separate [clonal-in-bulk-normal] from [subclonal-in-bulk-normal]
# do this by simulating clonal mutations, and see how low their vaf goes. 
# simulate mtrs from the observed depths using the binomial distribution, and use this to calculate the vaf of simulated mutations.
# from each of these simulations, find the lowest vaf achieved.
# run the simulation 1,000 times. 
# even a small number of germline mutations getting put in the subclonal category could be quite damaging for the structure at the top of the tree.
# so choose a cutoff which means that in 99% of simulations no germline mutations would have been called as subclonal_in_the_bulk.

bulk_norm_deps <- filtmat$norm_dep[filtmat$norm_dep>norm_cutoff_min & filtmat$norm_dep<norm_cutoff_max]
sim_vafs_mins <- c()
for (i in 1:1000) {
  if ((i %% 100)==0) print(i)
  these_mtrs <- sapply(bulk_norm_deps, function(this_depth) rbinom(n=1, size=this_depth, prob=0.5))
  these_vafs <- these_mtrs/bulk_norm_deps
  sim_vafs_mins <- c(sim_vafs_mins, min(these_vafs))
}
hist(sim_vafs_mins, xlim=c(0,1), 50)
vafcutoff <- sort(sim_vafs_mins, decreasing = F)[round(length(sim_vafs_mins)*0.01)]

# now separate our mutations into those that are clonal in the bulk or not
filtmat$clonal_in_normals <- 'yes'
filtmat$clonal_in_normals[filtmat$norm_vaf<vafcutoff & filtmat$norm_dep>norm_cutoff_min & filtmat$norm_dep<norm_cutoff_max] <- 'no'
# any that were poorly covered in the normals it is difficult to comment on, and I don't want to use them for tree building.
filtmat$clonal_in_normals[filtmat$norm_dep<=norm_cutoff_min | filtmat$norm_dep >= norm_cutoff_max] <- 'poor_coverage_in_normal'
table(filtmat$clonal_in_normals) 

#####
# now take the ones that were subclonal in normals
# use these to build the phylogeny - but first, need to apply some filters.

## apply Tim's betabinomial filter to each pileup.
estimateRho_gridml = function(NV_vec,NR_vec) {
  #Make sure depth is non-zero
  NV_vec=NV_vec[NR_vec>0]
  NR_vec=NR_vec[NR_vec>0]
  
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rhovec = 10^seq(-6,-0.05,by=0.05) # rho will be bounded within 1e-6 and 0.89
  mu=sum(NV_vec)/sum(NR_vec)
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=NV_vec, size=NR_vec, rho=rhoj, prob=mu, log=T)))
  return(rhovec[ll==max(ll)][1])
}
beta.binom.filter = function(NR,NV,cutoff=0.1, binom.pval=F,pval.cutoff=0.05){
  # Function to apply beta-binomial filter for artefacts. Works best on sets of
  # clonal samples (ideally >10 or so). As before, takes NV and NR as input. 
  # Optionally calculates pvalue of likelihood beta-binomial with estimated rho
  # fits better than binomial. This was supposed to protect against low-depth variants,
  # but use with caution. Returns logical vector with good variants = TRUE
  
  rho_est = pval = rep(NA,nrow(NR))
  for (k in 1:nrow(NR)){
    rho_est[k]=estimateRho_gridml(NV_vec = as.numeric(NV[k,]),
                                  NR_vec = as.numeric(NR[k,]))
    if (k%%1000==0){
      print(k)
    }
  }
  flt_rho=rho_est>=cutoff
  return(rho_est)
}

NR30 <- aa30[,deps]
NV30 <- aa30[,mtrs]
NR10 <- aa10[,deps]
NV10 <- aa10[,mtrs]

betabinomout30 <- beta.binom.filter(NR=NR30, NV=NV30, cutoff=0.1, binom.pval = F, pval.cutoff = 0.05)
betabinomout10 <- beta.binom.filter(NR=NR10, NV=NV10, cutoff=0.1, binom.pval = F, pval.cutoff = 0.05)
filtmat$betabinom30 <- betabinomout30
filtmat$betabinom10 <- betabinomout10

filtmat$pass_betabinom <- filtmat$betabinom10>=0.05 & filtmat$betabinom30 >=0.05
filtmat$pass_betabinom[is.na(filtmat$pass_betabinom)] <- FALSE
filtmat$assign_to_tumour_tree <- FALSE
filtmat$assign_to_tumour_tree[filtmat$clonal_in_normals=='no' & filtmat$pass_betabinom==TRUE] <- TRUE
write.table(filtmat, 'mutations_with_filtering_status_after_betabinom.txt', sep='\t', col.names = T, quote=F)

# call mutations for the purpose of building the tree.
gdsom <- aa30[aa30$mutkey %in% filtmat$mutkey[filtmat$assign_to_tumour_tree==TRUE],]
write.table(gdsom, 'PD37590_mq30bq30_muts_pass_germline_and_betabinom_filters.txt',  ######################### change line with different patients
            sep='\t', col.names = T, row.names = F, quote=F)

# build the tree outside R. Phylip seqboot, then mix (Camin-Sokal), the consense.
# get the file into the correct format for phylip.
ND <- gdsom[,deps]
NV <- gdsom[,mtrs]

treemat <- NV
treemat[NV==0] <- 0
treemat[NV>3] <- 1
treemat[NV > 0 & NV <=3] <- NA
treemat[ND <10 & NV==0] <- NA
treemat$present <- rowSums(treemat, na.rm=T)
treemat$absent <- apply(treemat[,mtrs], 1, function(row) length(which(row==0 & !is.na(row))))
gdtreemat <- treemat[treemat$present>1 & treemat$absent>0,mtrs]
colnames(gdtreemat) <- gsub('_MTR', "", colnames(gdtreemat))
write.table(gdtreemat, 'informative_calls_matrix_for_treebuilding.txt', sep='\t', col.names=T, row.names=F, quote=F)

infile <- gdtreemat
infile[is.na(infile)] <- '?'
# add in a root of several mutations so that with bootstrapping it is still likely to be present
toadd <- matrix(nrow=10, ncol=ncol(infile), rep(1,ncol(infile)*10))
colnames(toadd) <- colnames(infile)
infile <- as.data.frame(rbind(toadd, infile))
infile$root000000 <- rep(0, nrow(infile))
infile <- t(infile)
chars_to_add <- unique(10 - nchar(rownames(infile))[1:length(samps)])
if (chars_to_add>0) {
  rownames(infile)[1:length(samps)] <- paste0(rownames(infile)[1:length(samps)], paste0(rep('0', chars_to_add), collapse=''))
}
if (chars_to_add<0) {
  rownames(infile)[1:length(samps)] <- substr(rownames(infile)[1:length(samps)], (unique(nchar(rownames(infile)[1:length(samps)]))-9), unique(nchar(rownames(infile)[1:length(samps)])))
}
dir.create('phylip_mtr4_dep10')
write.table(infile, 'phylip_mtr4_dep10/infile', sep='', col.names = F, row.names=T, quote=F)

# manually add the dimensions of the file to the sample. 
dim(infile)
# I have installed phylip here:  ~/Applications_installed_by_Henry/phylip-3.695. run seqboot, then mix, then consense

# I have built the tree in the terminal.
# now read it back in.
tree <- read.tree(paste0(getwd(), '/phylip_mtr4_dep10/consense/outtree'))
plot(tree)

# assign mutations to branches.
# use Nick Williams' script to assign mutations to tree. https://github.com/NickWilliamsSanger/treemut
rtree <- tree
if (chars_to_add > 0) {
  rtree$tip.label[rtree$tip.label!='root000000'] <- gsub('00','', rtree$tip.label[rtree$tip.label!='root000000'])
}
if (chars_to_add < 0) {
  rtree$tip.label[rtree$tip.label!='root000000'] <- paste0('PD375', rtree$tip.label[rtree$tip.label!='root000000'])
}
rtree$tip.label[rtree$tip.label=='root000000'] <- 'zeros'
rtree$edge.lengths <- rep(1, length(rtree$edge))
df <- reconstruct_genotype_summary(rtree)

# assign only high quality mutations, supported by at least 4 mtrs in at least one sample and vaf 0.1, plus passing the betabinomial and germline filters
toassign <- gdsom[,mtrs]
toassign[gdsom[,mtrs]<=3] <- 0
toassign[gdsom[,mtrs]>3] <- 1
toassign[gdsom[,vafs]<0.1] <- 0
colnames(toassign) <- gsub("_MTR", "", colnames(toassign))
sort(colSums(toassign)) 
head(toassign)
toassign$present <- rowSums(toassign)
table(toassign$present)
rownames(toassign) <- gdsom$mutkey
toassign$mutkey <- gdsom$mutkey
toassign <- toassign[toassign$present>0,]
write.table(toassign, 'mutations_to_assign_to_tree_4mtrs_vaf01_to_call.txt', sep='\t', col.names = T, row.names = F, quote=F)

# now need data frames of MTRs, DEPs, and a vector of error rates per sample
keepers <- which(gdsom$mutkey %in% toassign$mutkey)
nv <- gdsom[keepers,mtrs]
nr <- gdsom[keepers,deps]
colnames(nv) <- gsub('_MTR', '', colnames(nv))
colnames(nr) <- gsub('_DEP', '', colnames(nr))
errors <- c(rep(1e-02, length(mtrs)), 1e-06)
assigndat <- list(nv, nr, errors)
names(assigndat) <- c('mtrs', 'deps', 'errors')
res=assign_to_tree(tree=rtree, mtr=as.matrix(assigndat$mtrs), depth=as.matrix(assigndat$deps), error_rate = assigndat$errors, maxits=10)
tree_estimated <- res$tree

pdf("PD37590_tree_built_with_muts_not_clonal_in_normal_in_vitro_not_filtered.pdf")
plot(tree_estimated)
axisPhylo(side=1, backward=F)
dev.off()
write.tree(tree_estimated, 'PD37590_tree_built_with_muts_not_clonal_in_normal_in_vitro_not_filtered.tree')

# write out which mutations are assigned to which branch
assignments <- res$summary
# the label for the branch and the label for the node follow different systems.
# every the row number of the tree$edge matrix is the edge
# the numbers in this row are the nodes upstream and downstream of the edge.
# so for every edgelabel, I need to find the downstream node, and find its descendants.
edgelabs <- unique(assignments$edge_ml)
nodelabs <- sapply(edgelabs, function(cell) tree_estimated$edge[cell,2])
# now get the descendants of every node.
getDescendants.fn <- function(tree, node) {
  if (node <= length(tree$tip.label)) {
    return(tree$tip.label[node])
  }
  allnodes <- sort(unique(as.vector(tree$edge))[unique(as.vector(tree$edge)) > length(tree$tip.label)])
  tdescs <- tree$edge[tree$edge[,1]==node,2]
  tnodes <- tdescs[tdescs %in% allnodes]
  ttips <- tdescs[!tdescs %in% allnodes]
  # the number of nodes could be 0, 1, or 2.
  # I want to make a list of all the nodes to explore
  # if the descendant is a node, keep on going. If it's a tip, can stop.
  while(length(tnodes)>0) {
    thisnode <- tnodes[1]
    newdescs <- tree$edge[tree$edge[,1]==thisnode,2]
    nodestoadd <- newdescs[newdescs %in% allnodes]
    newtips <- newdescs[!newdescs %in% allnodes]
    # add the new descendant nodes to the list. 
    # now remove the node that I have just explored from the list. 
    tnodes <- tnodes[tnodes != thisnode]
    tnodes <- c(tnodes, nodestoadd)
    ttips <- c(ttips, newtips)
  }
  return(sort(tree$tip.label[ttips]))
}
tiplabs <- sapply(nodelabs, function(mynode) paste(getDescendants.fn(tree=tree_estimated, node=mynode), collapse=";"))
tipkey <- as.data.frame(cbind(edgelabs, nodelabs, tiplabs))
colnames(tipkey) <- c('edge', 'downstream_node', 'descendant_tips')
myres <- (res$summary)
myres$edge <- myres$edge_ml

# add the mutations.
ad <- as.data.frame(cbind(toassign, myres))

# now link up with the key
mmres <- merge(ad, tipkey)
mmres$edge <- NULL
write.table(mmres, 'assigned_mutations.txt', sep='\t', col.names=T, row.names = F, quote=F)
sort(table(mmres$descendant_tips))

# as these are not organoids, there are no in vitro mutations to filter out. But adding a few lines to harmonise with cases for which I have had to remove in vitro mutations.
mmres$invitro <- FALSE

# replot tree with in vitro mutations removed.
atree <- tree_estimated
atree$tip.label[atree$tip.label=='zeros'] <- 'fertilised_egg'
atree$edge.length

tipcols <- rep('grey', length(atree$tip.label))
tipcols[atree$tip.label %in% c('fertilised_egg')] <- 'black'

pdf('PD37590_phylogeny_in_vitro_mutations_removed_no_trunk.pdf')
plot(atree, font=2, direction='downwards', tip.color=tipcols, label.offset=10, cex=1)
axisPhylo(side=2, backward = F)
dev.off()



#####################################################################
# now that the phylogeny is built from the non-clonal-in-bulk-normal mutations, I need to go back to the clonal-in-bulk-normal-mutations
# there may be some early embryonic mutations that appeared to be clonal in the bulk normal. I need to fish these back out and inspect them.

# look at the vafs in the bulk normal sample.
# Plot vaf in blood against VAF in tumour.
# work out the vafs, summing across the tumour organoids, of all of the mutations.
tum_orgs <- samps
filtmat$tumour_mtr <- rowSums(aa30[,paste0(tum_orgs, '_MTR')])
filtmat$tumour_dep <- rowSums(aa30[,paste0(tum_orgs, '_DEP')])
filtmat$tumour_vaf <- filtmat$tumour_mtr/filtmat$tumour_dep

#	Take as embryonic those that are clearly off the curve. 
par(mfrow=c(2,1))
hist(filtmat$tumour_dep[grep('X|Y',filtmat$mutkey, invert = T)], 100, xlim=c(0,1000), main='autosomal mutations', xlab='depth')
hist(filtmat$tumour_dep[grep('X|Y',filtmat$mutkey)], 100, xlim=c(0,1000), main='sex chr mutations', xlab='depth')

filtmat$cov_in_tumour <- 'bad'
filtmat$cov_in_tumour[filtmat$tumour_dep>300 & filtmat$tumour_dep<600] <- 'good'

write.table(filtmat, 'filters_matrix.txt', sep='\t', col.names = T, row.names = F, quote=F)

# plot the vafs
# highlight the mutations that are in the trunk of the tumour
tum_trunk_muts <- mmres[mmres$descendant_tips==paste(sort(samps), collapse=';'),'mutkey']
tum_subclonal_muts <- mmres[!mmres$descendant_tips %in% paste(sort(samps), collapse=';'),'mutkey']

pdf('tumour_vaf_vs_normals_vaf_gd_coverage_regions.pdf')
plot(tumour_vaf~norm_vaf, data=filtmat[filtmat$clonal_in_normals!='poor_coverage_in_normal',], 
     pch=16, cex=0.7, xlim=c(0,1), ylim=c(0,1), main='Tumour vs normal vafs \n Good coverage areas')
points(tumour_vaf~norm_vaf, data=filtmat[filtmat$mutkey %in% tum_trunk_muts,], 
       pch=16, cex=1.2, col='red', main='all muts good coverage')
points(tumour_vaf~norm_vaf, data=filtmat[filtmat$mutkey %in% tum_subclonal_muts,], 
       pch=16, cex=1.2, col='brown', main='all muts good coverage')
points(tumour_vaf~norm_vaf, data=filtmat[filtmat$pass_betabinom=='FALSE' & filtmat$clonal_in_normals=='no',], 
       pch=16, cex=1.2, col='blue', main='all muts good coverage')
abline(v=vafcutoff, col='red', lwd=2)
legend('topleft', pch=1, cex=1, col=c('red', 'brown', 'blue'),
       legend=c('Tumour trunk', 'Tumour branches', 'Failed betabinom'), box.lty=0)
dev.off()

# this looks like all the useful mutations have been captured.

# inspect the mutations that appear clonal in the tumour and subclonal in the blood that have been filtered out by the beta binomial.
in_tum_failed_betabinom <- filtmat[filtmat$pass_betabinom==FALSE & filtmat$tumour_vaf>0.4 & filtmat$clonal_in_normals=='no',]
dim(in_tum_failed_betabinom)

write.csv(in_tum_failed_betabinom, 'subclonal_in_norm_clonal_in_tumour_failed_betabinom_to_check.csv', quote=F)

# of these 6, after manual review on jbrowse, none are good truncal mutations. Don't need to add anything back.


##############
# finalise set of mutations
btree <- atree
plot(btree)

write.tree(btree, 'PD37590_tree_failed_betabinom_salvaged.tree')

pdf('PD37590_phylogeny_failed_betabinom_salvaged.pdf')
plot(btree, font=2, direction='downwards', tip.color=tipcols, label.offset=10, cex=1)
axisPhylo(side=2, backward = F)
dev.off()

# make a file of all the mutations in the original pileup, which filters they passed or failed, and if they passed where they're assigned 
# add trinucleotides to this file
# write out all mutations
aa30$pass_germline_filter <- 'fail'
aa30$pass_betabinomial_filter <- 'fail'

aa30$pass_germline_filter[aa30$mutkey %in% filtmat$mutkey[filtmat$clonal_in_normals=='no']] <- 'pass'
aa30$pass_betabinomial_filter[aa30$mutkey %in% filtmat$mutkey[filtmat$pass_betabinom==TRUE]] <- 'pass'

mmres_to_merge <- mmres[,c('mutkey', 'edge_ml', 'descendant_tips', 'invitro')]
muts_with_status <- merge(aa30, mmres_to_merge, all.x=T, all.y=T)

# sort out the ones I salvaged.
table(mmres$edge_ml)

# now calculate trinucleotide context for every mutation. 
gr <- GRanges(seqnames = paste0("chr", muts_with_status$Chrom), IRanges(start = muts_with_status$Pos, width=1), ref=muts_with_status$Ref, alt=muts_with_status$Alt)
if (all(substr(seqlevels(gr), 1, 3) != "chr")) {
  gr <- renameSeqlevels(gr, paste0("chr", seqlevels(gr)))
}
genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
seqinfo(gr) <- seqinfo(genome)[seqlevels(gr)]
gr <- sort(gr)
bases <- c("A", "C", "G", "T")
trinuc_levels <- paste0(rep(bases, each = 16), rep(rep(bases, each = 4), 4), rep(bases, 16))
get_trinuc <- function(seqname) {
  pos <- start(gr[seqnames(gr) == seqname])
  view <- Views(genome[[seqname]], start = pos - 1, end = pos + 1)
  ans <- factor(as.character(view), levels = trinuc_levels, labels = 1:64)
  return(as.numeric(ans))
}
trinuc <- sapply(seqlevels(gr), get_trinuc)
gr$trinuc <- factor(unlist(trinuc, use.names = FALSE), levels = 1:64, labels = trinuc_levels)
remove(trinuc)
gr$REF <- gr$ref
gr$ALT <- gr$alt
gr$context <- gr$trinuc
torc <- which(gr$ref %in% c("A", "G"))
gr$REF[torc] <- as.character(reverseComplement(DNAStringSet(gr$REF[torc])))
gr$ALT[torc] <- as.character(reverseComplement(DNAStringSet(gr$ALT[torc])))
gr$context[torc] <- as.character(reverseComplement(DNAStringSet(gr$context[torc])))
gr$class <- paste(gr$REF, gr$ALT, "in", gr$context, sep = ".")
class_levels <- paste(rep(c("C", "T"), each = 48), rep(c("A", "G", "T", "A", "C", "G"), each = 16), "in", paste0(rep(rep(bases, each = 4), 6), rep(c("C", "T"), each = 48), rep(bases, 24)), sep = ".")
gr$class <- factor(gr$class, levels = class_levels)

# match them up with one another using coordinates - I think that they may have been reshuffled.
grdf <- as.data.frame(gr)
grdf$seqnames <- gsub('chr', '', grdf$seqnames)
grdf$mutkey <- paste0(grdf$seqnames, '_', grdf$start, '_', grdf$ref, '_', grdf$alt)
muts <- merge(muts_with_status, grdf[,c('mutkey', 'class')])

write.table(muts, 'PD37590_pileup_with_filters_assignments_trinucs.txt', sep='\t', col.names=T, row.names=F, quote=F)

# write out mutations to check
edgestocheck <- names(sort(table(muts$edge_ml[muts$invitro==FALSE & !is.na(muts$invitro)])))
edge <- edgestocheck[1]
num_to_check <- 2
tocheck <- data.frame()
for (edge in edgestocheck) {
  tedge <- muts[muts$edge_ml==edge & !is.na(muts$edge_ml),c('mutkey', 'edge_ml', 'descendant_tips')]
  if (nrow(tedge)>num_to_check) {
    tedge <- tedge[sample(1:nrow(tedge), num_to_check, replace = F),]
  }
  tocheck <- as.data.frame(rbind(tocheck, tedge))
}
write.csv(tocheck, 'mutations_to_check_by_branch.csv', quote=F)

# plot trinucleotides per filter and per branch
trincols <- rep(c('blue', 'black', 'red', 'grey', 'darkseagreen', 'pink'), each=16)
muts$invitro <- FALSE

pdf('PD37590 trinucleotide plots by filter and by branch.pdf')
par(mfrow=c(2,1))
barplot(table(muts$class[muts$pass_germline_filter=='pass']), col=trincols, border = NA, las=2, cex.names = 0.2, main='pass germline filter')
barplot(table(muts$class[muts$pass_germline_filter=='fail']), col=trincols, border = NA, las=2, cex.names = 0.2, main='fail germline filter')
barplot(table(muts$class[muts$pass_germline_filter=='pass' & muts$pass_betabinomial_filter=='pass']), col=trincols, border = NA, las=2, cex.names = 0.2, main='pass germline and betabinom filters')
barplot(table(muts$class[muts$pass_germline_filter=='pass' & muts$pass_betabinomial_filter=='fail']), col=trincols, border = NA, las=2, cex.names = 0.2, main='pass germline but fail betabinom filters')
par(mfrow=c(1,1))
barplot(table(muts$class[muts$pass_germline_filter=='pass' 
                         & (muts$pass_betabinomial_filter %in% c('pass',"salvaged_by_jbrowse"))
                         & muts$invitro==FALSE 
                         & !is.na(muts$invitro) 
                         & !is.na(muts$edge_ml)]), col=trincols, border = NA, 
        las=2, cex.names = 0.2, main='All mutations assigned to the tree')
# now plot per branch.
par(mfrow=(c(3,1)))
for (edge in edgestocheck) {
  barplot(table(muts$class[muts$pass_germline_filter=='pass' 
                           & (muts$pass_betabinomial_filter %in% c('pass',"salvaged_by_jbrowse"))                           
                           & muts$invitro==FALSE 
                           & !is.na(muts$invitro) 
                           & !is.na(muts$edge_ml) &
                             muts$edge_ml==edge]), col=trincols, border = NA, 
          las=2, cex.names = 0.2, cex.main=0.7,
          main=paste0('muts in ', muts$descendant_tips[muts$edge_ml==edge & !is.na(muts$edge_ml)][1]))
}
dev.off()



# have a look at mutations by vaf
treemuts <- muts[muts$pass_germline_filter=='pass' & muts$pass_betabinomial_filter=='pass' & muts$invitro==FALSE & !is.na(muts$invitro),]
dim(treemuts)
# calculate the max vaf and max mtrs of every mutation.
treemuts$maxmtr <- apply(treemuts[,paste0(samps, '_MTR')], 1, function(row) max(row))
treemuts$maxvaf <- apply(treemuts[,paste0(samps, '_VAF')], 1, function(row) max(row))
treemuts$mtrbins <- cut(treemuts$maxmtr, 5)
treemuts$vafbins <- cut(treemuts$maxvaf, 5)

pdf('mutations assigned to tree by mtr and vaf.pdf')
layout(matrix(1:10, nrow=5, ncol=2, byrow=F))
par(mar=c(1,4,4,1))
for (mtrbin in sort(unique(treemuts$mtrbins))) {
  tbinclass <- treemuts$class[treemuts$mtrbins == mtrbin]
  barplot(table(tbinclass), col=trincols, border = NA, las=2, cex.names = 0.2, main=paste0('mtr bin ', mtrbin), names.arg='')
}
for (vafbin in sort(unique(treemuts$vafbins))) {
  tbinclass <- treemuts$class[treemuts$vafbins == vafbin]
  barplot(table(tbinclass), col=trincols, border = NA, las=2, cex.names = 0.2, main=paste0('vaf bin ', vafbin), names.arg='')
}
dev.off()


# look at the clonality fo samples
# muts <- read.csv('PD37590_pileup_with_filters_assignments_trinucs.txt', sep='\t', header=T, stringsAsFactors = F)
head(muts)
dim(muts)

# for every branch plot clonal and subclonal mutations
samp <- samps[1]

pdf('PD37590_vafs_shared_vs_private.pdf')
par(mfrow=c(2,1))
for (samp in samps) {
  tmuts <- muts[grep(samp, muts$descendant_tips),]
  private <- tmuts[tmuts$descendant_tips==samp,]
  shared <- tmuts[tmuts$descendant_tips!=samp,]
  
  hist(shared[,paste0(samp, "_VAF")], 50, xlim=c(0,1), main=paste0(samp, ' shared muts autosomes'))
  hist(private[,paste0(samp, "_VAF")], 50, xlim=c(0,1), main=paste0(samp, ' private muts autosomes'))
}
dev.off()

