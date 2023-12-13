# 2023.12.04
# wilms organoids 
# phylogeny construction

library('GenomicRanges')
library('BSgenome.Hsapiens.UCSC.hg19')
library("dplyr")
library("VGAM")
library('ape')
library('treemut')


setwd('/Users/hl11/Documents/Research/Projects/wilms_organoids/phylogenies_same_method_dec23/')

patient <- 'PD49348' # can code this so can come in from the command line
setwd(paste0(getwd(), '/', patient))

# read in the pileup
aa30 <- read.csv(paste0(patient,'_pileup_mq30bq30.tsv'), sep='\t', header=T, stringsAsFactors = F)
aa30$mutkey <- paste0(aa30$Chrom, '_', aa30$Pos, '_', aa30$Ref, '_', aa30$Alt)

# keep in the single cell derived organoids ('orgs). b and k are normal organoids. However, k is polyclonal, so I am removing it.
tum_orgs <- c('PD49348a', 'PD49348e', 'PD49348f', 'PD49348g', 'PD49348h', 'PD49348i', 'PD49348j')
norm_orgs <- c('PD49348b', 'PD49348d')
samps <- c(tum_orgs, norm_orgs)
normals <- c("PD49348v") # we do have bulk normal kidney, but I think leave it out as it could lead to me losing embryonic kidney-specific mtuations

# get useful vectors of column names to work on the pileup.
mtrs <- paste0(samps, '_MTR')
deps <- paste0(samps, '_DEP')
vafs <- paste0(samps, '_VAF')
ofs <- paste0(samps, '_OFS')

# calculate the vaf of every mutation across all normal samples.
filtmat <- aa30[,c('mutkey', paste0(rep(normals, each=2), c('_MTR', '_DEP')))]
filtmat$norm_mtr <- (filtmat[,grep('MTR', colnames(filtmat))])
filtmat$norm_dep <- (filtmat[,grep('DEP', colnames(filtmat))])
filtmat$norm_vaf <- filtmat$norm_mtr/filtmat$norm_dep

# only keep mutations that are well covered in the bulk normals. If they are poorly covered, they could be dodgy sites, or we can't confidently exclude them from being germline.
# set a cutoff for good depth for the normal samples. This is done manually based on the histogram for this sample.
hist(filtmat$norm_dep, 100, xlim=c(0,500))
abline(v=100, col='red')
abline(v=350, col='red')
norm_cutoff_min <- 100
norm_cutoff_max <- 350

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
print(vafcutoff)

# now separate our mutations into those that are clonal in the bulk normal or not
filtmat$clonal_in_normals <- 'yes'
filtmat$clonal_in_normals[filtmat$norm_vaf<vafcutoff & filtmat$norm_dep>norm_cutoff_min & filtmat$norm_dep<norm_cutoff_max] <- 'no'
# any that were poorly covered in the normals it is difficult to comment on, and I don't want to use them for tree building.
filtmat$clonal_in_normals[filtmat$norm_dep<=norm_cutoff_min | filtmat$norm_dep >= norm_cutoff_max] <- 'poor_coverage_in_normal'
table(filtmat$clonal_in_normals) 


# now decide if enriched in tumour compared to normal
filtmat$tumour_mtr <- rowSums(aa30[,paste0(tum_orgs, '_MTR')])
filtmat$tumour_dep <- rowSums(aa30[,paste0(tum_orgs, '_DEP')])
filtmat$tumour_vaf <- filtmat$tumour_mtr/filtmat$tumour_dep

#	Take as embryonic those that are clearly off the curve. 
hist(filtmat$tumour_dep, 100, xlim=c(0,300))

par(mfrow=c(2,1))
hist(filtmat$tumour_dep[grep('X|Y',filtmat$mutkey, invert = T)], 100, xlim=c(0,300), main='autosomal mutations', xlab='depth')
hist(filtmat$tumour_dep[grep('X|Y',filtmat$mutkey)], 100, xlim=c(0,300), main='sex chr mutations', xlab='depth')
# this is a female, so there is no difference between sex chroms and autosomes. 

# don't need to worry about cn states as clonally diploid
filtmat$cn_state <- '2_1'

filtmat$cov_in_tumour <- 'bad'
filtmat$cov_in_tumour[(!grepl('X|Y',filtmat$mutkey)) & filtmat$tumour_dep>50 & filtmat$tumour_dep<300] <- 'good'
filtmat$cov_in_tumour[(grepl('X|Y',filtmat$mutkey)) & filtmat$tumour_dep>50 & filtmat$tumour_dep<300] <- 'good'


#### I will accept mutations that:
# feature at below the vaf cutoff in the normal - this removes germline
# are enriched in the tumour compared to the normal, if they have high vaf  - this removes germline & recurrent artefacts but salvages embryonic mutations that are present at high-ish vaf in the normal but are enriched in the tumour.
#### BUT: I don't want to apply this to mutations that are at low vaf in the tumour, e.g. just called in a single tumour sample.
#### AND: it must not apply, of course, to mutations that I am trying to call in the normal samples. 
########## For this latter condition I am seeing if applying a two sided test, along with the low vaf derogation, fixes the problem.
# are not in a region of CN change - this removes germline muts that are higher vaf in the tumour because of a cn change.
filtmat$norm_vaf[is.na(filtmat$norm_vaf)] <- 1 # set NAs as 1 as these are regions that I don't want to use
filtmat$norm_vaf[filtmat$norm_vaf==0] <- 1e-06
filtmat[filtmat$tumour_dep<filtmat$tumour_mtr,]
filtmat$tumour_dep[filtmat$tumour_dep==0] <- 1
filtmat$tumour_dep <- as.integer(filtmat$tumour_dep)
filtmat$tumour_mtr <- as.integer(filtmat$tumour_mtr)
filtmat$tum_vs_normal_binom_pval <- apply(filtmat, 1, function(row) {
  binom.test(x=as.numeric(row['tumour_mtr']), n=as.numeric(row['tumour_dep']), p=as.numeric(row['norm_vaf']), alternative = 'two.sided')$p.val
})

# now apply a multiple testing correction
filtmat$tum_vs_normal_binom_qval <- p.adjust(p=filtmat$tum_vs_normal_binom_pval, method='hochberg')

# now apply a derogation for the examples listed above (low vaf)
filtmat$tum_vs_normal_binom_qval[filtmat$norm_vaf <0.1 & filtmat$tumour_vaf < 0.1] <- NA
length(which(is.na(filtmat$tum_vs_normal_binom_qval)))

write.table(filtmat, paste0(patient, '_filters_outcome_with_cn_state_and_binom_test.txt'), sep='\t', col.names = T, row.names=F, quote=F)


pdf(paste0(patient, '_tumour_vaf_vs_normals_vaf_gd_coverage_regions.pdf'))
plot(tumour_vaf~norm_vaf, data=filtmat[filtmat$clonal_in_normals!='poor_coverage_in_normal',], 
     pch=16, cex=0.7, xlim=c(0,1), ylim=c(0,1), main='Tumour vs normal vafs \n Good coverage areas')
points(tumour_vaf~norm_vaf, data=filtmat[filtmat$pass_betabinom=='FALSE' & filtmat$clonal_in_normals=='no',], 
       pch=1, cex=1.2, col='blue')
abline(v=vafcutoff, col='red', lwd=2)
abline(h=vafcutoff*2, col='red', lwd=2, lty=3)
legend('topleft', pch=1, cex=1, col=c('blue'),
       legend=c('Failed betabinom'), box.lty=0)

plot(tumour_vaf~norm_vaf, data=filtmat[filtmat$clonal_in_normals!='poor_coverage_in_normal',], 
     pch=16, cex=0.7, xlim=c(0,1), ylim=c(0,1), main='Tumour vs normal vafs \n Good coverage areas')
points(tumour_vaf~norm_vaf, data=filtmat[filtmat$tum_vs_normal_binom_qval<0.01 & filtmat$clonal_in_normals=='no' & filtmat$cn_state=='2_1',], 
       pch=1, cex=1.2, col='orange')
abline(v=vafcutoff, col='red', lwd=2)
abline(h=vafcutoff*2, col='red', lwd=2, lty=3)
legend('topleft', pch=1, cex=1, col=c('orange'),
       legend=c('Enriched tum vs norm'), box.lty=0)

plot(tumour_vaf~norm_vaf, data=filtmat[filtmat$clonal_in_normals!='poor_coverage_in_normal',], 
     pch=16, cex=0.7, xlim=c(0,1), ylim=c(0,1), main='Tumour vs normal vafs \n Good coverage areas')
points(tumour_vaf~norm_vaf, data=filtmat[filtmat$cn_state != '2_1' & filtmat$clonal_in_normals=='no',], 
       pch=1, cex=1.2, col='red')
abline(v=vafcutoff, col='red', lwd=2)
abline(h=vafcutoff*2, col='red', lwd=2, lty=3)
legend('topleft', pch=1, cex=1, col=c('red'),
       legend=c('In non-diploid regions'), box.lty=0)
dev.off()


# add a few more filters
## the mutation should be called in at least one sample at a decent vaf and mtrs
# as part of this, apply Tolly's directionality filter. I can use the MDR column. MDR=3 means that it was present on both forward and reverse. 
# if this doesn't remove enough, I could demand a certain skew ratio - could do this based on all the mtrs

num_called <- apply(aa30[,paste0(rep(samps, each=3), c('_MTR', '_VAF', '_MDR'))], 1 , function(row) {
  as.numeric(row[(1:length(samps)*3 - 2)] >= 4) + as.numeric(row[(1:length(samps)*3 -1)]>= 0.3) + as.numeric(row[(1:length(samps)*3)]== 3)
})
num_called[1:5,1:10]
nc <- apply(num_called, 2, function(row) length(which(row==3 & !is.na(row))))
hist(nc)
filtmat$samples_called <- nc
head(filtmat)

#### check how many mutations fail the skew filter
#### I need to do this just for the mutant base each time
forward_mtrs <- apply(aa30,1, function(row) {
  sum(as.numeric(row[paste0(samps, "_F", row['Alt'], 'Z')]))
})
reverse_mtrs <- apply(aa30,1, function(row) {
  sum(as.numeric(row[paste0(samps, "_R", row['Alt'], 'Z')]))
})

for_rev <- as.data.frame(cbind(forward_mtrs, reverse_mtrs))
filtmat$for_rev_imbal_pval <- apply(for_rev, 1, function(row) binom.test(x=c(row[1], row[2]), p=0.5, alternative = 'two.sided')$p.val)
filtmat$for_rev_imbal_qval <- p.adjust(filtmat$for_rev_imbal_pval, method='hochberg')
hist(filtmat$for_rev_imbal_qval)
dim(filtmat[filtmat$for_rev_imbal_qval<0.05,])

### Now apply the filters to decide which mutations to use.
filtmat$assign_to_tumour_tree <- FALSE
filtmat$assign_to_tumour_tree[filtmat$clonal_in_normals=='no' & (filtmat$tum_vs_normal_binom_qval<0.01 | is.na(filtmat$tum_vs_normal_binom_qval)) & filtmat$samples_called>=1 & !grepl(';', filtmat$cn_state) & filtmat$for_rev_imbal_qval>0.05] <- TRUE
table(filtmat$assign_to_tumour_tree)

# call mutations for the purpose of building the tree.
gdsom <- aa30[aa30$mutkey %in% filtmat$mutkey[filtmat$assign_to_tumour_tree==TRUE],]
dim(gdsom)
write.table(gdsom, paste0(patient, '_muts_pass_filters_to_assign_to_tree.txt'), sep='\t', col.names = T, row.names = F, quote=F)

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
dim(infile) #  

# I have installed phylip here:  ~/Applications_installed_by_Henry/phylip-3.695/exe. run seqboot, then mix, then consense
# I have built the tree in the terminal.
# now read it back in.
tree <- read.tree(paste0(getwd(), '/phylip_mtr4_dep10/consense/outtree'))
dev.off()
plot(tree)

# assign mutations to branches.
# use Nick Williams' script to assign mutations to tree. https://github.com/NickWilliamsSanger/treemut
rtree <- tree
if (chars_to_add > 0) {
  rtree$tip.label[rtree$tip.label!='root000000'] <- gsub('00','', rtree$tip.label[rtree$tip.label!='root000000'])
}
if (chars_to_add < 0) {
  rtree$tip.label[rtree$tip.label!='root000000'] <- paste0('PD548', rtree$tip.label[rtree$tip.label!='root000000'])
}
rtree$tip.label[rtree$tip.label=='root000000'] <- 'zeros'
rtree$edge.lengths <- rep(1, length(rtree$edge))
df <- reconstruct_genotype_summary(rtree)

# assign only high quality mutations. Genotype on the basis of 4 mtrs and vaf 0.1.
toassign <- gdsom[,mtrs]
toassign[gdsom[,mtrs]<=3] <- 0
toassign[gdsom[,mtrs]>3] <- 1
toassign[gdsom[,vafs]<0.1] <- 0
colnames(toassign) <- gsub("_MTR", "", colnames(toassign))
sort(colSums(toassign)) 
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

pdf(paste0(patient, "_tree_built_with_muts_not_clonal_in_normal_in_vitro_not_filtered.pdf"))
plot(tree_estimated)
axisPhylo(side=1, backward=F)
dev.off()
write.tree(tree_estimated, paste0(patient, '_tree_built_with_muts_not_clonal_in_normal_in_vitro_not_filtered.tree'))


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


# filter out in vitro mutations.
# for every mutation, I can calculate the probability of seeing it at the observed vaf if it is actually clonal.
# I should cast out mutations there there is a low probability of seeing it at this vaf. 
# NB this is the other way round to a pvalue, as I am setting a genuine clonal mutation to being the null hypothesis. And the higher the probability, the more likely the mutation is to be real.
mmres$invitro <- FALSE
invitro_pval <- 0.1

pdf('Filtering_out_in_vitro_mutations.pdf', height=10, width=10)
layout(matrix(c(1:6),nrow=3, ncol=2, byrow=F))
for (samp in samps) {
  print(samp)
  tmres <- mmres[grep(samp, mmres$descendant_tips),]
  tprivate <- tmres[tmres$descendant_tips==samp,]
  tshared <- tmres[tmres$descendant_tips!=samp,]
  
  # all the organoids should be clonal
  tvaf <- 0.5
  
  # work out the probability of a private mutation being clonal.
  tmtrs <- aa30[aa30$mutkey %in% tprivate$mutkey, paste0(samp, '_MTR')]
  tdeps <- aa30[aa30$mutkey %in% tprivate$mutkey, paste0(samp, '_DEP')]
  tmuts <- as.data.frame(cbind(tmtrs, tdeps))
  tmuts$tvafs <- tmtrs/tdeps
  tmuts$mutkey <- aa30[aa30$mutkey %in% tprivate$mutkey, 'mutkey']
  
  # for every mutation, work out the probability that it could come from the mean vaf.
  tmuts$prob_vaf_this_low <- apply(tmuts[,c(1,2)], 1, function(row) {
    trands <- rbinom(size=row[2], prob=tvaf, n=10000)
    return(length(which(trands <= row[1]))/length(trands))
  })
  
  if (nrow(tshared) >0) {
    hist(aa30[aa30$mutkey %in% tshared$mutkey, paste0(samp, '_VAF')], xlim=c(0,1), 50, 
         main=paste0(samp, ' shared mutations'), xlab='VAF')
    abline(v=tvaf, col='red', lwd=2)
  } else {
    plot(1~1, type='n', axes=F, xlab='', ylab='')
    text(1,1, 'No shared mutations')
  }
  hist(aa30[aa30$mutkey %in% tprivate$mutkey, paste0(samp, '_VAF')], xlim=c(0,1), 50, 
       main=paste0(samp, ' private mutations \n ', length(tprivate$mutkey), ' subs unfiltered'), xlab='VAF')
  abline(v=tvaf, col='red', lwd=2)
  hist(aa30[aa30$mutkey %in% tmuts$mutkey[tmuts$prob_vaf_this_low>invitro_pval], paste0(samp, '_VAF')], xlim=c(0,1), 50, 
       main=paste0(samp, ' private mutations \n ', length(tmuts$mutkey[tmuts$prob_vaf_this_low>invitro_pval]), ' subs after filtering'), xlab='VAF')
  abline(v=tvaf, col='red', lwd=2)
  
  # now save which ones passed the filter
  mmres$invitro[(mmres$mutkey %in% tmuts$mutkey[tmuts$prob_vaf_this_low<invitro_pval])] <- TRUE
}
dev.off()

# replot tree with in vitro mutations removed.
atree <- tree_estimated
edgeml <- mmres$edge_ml[1]
for (edgeml in unique(mmres$edge_ml)) {
  atree$edge.length[as.numeric(edgeml)] <- nrow(mmres[mmres$edge_ml==edgeml & mmres$invitro==FALSE,])
}
atree$tip.label[atree$tip.label=='zeros'] <- 'fertilised_egg'
atree$edge.length

tipcols <- rep('red', length(atree$tip.label))
tipcols[atree$tip.label %in% norm_orgs] <- 'darkgrey'
tipcols[atree$tip.label %in% c('fertilised_egg')] <- 'black'

pdf(paste0(patient, '_phylogeny_in_vitro_mutations_removed_no_trunk.pdf'))
plot(atree, font=2, direction='downwards', tip.color=tipcols, label.offset=10, cex=1)
axisPhylo(side=2, backward = F)
dev.off()


# now calculate trinucleotide context for every mutation.
mmres_to_merge <- mmres[,c('mutkey', 'edge_ml', 'descendant_tips', 'invitro')]
muts_with_status <- merge(aa30, mmres_to_merge, all.x=T, all.y=T)

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
intersect(colnames(muts), colnames(filtmat))
muts <- merge(muts, filtmat[,!grepl(patient, colnames(filtmat))])
table(muts$descendant_tips)

write.table(muts, paste0(patient, '_pileup_with_filters_assignments_trinucs.txt'), sep='\t', col.names=T, row.names=F, quote=F)


# write out mutations to check
edgestocheck <- names(sort(table(muts$edge_ml[muts$invitro==FALSE & !is.na(muts$invitro)])))
edge <- edgestocheck[1]
num_to_check <- 5
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

pdf(paste0(patient, ' trinucleotide plots by filter and by branch.pdf'))
par(mfrow=c(2,1))
barplot(table(muts$class[muts$clonal_in_normals=='no']), col=trincols, border = NA, las=2, cex.names = 0.2, main='pass germline filter')
barplot(table(muts$class[muts$clonal_in_normals=='yes']), col=trincols, border = NA, las=2, cex.names = 0.2, main='fail germline filter')
barplot(table(muts$class[(muts$tum_vs_normal_binom_qval<0.01 | is.na(muts$tum_vs_normal_binom_qval)) & muts$clonal_in_normals=='no']), col=trincols, border = NA, las=2, cex.names = 0.2, main='pass tum vs norm filter')
barplot(table(muts$class[muts$tum_vs_normal_binom_qval>0.01 & !is.na(muts$tum_vs_normal_binom_qval & muts$clonal_in_normals=='no']), col=trincols, border = NA, las=2, cex.names = 0.2, main='fail tum vs norm filter')
        barplot(table(muts$class[muts$for_rev_imbal_qval<0.05]), col=trincols, border = NA, las=2, cex.names = 0.2, main='fail strand imbalance test')
        par(mfrow=c(1,1))
        barplot(table(muts$class[muts$assign_to_tumour_tree==TRUE]), col=trincols, border = NA, 
                las=2, cex.names = 0.2, main='All mutations assigned to the tree')
        # now plot per branch.
        par(mfrow=(c(3,1)))
        for (edge in edgestocheck) {
          barplot(table(muts$class[!is.na(muts$edge_ml) &
                                     muts$edge_ml==edge]), col=trincols, border = NA, 
                  las=2, cex.names = 0.2, cex.main=0.7,
                  main=paste0('muts in ', muts$descendant_tips[muts$edge_ml==edge & !is.na(muts$edge_ml)][1]))
        }
        # plot the trunk vs all the branches.
        
        
        # have a look at mutations by vaf
        treemuts <- muts[!is.na(muts$edge_ml),]
        # calculate the max vaf and max mtrs of every mutation.
        treemuts$maxmtr <- apply(treemuts[,paste0(samps, '_MTR')], 1, function(row) max(row))
        treemuts$maxvaf <- apply(treemuts[,paste0(samps, '_VAF')], 1, function(row) max(row))
        treemuts$mtrbins <- cut(treemuts$maxmtr, 5)
        treemuts$vafbins <- cut(treemuts$maxvaf, 5)
        
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
        
        par(mfrow=c(2,1))
        for (samp in samps) {
          tmuts <- muts[grep(samp, muts$descendant_tips),]
          private <- tmuts[tmuts$descendant_tips==samp,]
          shared <- tmuts[tmuts$descendant_tips!=samp,]
          
          hist(shared[!shared$Chrom %in% c("X", "Y"),paste0(samp, "_VAF")], 50, xlim=c(0,1), main=paste0(samp, ' shared muts autosomes'))
          if (nrow(private[!private$Chrom %in% c("X", "Y"),])==0) {
            plot(1~1, type='n', axes=F, xlab='', ylab='', main=paste0(samp, ': no private subs'))
            next
          }
          hist(private[!private$Chrom %in% c("X", "Y"),paste0(samp, "_VAF")], 50, xlim=c(0,1), main=paste0(samp, ' private muts autosomes'))
          
          if (nrow(shared[shared$Chrom %in% c("X", "Y"),])>0 & nrow(private[private$Chrom %in% c("X", "Y"),])>0) {
            hist(shared[shared$Chrom %in% c("X", "Y"),paste0(samp, "_VAF")], 50, xlim=c(0,1), main=paste0(samp, ' shared muts sex chroms'))
            hist(private[private$Chrom %in% c("X", "Y"),paste0(samp, "_VAF")], 50, xlim=c(0,1), main=paste0(samp, ' private muts sex chroms'))
          }
        }
dev.off()
        
