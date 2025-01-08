# Henry Lee-Six
# 2024.12.13
# wilms organoids projects
# simulate the true mutation burden of every cell in the tumour and the number of mutations captured by wgs and nanoseq
# do this with different parameters that read in as standard input

myline <- commandArgs(trailingOnly = TRUE)[1] # as an example, this might be '1', to use the parameters in line 1 of the parameters file.
print(myline)

# read in matrix of variables
vm = read.csv('variables_for_simulation.txt', header=T, stringsAsFactors = F, sep='\t')
# an example header and line from this file might be as follows:
# trunk_length	cell_number	bottleneck_prop	bottleneck_time	muts_per_doubling_div
# 10	65536	0.1	3	3

trow = vm[as.numeric(myline),]
print(trow)

sim_muts.fn = function(trunklen, cellnum, bottleneck_props, bottleneck_gens) {
  # make the subclonal mutations
  # matrix where every row is a cell and every column is a generation or mutation
  # iterate through generations.
  # if it's a doubling division, I need to increase the size of the grid.
  # if it's a bottleneck, I need to reduce it. 
  
  # make a starting grid.
  cellcounter = 1
  gencounter=1
  cellmat=(matrix(nrow=cellcounter, ncol=gencounter, c('0_1')))
  
  for (gen in 1:1000) {
    # print(paste0('generation: ', gen))
    # print(paste0('cells :', cellcounter))
    # stop when I have reached the number of cells I want
    if (cellcounter>cellnum) {
      break
    }
    
    if (gen %in% bottleneck_gens) {
      # print('bottleneck')
      keepers = sample(1:nrow(cellmat), size=ceiling(bottleneck_props*nrow(cellmat)), replace=F)
      newcellmat = matrix(cellmat[keepers,], nrow=length(keepers))
      cellcounter = nrow(cellmat)
      
      # need to add on a column for generation even though the bottleneck means that there hasn't been growth
      cellmat = cbind(cellmat, rep('dummy', nrow(cellmat)))
      
    } else {
      # make a new generation of cells to tack on. 
      cellcounter = cellcounter*2
      newgen = paste0(rep(gen, cellcounter), '_', 1:cellcounter)
      
      # double the size of the grid, and double every previous cell.
      newcellmat=matrix(nrow=cellcounter, ncol=gen+1)
      newcellmat[,gen+1] = newgen
      for (i in 1:gen) {
        newcellmat[,i] = rep(cellmat[,i], each=2)
      }
      cellmat=newcellmat
    }
  }
  
  # make the truncal mutations 
  # do this as a matrix. 
  trunk_muts = matrix(nrow=nrow(cellmat), ncol=trunklen, rep(paste0('t_', rep(1:trunklen)), each = nrow(cellmat)))
  allmuts = as.data.frame(cbind(trunk_muts, cellmat))
  
  return(allmuts)
}

# now, for different trunk lengths and models, calculate the wgs and nanoseq burdens. 
run_sim.fn = function(trunklen, cellnum, bottleneck_props, bottleneck_gens, muts_per_cell_division) {
  tsim = sim_muts.fn(trunklen=trunklen, cellnum=cellnum, bottleneck_props = bottleneck_props, bottleneck_gens = bottleneck_gens)
  
  # truth
  all_muts = unique(unlist(tsim))
  all_muts = all_muts[all_muts!= 'dummy']
  truncal_muts = length(which(substr(all_muts,1,1)=='t'))
  subclonal_muts = length(which(substr(all_muts,1,1)!='t')) * muts_per_cell_division
  
  total_muts = truncal_muts + subclonal_muts
  
  # wgs
  vafs = (table(as.vector(unlist(tsim)))/nrow(tsim))*0.5
  vafs=vafs[names(vafs)!='dummy']
  subclon_ids = substr(names(vafs),1,1)!= 't'
  trunk_ids = substr(names(vafs),1,1)== 't'
  vafs = as.vector(vafs)
  vafs = c(vafs[trunk_ids], rep(vafs[subclon_ids], each=muts_per_cell_division))
  deps=rpois(n=length(vafs), lambda=30)
  wgsmat=as.data.frame(cbind(vafs, deps))
  wgsmat$mtr = apply(wgsmat, 1, function(row) {
    rbinom(n=1, size=row['deps'], prob=row['vafs'])
  })
  wgsmat$obsvaf = wgsmat$mtr/wgsmat$deps
  wgscalled=nrow(wgsmat[wgsmat$mtr>=4 & wgsmat$obsvaf>=0.1,])
  
  # rens
  # to calculate the nanoseq burden, I just want the average mutation burden per cell.
  # as all the cells have the same number of mutations, it's just the number of columns minus the dummies.
  one_cell = tsim[1,]
  nsburden = length(which(substr(one_cell,1,1)=='t')) + length(which(substr(one_cell,1,1)!='t'))*muts_per_cell_division
  
  return(c(total_muts, wgscalled, nsburden))
}

tout = run_sim.fn(trunklen=as.numeric(trow['trunk_length']), 
                  cellnum=as.numeric(trow['cell_number']), 
                  bottleneck_props = as.numeric(trow['bottleneck_prop']), 
                  bottleneck_gens = as.numeric(trow['bottleneck_time']), 
                  muts_per_cell_division = as.numeric(trow['muts_per_doubling_div']))

towrite = as.character(c(as.numeric(trow), tout))
filename = paste0(paste(as.numeric(trow), collapse='_'), '_out.txt')

writeLines(towrite, filename)
