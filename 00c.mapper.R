
parameters <- read.table('00b.runParameters.txt', sep='=', header=FALSE, stringsAsFactors=FALSE, comment.char='#')
for(i in 1:nrow(parameters)){
  assign(parameters[i, 1], as.numeric(eval(parse(text=strsplit(as.character(parameters[i,2]), ',')))), envir=.GlobalEnv)
}

################################################################################
#############################   make jobMap   ###############################

#all jobs to be used for deltaK
jobMap <- expand.grid(
  strucRep = 1:STRUCREPS,
  K = STRUCMINK:STRUCMAXK,
  samSecT = SAMSEC,
  tarFst = TARFST,
  strucScenario = c('pan', 'sub'),
  invScenario = c('Allo', 'Sec'), 
  simRep = sprintf('%02d', 1:NUMSIMREPS), 
  stringsAsFactors = FALSE
)

jobMap$arrayID <- sprintf('%04d', seq_len(nrow(jobMap)))

#generate a random number to use as base for STRUCTURE job seeds
seed_base <- sample(1:1000000000, 1)
cat(paste0('Seed base for STRUCTURE jobs: ', seed_base))
jobMap$seed <- seed_base + as.integer(jobMap$arrayID)

jobMap$tAncFP <- paste0('11b.true_ancestry/tAnc_', jobMap$strucScenario, 
                        jobMap$invScenario, '_simRep', jobMap$simRep, '_tarF', 
                        jobMap$tarFst, '_samT', sprintf('%02d', jobMap$samSecT), '.txt')
jobMap$snpFP <- paste0('11c.SNP_files/SNP_', jobMap$strucScenario, jobMap$invScenario, 
                       '_simRep',  jobMap$simRep, '_tarF', jobMap$tarFst, 
                       '_samT', sprintf('%02d', jobMap$samSecT), '.txt')
jobMap$status <- 'notSched'
jobMap$jobID <- NA
jobMap$parentJobID <- NA
jobMap$runTime <- NA

#reorder, this is important because the columns are read in a specific order
jobMap <- jobMap[, c(
  'arrayID', 'status', 'jobID', 'parentJobID', 'simRep', 'invScenario', 'strucScenario', 
  'tarFst', 'samSecT', 'K', 'strucRep', 'seed', 'tAncFP', 'snpFP', 'runTime'
)]

write.table(jobMap, '20b.jobMap.txt', sep = '\t', quote = FALSE, row.names = FALSE)

################################################################################
############################   make popIDs table   #############################

options(scipen = 999)
numSets <- 1 + max(ceiling(SAMSEC/SAMALLO))
sampTimes <- seq(TBUR, TBUR+TALLO, by = SAMALLO)
popIDs <- data.frame(
  tick = sampTimes,
  pan = rep(5 + (1:numSets), length.out = length(sampTimes)),
  sub3 = rep(5 + numSets + (1:numSets), length.out = length(sampTimes)),
  sub4 = rep(5 + (2*numSets) + (1:numSets), length.out = length(sampTimes)),
  sub5 = rep(5 + (3*numSets) + (1:numSets), length.out = length(sampTimes))
)
write.table(popIDs, '21a.popIDs.txt', sep = '\t', quote = FALSE, row.names = FALSE)

################################################################################
#######################   modify mainParams if needed   ########################

lines <- readLines('02b.mainParams.txt')
update_needed <- FALSE

bur_idx <- grep('#define BURNIN', lines)
if(strsplit(lines[bur_idx], ' ')[[1]][3] != as.character(NUMBURNIN)){
  update_needed <- TRUE
  bur_rep <- strsplit(lines[bur_idx], ' ')[[1]]
  bur_rep[3] <- as.character(NUMBURNIN)
  lines[bur_idx] <- paste(bur_rep, collapse = ' ')
}

ite_idx <- grep('#define NUMREPS', lines)
if(strsplit(lines[ite_idx], ' ')[[1]][3] != as.character(NUMITERATIONS)){
  update_needed <- TRUE
  ite_rep <- strsplit(lines[ite_idx], ' ')[[1]]
  ite_rep[3] <- as.character(NUMITERATIONS)
  lines[ite_idx] <- paste(ite_rep, collapse = ' ')
}

if(update_needed) {
  writeLines(lines, '02b.mainParams.txt')
  cat('UPDATED mainparams.txt\n')
}
