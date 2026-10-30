#this script reads the outputs of all the STRUCTURE runs and processes them
#final output is the 23.fullResults_date.Rdata file
#intermediate Rdata files 23a-d are also written in case the script fails partway

#this script requires the following files and folders in the working directory:
#00b.runParameters.txt, 10.SLURM_outputs/, 11b.true_ancestry/, 11c.SNP_files/,
# 12.STRUC_outputs/, 20a.mainLog.txt, 20b.jobMap.txt
#The numbered sections must be run in order. Section 2 takes the longest,
# as it reads through STRUCTURE output files of every single run.

######### 0.initialization #####################################################

#this section reads run parameters and sets up tools for logging and environment-cleaning
{
#time at which analysis starts, for tracking
startTime <- Sys.time()
#message-logging function
logMsg <- function(msg, time_since = startTime){
  elapsed_sec <- as.numeric(Sys.time() - time_since, units = 'secs')
  elapsed <- ifelse(elapsed_sec < 60,
                    paste0(round(elapsed_sec, 0), 's '),
                    paste0(round(elapsed_sec/60, 0), 'm '))
  cat(paste0(elapsed, msg, '\n'))
  flush.console()
}

#names of objects not to be deleted when environment is cleared at the end of each section
notDelete <- c('notDelete', 'startTime', 'logMsg',
               'pMatch', 'uMatch',
               'tRes', 'pRes', 'indData',
               'empFst', 'fData', 'selT',
               'totIndSampSize',
               'simSeeds', 'pySeeds')

#load run parameters
parameters <- read.table('00b.runParameters.txt', sep='=', header=FALSE, stringsAsFactors=FALSE, comment.char='#')
for(i in 1:nrow(parameters)){
  assign(parameters[i, 1], as.numeric(eval(parse(text=strsplit(as.character(parameters[i,2]), ',')))), envir=.GlobalEnv)
  notDelete <- c(notDelete, parameters[i, 1])
}
totIndSampSize <- RECINDSAMPSIZE * (1 + INVINDSAMPSIZERATIO)

#rows of data are grouped by keys
#pKey groups rows by unique combination of these 4 parameters:
# invasion scenario: invasion/no-invasion (invScenario)
# population structure scenario: panmictic/substructured (strucScenario)
# target Fst level (tarFst)
# sampling time following start of invasion (samSecT)
#uKey further splits pKey groups by unique combination of:
# simulation replicate (simRep)
#so within a uKey group, rows differ by these STRUCTURE analysis parameters:
# K-value (K)
# replicate run (strucRep)

#pMatch: same simulation parameters and scenarios; different simulation replicates, K, strucRep
pMatch <- c('strucScenario', 'invScenario', 'tarFst', 'samSecT')
#uMatch is the same as pMatch, but further split by simulation replicate; different K, strucRep
uMatch <- c('simRep', pMatch)

rm(list = ls()[!ls() %in% notDelete])
}

######### 1.load jobMap ########################################################

#this section loads the jobMap table, which is used as the basis for the main results table
{
logMsg('start 1.load jobMap')

tRes <- read.table('20b.jobMap.txt', header = T, 
             colClasses = c('character', 'character', 'character', 'character', 'character', #arrayID, status, jobID, parentJobID, simRep
                            'character', 'character', 'numeric', 'integer', 'integer', 'integer', #inv, struc, tarFst, samSecT, K, strucRep
                            'integer', 'character', 'character', 'integer')) #seed, tAncFP, snpFP, runTime
tRes$uKey <- do.call(paste, c(tRes[, uMatch], sep = '_'))
tRes$pKey <- do.call(paste, c(tRes[, pMatch], sep = '_'))
numCols <- c('numPolySites', 'LK', 'varLn', 'meanAlpha')
tRes[numCols] <- replicate(length(numCols), rep(NA_real_, nrow(tRes)))
listCols <- c('netNucDist', 'strucHe', 'strucFst', 'qAncTab', 'tAncTab')
tRes[listCols] <- replicate(length(listCols), vector('list', nrow(tRes)), simplify=F)

logMsg('good job 1.load jobMap')
rm(list = ls()[!ls() %in% notDelete])
}

######### 2.read output files for each jobID, calculate tRes for K=2 runs ######

#this section reads the output of every STRUCTURE run, and processes data for later analysis
#can take several hours, depending on number of runs
{
logMsg('start 2.make tRes, reading output files')

#getting number of polymorphic sites: doing it by SNP file because many rows share SNP files
for(snpFP in unique(tRes$snpFP)){
  rows <- which(tRes$snpFP == snpFP)
  tRes$numPolySites[rows] <- length(strsplit(readLines(snpFP, n=1), ' ')[[1]]) - 1
}
logMsg('got number of poly sites')
  
#getting seeds of SLiM simulation jobs and msprime mutation simulations
# from simulation log files
simSeeds <- c(); pySeeds <- c()
for(simRep in paste(sprintf("%02d", 1:NUMSIMREPS))){
  simLog <- readLines(paste0('10.SLURM_outputs/1_simRep', simRep, '.txt'))
  simSeeds <- c(simSeeds, simLog[grep('random seed:', simLog)+1])
  pySeedLine <- strsplit(grep('seed for msprime', simLog, value = T), ' ')[[1]]
  pySeeds <- c(pySeeds, pySeedLine[length(pySeedLine)])
}
logMsg('got SLiM and msprime seeds')

#get 10% intervals of number of rows in tRes, to not clog up log file
numRows_tRes <- nrow(tRes)
print_quants <- as.integer(quantile(seq(0, numRows_tRes), probs = seq(0.1, 1, 0.1)))

logMsg(paste0('starting to fill tRes ', numRows_tRes, ' rows'))
for(row in 1:nrow(tRes)){
  
  for(varName in c('arrayID', 'K')){
    assign(varName, tRes[row, varName])
  }
  
  if(row %in% print_quants){
    logMsg(paste0('Reading row ', row, '/', numRows_tRes, ', arrayID ', arrayID))
  }
  
  tryCatch({
  
  #get data from STRUCTURE full result file
  fFile <- tryCatch(
    suppressWarnings(readLines(con = paste0('12.STRUC_outputs/arrayID', arrayID, '_f.txt'))),
    error = function(e) {
      logMsg(paste0('failed to find: ', paste0('12.STRUC_outputs/arrayID', arrayID, '_f.txt'), ', trying _f version'))
      tryCatch(
        suppressWarnings(readLines(con = paste0('12.STRUC_outputs/arrayID', arrayID, '_f'))),
        error = function(e2) {
          logMsg(paste0('failed to find: ', paste0('12.STRUC_outputs/arrayID', arrayID, '_f')))
        }
      )
    }
  )
  
  tRes$LK[row] <- as.numeric(sub('.*= *', '', grep('Estimated Ln Prob of Data   = ', fFile, value = T)))
  tRes$varLn[row] <- as.numeric(sub('.*Variance of ln likelihood *= *', '', grep('Variance of ln likelihood', fFile, value = T)))
  for(clusterK in 1:K){
    tRes$strucHe[[row]] <- c(tRes$strucHe[[row]], as.numeric(sub('.*: *', '', fFile[grep('expected heterozygosity', fFile) + clusterK])))
    tRes$strucFst[[row]] <- c(tRes$strucFst[[row]], as.numeric(sub('.*= *', '', grep(paste0('Mean value of Fst_', clusterK), fFile, value=T))))
  }
  if(K > 1){
    tRes$meanAlpha[row] <- as.numeric(sub('.*= *', '', grep('alpha', fFile, value=T)))
    tRes$netNucDist[[row]] <- 
      read.table(text = fFile[(grep('Net nucleotide distance', fFile) + 4):(grep('Net nucleotide distance', fFile) + 3 + K)]
    )[, 2:(1 + K)]
  }
  
  if(K == 2){
    
    #get q-values from STRUCTURE full result file
    tab <- do.call(rbind, strsplit(
      sub('.*(p.*)', '\\1', gsub('([,():])', ' ', 
               fFile[(grep('Inferred clusters', fFile)+1):(grep('Inferred clusters', fFile) + totIndSampSize)]
               )), '\\s+'))[, -2]
    tRes$qAncTab[[row]] <- setNames(data.frame(
      tab[, 1],
      lapply(as.data.frame(tab[, 2 : (1 + K)]), as.numeric)),
      c('indID', paste0('qAnc', 1:K)))
    
    #get true ancestry of individuals
    tRes$tAncTab[[row]] <- setNames(read.table(tRes$tAncFP[row], header=F, stringsAsFactors=F), 
                                    c('indID', 'tAncInv', 'tAncRec'))
  }
  rm(fFile)
  },
  error = function(e) {
    message('Error at row ', row, '/', nrow(tRes), ', arrayID ', arrayID, ': ', conditionMessage(e))
    quit(status = 1)
  })
}

#Add isMaxLK column: for every uMatch+K, strucRep row with highest LK gets 1,
# the rest get 0. For each set of replicate STRUCTURE runs, this identifies
# the highest-LK one. These are the STRUCTURE runs from which we take Q-values
# for later analyses.
tRes$isMaxLK <- 0
for(uK in unique(tRes$uKey)){
  for(k in STRUCMINK:STRUCMAXK){
    idx <- which(tRes$uKey == uK & tRes$K == k)
    tRes$isMaxLK[idx[which.max(tRes$LK[idx])]] <- 1
  }
}

save(tRes, file = paste0('23a.tRes_', Sys.Date(), '.Rdata'))

logMsg('good job 2.make tRes, saved 23a.tRes')
rm(list = ls()[!ls() %in% notDelete])
}

######### 3.make pRes ##########################################################

#this section makes the pRes table, which has 1 row per unique combination of simulation parameters (pKey)
#then calculates the deltaK K-selection results, and identifies for which simReps deltaK selects K=2
#then gathers individual-ID results from the max-LK STRUCTURE K=2 runs at the simReps for which deltaK selected K=2
{
logMsg('start 3.make pRes, calculate deltaK, get relevant indID results')
pRes <- tRes[!duplicated(tRes$pKey), c(pMatch, 'pKey')]

#add necessary columns to pRes
listCols <- c('LK',  'sd',  'LK2', 'LK3', 'dK', 'deltaK', 'simRepsK2')
pRes[listCols] <- replicate(length(listCols), vector('list', nrow(pRes)), simplify=F)
pRes$freqK <- NA #most frequently chosen K across simReps

#calculate deltaK for each uKey-group of rows in tRes (all K and strucReps for a combination of parameter values)
for(uKey in unique(tRes$uKey)){
  for(varName in c('pKey', uMatch)){
    assign(varName, tRes[which(tRes$uKey == uKey)[1], varName])
  }
  
  redRes <- tRes[tRes$uKey == uKey, c('K', 'LK')]
  
  dK_table <- aggregate(LK ~ K, data = redRes, FUN = mean)
  dK_table$sd <- aggregate(LK ~ K, data = redRes, FUN = sd)[,2]
  for(k in (STRUCMINK+1):STRUCMAXK){
    dK_table[dK_table$K==k, 'LK2'] <- dK_table[dK_table$K==k, 'LK'] - dK_table[dK_table$K==k-1, 'LK']
  }
  for(k in (STRUCMINK+1):(STRUCMAXK-1)){
    dK_table[dK_table$K==k, 'LK3'] <- abs(dK_table[dK_table$K==k+1, 'LK2'] - dK_table[dK_table$K==k, 'LK2'])
    dK_table[dK_table$K==k, 'dK'] <- dK_table[dK_table$K==k, 'LK3'] / dK_table[dK_table$K==k, 'sd']
  }
  #most supported value of K, according to deltaK
  chosen_k <- dK_table$K[which.max(dK_table$dK)]
  
  #add data from this uKey to pRes
  pRow <- which(pRes$pKey == pKey)
  for(m in listCols[!listCols %in% c('deltaK', 'simRepsK2')]){
    pRes[[m]][[pRow]] <- rbind(pRes[[m]][[pRow]], dK_table[, m])
  }
  pRes$deltaK[[pRow]] <- c(pRes$deltaK[[pRow]], chosen_k)
  
  #if deltaK selected K=2, record this simRep number for this pRes row
  if(chosen_k == 2){pRes$simRepsK2[[pRow]] <- c(pRes$simRepsK2[[pRow]], simRep)}
}

#get most frequently chosen K across simReps
for(row in 1:nrow(pRes)){
  tab <- table(factor(pRes$deltaK[[row]], levels = (STRUCMINK+1):(STRUCMAXK-1)))
  pRes$freqK[row] <- ifelse( #else triggered by a tie
    sum(tab == max(tab)) == 1,
    as.numeric(names(which.max(tab))),
    NA
  )
}

save(pRes, file = paste0('23b.pRes_', Sys.Date(), '.Rdata'))

logMsg('good job 3.make pRes, saved 23b.pRes')
rm(list = ls()[!ls() %in% notDelete])
}

######### 4.make indData: individual-level true/inferred ancestry ##############

#this section merges data from the tRes$qAncTab and tRes$tAncTab to get, for
# each uKey, a table with one row per sampled individual, showing its true
# ancestry proportions and its Q-values from the maxLK STRUCTURE K=2 run
#indData is a list of tables, with one table per uKey
{
logMsg('start 4.make indData')

#make indData, which is required for sr/cID/heatmap results
#takes <30sec; data is already all in tRes, this just reorganizes it
indData <- list()
for(simRep in 1:NUMSIMREPS){
  for(tarFst in TARFST){
    for(samSecT in SAMSEC){
      for(invScenario in c('Allo', 'Sec')){
        for(strucScenario in c('pan', 'sub')){
          uKey <- paste(sprintf('%02d', simRep), strucScenario, invScenario, tarFst, samSecT, sep = '_')
          tRow <- which(tRes$uKey == uKey & tRes$K == 2 & tRes$isMaxLK == 1)
          logMsg(paste0('processing uKey ', uKey, ', tRow ', tRow))
          
          #merge true ancestry and Q-values into one table, matching by indID because they're in different orders
          iData <- merge(tRes$qAncTab[[tRow]], tRes$tAncTab[[tRow]], by = 'indID')
          
          #get popIDs: populations from which individuals were sampled
          iData$pop <- sub('^p(\\d+)i.*', '\\1', iData$indID)
          
          #get true identity category (either pure invader, hybrid, or pure recipient) for each individual
          iData$trueID <- ifelse(iData$tAncInv==1, 'inv', ifelse(iData$tAncRec==1, 'rec', 'hyb'))
          
          #put table into the indData list
          indData[[uKey]] <- iData
        }
      }
    }
  }
}

save(indData, file = paste0('23c.indData_', Sys.Date(), '.Rdata'))

logMsg('good job 4.make indData, saved 23c.indData')
rm(list = ls()[!ls() %in% notDelete])
}

######### 5.get success rates and correct ID rates across 0-0.5 qT range #######

#this section gets 2 results for each scenario, across the qT range of 0.01-0.49:
#The Hybrid Signal's success rate: proportion of K=2 simReps where it gives
# the correct answer, either by detecting 0 or 1+ hybrids. This is saved as
# an empirical cdf in pRes$sr_ecdf.
#The correct ID rate: the proportion of individuals correctly labelled as pure
# invader, hybrid, or pure recipient (native). This is saved as a vector of 49
# values in pRes$cID.

{
logMsg('start 5.get SR ecdf and cID across 0-0.5 qT range')

pRes$sr_ecdf <- vector('list', nrow(pRes))
pRes$cID <- vector('list', nrow(pRes))

for(pRow in seq(nrow(pRes))){
  pKey <- pRes$pKey[[pRow]]
  logMsg(paste0('Processing pKey ', pKey, ', row ', pRow))
  
  #simReps that got K=2
  k2sr <- pRes$simRepsK2[[pRow]]
  
  ##getting the SR empirical cdf
  #get all the K=2 uKeys
  uKeys <- paste(k2sr, pKey, sep = '_')
  #what is the closest q-value to 0.5 in each simRep?
  # this is the individual that requires the highest Q-threshold to be labelled a pure of either pop
  maxQ <- unlist(sapply(uKeys, function(n) 0.5 - min(abs(0.5 - indData[[n]]$qAnc1))))
  #get ecdf
  pRes$sr_ecdf[[pRow]] <- ecdf(maxQ)
  
  ##getting the cID
  cIDTab <- matrix(NA, nrow = length(k2sr), ncol = 49)
  for(simRep_pos in seq_along(k2sr)){
    simRep <- k2sr[simRep_pos]
    iData <- indData[[paste0(simRep, '_', pKey)]]
    
    #true ancestry categories
    tAnc <- iData$trueID
    
    #determine which cluster most likely corresponds to inv
    invCl <- which.max(c(cor(iData$qAnc1, iData$tAncInv), cor(iData$qAnc2, iData$tAncInv)))
    #get those q-values
    qAncInv <- iData[[paste0('qAnc', invCl)]]
    
    for(qT_pos in seq_along(seq(0.01, 0.49, by = 0.01))){
      qT <- seq(0.01, 0.49, by = 0.01)[qT_pos]
      
      #get category of each individual under this q-threshold
      qAncCat <- ifelse(qAncInv < qT, 'rec', ifelse(qAncInv > 1 - qT, 'inv', 'hyb'))
      #get success rate
      cIDTab[simRep_pos, qT_pos] <- mean(qAncCat == tAnc)
    }
  }
  #get mean success rate across the K=2 selected simReps, at each q-threshold
  pRes$cID[[pRow]] <- colMeans(cIDTab)
}

save(pRes, file = paste0('23b.pRes_', Sys.Date(), '.Rdata'))

logMsg('good job 5.get SR ecdf and cID across 0-0.5 qT range, saved new pRes')
rm(list = ls()[!ls() %in% notDelete])
}

######### 6.pre-invasion Fst ##################################################

#This section gets Fst data, as calculated from the pre-invasion tree-sequences.
#Note: different from the "empirical-POV" Fst calculated later from the SNP data.
#Also gets the Fst-selected divergence times at which invasion was started.
#Results in fData and selT, which are saved with the fullResults file in section 8.
{
logMsg('start 6.pre-invasion Fst')

fData <- read.table('21b.Fst.txt', row.names = 1)
fRef <- as.data.frame(do.call(rbind, strsplit(rownames(fData), '_')))
colnames(fRef) <- c('simRep', 'metric', 'pop')

#getting Fst-selected time for each simRep/tarFst/strucScenario combination
sdata <- fData[fRef$metric == 'hFst', ]
selT <- expand.grid(tarFst = TARFST, strucScenario = c('pan', 'sub'))
selT$t <- vector('list', nrow(selT))
for(tarFst in TARFST){
  for(strucScenario in c('pan', 'sub')){
    dat <- c()
    for(simRep in paste(sprintf("%02d", 1:NUMSIMREPS))){
      row <- which(rownames(sdata) == paste(simRep, 'hFst', strucScenario, sep = '_'))
      dat <- c(dat, as.numeric(SAMALLO * (which.min(abs(sdata[row, ] - tarFst)) - 1)))
    }
    selT$t[[which(selT$tarFst == tarFst & selT$strucScenario == strucScenario)]] <- dat
  }
}

#get mean of each metric-pop combo
for(mpC in unique(paste0(fRef$metric, '_', fRef$pop))){
  idx <- which(paste0(fRef$metric, '_', fRef$pop) == mpC)
  fData[paste0('mean_', mpC), ] <- colMeans(fData[idx, ])
  fRef <- rbind(fRef, data.frame(
    simRep = 'mean', metric = fRef$metric[idx[1]], pop = fRef$pop[idx[1]]
  ))
}

logMsg('good job 6.pre-invasion Fst')
rm(list = ls()[!ls() %in% notDelete])
}

######### 7.empirical-POV Fst ##################################################

#This section calculates Fst from an empirical study's point of view, using the
# population from which individuals were sampled as groups.
# It's slow because it has to read each SNP file.
# Results in empFst, which is saved to 23d.empFst here, and later saved with the
# fullResults file in section 8.
{
logMsg('start 7.get empirical-POV Fst')

#Fst with hierfstat
library(hierfstat)
empFst <- pRes[, c(pMatch, 'pKey')]
empFst[paste(1:NUMSIMREPS)] <- replicate(NUMSIMREPS, rep(NA_real_, nrow(empFst)))

#get data across all scenarios
for(simRep in 1:NUMSIMREPS){
  for(tarFst in TARFST){
    for(samSecT in SAMSEC){
      for(invScenario in c('Allo', 'Sec')){
        for(strucScenario in c('pan', 'sub')){
          
          ### getting the SNP file ###
          uKey <- paste(sprintf('%02d', simRep), strucScenario, invScenario, tarFst, samSecT, sep = '_')
          logMsg(paste('processing', uKey))
          tRow <- which(tRes$uKey == uKey & tRes$K == 2 & tRes$isMaxLK == 1)
          snp <- read.table(tRes$snpFP[tRow])
          
          
          ### Fst between 'sampling-site groups' from hierfstat ###
          pRow <- which(empFst$pKey == paste(strucScenario, invScenario, tarFst, samSecT, sep = '_'))
          snp[, 1] <- ifelse(grepl('p1i', snp[, 1]), 1, 2) #hierfstat requires popIDs in column 1
          
          #my data has 1 row per haploid; have to combine it so each row is 1 diploid
          #the two alleles at each locus for each individual have to be collated (e.g. instead of 2 and 1, 21)
          snp[snp == 0] <- 2 #replace all 0 with 2, so the next step (addition to collate alleles into genotype) works
          pops <- snp[seq(1, nrow(snp), by = 2), 1] #popID vector per diploid individual
          #allele matrices with each diploid individual split
          hap1 <- snp[seq(1, nrow(snp), by = 2), -1]
          hap2 <- snp[seq(2, nrow(snp), by = 2), -1]
          #add alleles to get genotypes
          geno <- (hap1 * 10) + hap2
          #diploid genotype dataset
          dipsnp <- data.frame(pops, geno)
          
          empFst[pRow, paste(simRep)] <- pairwise.WCfst(dipsnp)[1, 2]
        }
      }
    }
  }
}

save(empFst, file = paste0('23d.empFst_', Sys.Date(), '.Rdata'))

logMsg('good job 7.get empirical-POV Fst, saved 23d.empFst')
rm(list = ls()[!ls() %in% notDelete])
}

######### 8.save fullResults file ##############################################

logMsg('start 8.save 23.fullResults')
save(tRes, pRes, indData, empFst, fData, selT, simSeeds, pySeeds,
     file = paste0('23.fullResults_', Sys.Date(), '.Rdata'))
logMsg('good job 8.save 23.fullResults')
