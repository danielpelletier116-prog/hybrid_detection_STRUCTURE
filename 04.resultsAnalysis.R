#This script analyses the output of the simulations and 03.STRUCTUREanalysis.R
#Requires only 2 files in the working directory:
# 23.fullResults_"date".Rdata
# 00b.runParameters.txt
#It outputs main plots in 13a.plots, QvQ plots in 13b.QvQplots.
#The 0.initialization section must be run first, then any section can be
# run individually to produce plots and output results.
#After running section 0, to clean environment between sections, use:
# rm(list = ls()[!ls() %in% notDelete])

######### 0.initialization #####################################################

#this section reads run parameters and data files, and sets up tools for logging and environment-cleaning
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
               'pMatch',
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

#read the latest-date version of 23.fullResults in the working directory
logMsg('load 23.fullResults for 04.resultsAnalysis.R')
load(grep('23.fullResults', list.files(), value = T)[length(grep('23.fullResults', list.files()))])
logMsg('good job loading 23.fullResults')

rm(list = ls()[!ls() %in% notDelete])
}

######### 1.plot Fst ###########################################################

#plots Fig S1, and output Fst-related results for text
{
logMsg('start 1.plot Fst')
dir.create('13a.plots', showWarnings = FALSE)

#make fRef: reference table with info about each row of fData
fRef <- as.data.frame(do.call(rbind, strsplit(rownames(fData), '_')))
colnames(fRef) <- c('simRep', 'metric', 'pop')

#Fig S1: Fst and selT plot, with empirical-POV Fst and Fst between subpops
# pre-invasion Fst, and selT, calculated in 3.R section 6
# empFst calculated in 3.R section 7, wcFst between 'sample-site groups' with hierfstat
{
#get max selT, for x-range of plot
max_selT <- max(unlist(selT$t))

#plotting colours
cols <- c(
  hFst_pan = '#5B8FA8', hFst_sub = '#B1746F', demeHFst_sub = '#3DA873'
)

png(filename = paste0('13a.plots/s1.Fst.png'), width = 1600, height = 800)

par(mar = c(4.1, 4.1, 0.5, 0.5), cex = 2.2)

plotColsBase <- cols[paste(fRef$metric, fRef$pop, sep = '_')]
plotCols <- ifelse(fRef$simRep == 'mean', plotColsBase, adjustcolor(plotColsBase, alpha.f = 0.4))
#plot per-simRep and mean Fst accumulation
matplot(seq(0, max_selT + 3*SAMALLO, by = SAMALLO), t(fData[, 1:(4+(max_selT/SAMALLO))]), 
        type = 'l', xlab = 'Generations of divergence', ylab = 'Fst',
        lwd = ifelse(fRef$simRep == 'mean', 4, 1), lty = 1, bty = 'n', col = plotCols,
        xaxs = 'i', yaxs = 'i')

#add tarFst horizontal lines
abline(h = TARFST, lwd = 2, lty = 2, col = 'grey50')

#selT frequency histograms
#first, get max bin count across all tarFst
maxCount <- 0
for(tarFst in TARFST){
  pT <- selT$t[[which(selT$tarFst == tarFst & selT$strucScenario == 'pan')]]
  sT <- selT$t[[which(selT$tarFst == tarFst & selT$strucScenario == 'sub')]]
  xrange <- seq(from = min(pT, sT) - 5, to = max(pT, sT) + 5, by = 10)
  
  maxCount <- max(
    maxCount,
    max(hist(pT, breaks = xrange, plot = F)$counts),
    max(hist(sT, breaks = xrange, plot = F)$counts)
  )
}

#plot hists
for(tarFst in TARFST){
  pT <- selT$t[[which(selT$tarFst == tarFst & selT$strucScenario == 'pan')]]
  sT <- selT$t[[which(selT$tarFst == tarFst & selT$strucScenario == 'sub')]]
  xrange <- seq(from = min(pT, sT) - 5, to = max(pT, sT) + 5, by = 10)
  
  pH <- hist(pT, breaks = xrange, plot = F)
  sH <- hist(sT, breaks = xrange, plot = F)
  
  rect(pH$mids - 5, tarFst, pH$mids + 5, tarFst + (0.035*pH$counts/maxCount),
       col = adjustcolor(cols['hFst_pan'], alpha.f = 0.6))
  rect(sH$mids - 5, tarFst, sH$mids + 5, tarFst + (0.035*sH$counts/maxCount),
       col = adjustcolor(cols['hFst_sub'], alpha.f = 0.6))
}

#add empFst measures (empirical-POV Fst)
#for each tarFst, at each samSecT (added to mean selT), 4 points with mean empFst across 50 simReps
#pch 1 for Allo, 13 for Sec, colours for pan/sub
for(strucScenario in c('pan', 'sub')){
  for(tarFst in TARFST){
    #get average selT for this strucScenario/tarFst combo
    avSelT <- mean(selT$t[[which(selT$tarFst == tarFst & selT$strucScenario == strucScenario)]])
    
    for(invScenario in c('Allo', 'Sec')){
      for(samSecT in SAMSEC){
        pKey <- paste(strucScenario, invScenario, tarFst, samSecT, sep = '_')
        avEmpFst <- mean(as.numeric(empFst[empFst$pKey == pKey, 6:(5+NUMSIMREPS)]))
        
        points(
          avSelT + samSecT, avEmpFst, 
          pch = ifelse(invScenario == 'Allo', 1, 13), 
          col = cols[paste0('hFst_', strucScenario)], cex = 2, lwd = 3
          )
        
      }
    }
  }
}

dev.off()
}

#output mean selT (times to get to target Fst levels) for text
{
cat('Mean times to reach each target Fst level, for each pop structure scenario\n')
for(strucScenario in c('pan', 'sub')){
  cat(strucScenario, sep = '\n')
  for(tarFst in TARFST){
    meanT <- mean(selT$t[[which(selT$tarFst == tarFst & selT$strucScenario == strucScenario)]])
    cat(tarFst, meanT, 'gens\n')
  }
}
}

rm(list = ls()[!ls() %in% notDelete])
logMsg('good job 1.plot Fst')
}

######### 2.plot deltaK frequency ##############################################

#plots Fig S2: frequency of deltaK K-selections, for each scenario
{
logMsg('start 2.plot deltaK frequency')
dir.create('13a.plots', showWarnings = FALSE)

#Fig S2: frequency of deltaK K-selections, for each scenario
{
strucCols <- c(pan = '#5B8FA8', sub = '#B1746F')
prettyNames <- c(Allo = 'No-invasion', Sec = 'Invasion', pan = 'Panmictic', sub = 'Substructured')
invLty <- c(Allo = 2, Sec = 1)

#optional: show counts of each bar
show_counts <- F

#laymat and dimensions
{
#one row per tarFst, one col per samSecT/strucScenario/invScenario
layMat <- matrix(0, nrow = 3 + length(TARFST), ncol = 2 + (4*length(SAMSEC)))
layMat[1, 3:ncol(layMat)] <- max(layMat) + rep(c(1,2), each = 2*length(SAMSEC)) #strucScenario column labels
layMat[2, 3:ncol(layMat)] <- max(layMat) + rep(c(1:(2*length(SAMSEC))), each = 2) #samSecT column labels
layMat[2 + seq_along(TARFST), 1] <- max(layMat) + seq_along(TARFST) #tarFst row labels
layMat[3:(nrow(layMat)-1), 2] <- max(layMat) + 1 #first y-axis label
layMat[nrow(layMat), 3:ncol(layMat)] <- max(layMat) + 1 #x-axis label
layMat[nrow(layMat), 1:2] <- max(layMat) + 1 #invScenario lty legend
layMat[3:(nrow(layMat)-1), 3:ncol(layMat)] <- max(layMat) + 
  matrix(1:(length(SAMSEC)*length(TARFST)*4), nrow = length(TARFST), ncol = length(SAMSEC)*4, byrow = F)

widths = c(1.25, 1.5, rep(3.5, 4*length(SAMSEC)))
heights = c(1, 1, rep(3.5, length(TARFST)), 1.5)
}

png(filename = '13a.plots/s2.deltaK_freq.png', width = sum(widths)*250, height = sum(heights)*250)
layout(layMat, widths = widths, heights = heights)
par(mar = c(0,0,0,0), cex = 8, xpd = NA, mgp = c(3, 0.5, 0))

#plotting labels and legends
{
for(strucScenario in c('pan', 'sub')){
  plot.new()
  text(0.5, 0.7, cex = 1.2, prettyNames[strucScenario]); lines(c(0.05, 0.95), c(0.3, 0.3), lwd = 5)
}
for(samSecT in rep(SAMSEC, 2)){ #samSecT column labels
  plot.new()
  text(0.5, 0.7, cex = 1.2, paste(samSecT, 'gens.')); lines(c(0.05, 0.95), c(0.3, 0.3), lwd = 5)
}
for(i in seq_along(TARFST)){ #tarFst row labels
  plot.new(); text(0.75, 0.5, srt = 90, cex = 1.2, TARFST[i])
  if(i == ceiling(length(TARFST) / 2)){
    text(0.25, ifelse(length(TARFST)%%2==1, 0.5, 0.0), 
         'Target Fst level', srt = 90, font = 2, cex = 1.2)
  }
}
plot.new() #y-axis label
text(0.5, 0.5, srt = 90, cex = 1.2, 'Frequency of K-selection across simulation replicates')
plot.new() #x-axis label
text(0.5, 0.35, cex = 1.2, 'Most supported number of populations by ΔK')
plot.new() #invScenario lty legend
legend('bottomleft', bty = 'n', legend = prettyNames[names(invLty)], 
       lty = invLty, lwd = 8, col = 'black', ncol = 2, cex = 1.2
)
}

#deltaK frequency plots
freq_y_max <- max(sapply(pRes$deltaK, function(x) {
  tab <- table(x)
  max(tab)}))
for(strucScenario in c('pan', 'sub')){
  for(samSecT in SAMSEC){
    for(invScenario in c('Allo', 'Sec')){
      par(mar = c(0.5, ifelse(invScenario == 'Allo', 0.6, 0.3), 0.5, ifelse(invScenario == 'Allo', 0.3, 0.6)))
      for(tarFst in TARFST){
        row <- which(pRes$pKey == paste(strucScenario, invScenario, tarFst, samSecT, sep = '_'))
        
        cols <- rep('gray', length((STRUCMINK+1):(STRUCMAXK-1)))
        cols[((STRUCMINK+1):(STRUCMAXK-1)) == pRes$freqK[row]] <- strucCols[strucScenario]
        
        counts <- table(factor(pRes$deltaK[[row]], levels = (STRUCMINK+1):(STRUCMAXK-1)))
        bar_centers <- barplot(
          counts, ylim = c(0, freq_y_max+0.5),
          xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', col = cols, border = NA)
        
        #if option was selected, show count on each bar
        if(show_counts == TRUE){
          text(bar_centers[counts != 0], counts[counts != 0], counts[counts != 0])
        }
        
        #add box around plot
        box(lwd = 10, col = strucCols[strucScenario], lty = invLty[invScenario])
        
        #add x-axis
        if(tarFst == max(TARFST)){axis(1, at = bar_centers, labels = (STRUCMINK+1):(STRUCMAXK-1), pos = 0, tick = F)}
        #add y-axis
        if(samSecT == SAMSEC[1] & strucScenario == 'pan' & invScenario == 'Allo'){
          axis(2, at = c(0, ceiling(freq_y_max/2), freq_y_max), labels = c(0, ceiling(freq_y_max/2), freq_y_max), las = 1)
        }
      }
    }
  }
}
dev.off()
}

rm(list = ls()[!ls() %in% notDelete])
logMsg('good job 2.plot deltaK frequency')
}

######### 3.plot success rate, optimal Q-thresholds, and correct ID rate #######

#This section outputs results of the main tests, as plots and data for tables.
#Success rate (SR) of Hybrid Signal, over a range of Q-thresholds:
# whether STRUCTURE gives 0 hybrids in Allo, and >0 hybrids in Sec.
#Optimal Q-thresholds (optQT): Q-thresholds at which SR is maximized.
#Correct ID rate (cID): proportion of individuals correctly labelled.
#Using results from K=2 simulation replicates.
{
logMsg('start 3.main results section')
dir.create('13a.plots', showWarnings = FALSE)

##Fig S3: allo/sec success rate over Q-threshold range, pan only
{
#plotting stuff
{
invLty <- c(Allo = 2, Sec = 1)
prettyNames <- c(Allo = 'No-invasion', Sec = 'Invasion')
#layout matrix
qgtLm <- matrix(0, nrow = 2 + length(TARFST), ncol = 2 + length(SAMSEC))
qgtLm[1, 2 + seq_along(SAMSEC)] <- max(qgtLm) + seq_along(SAMSEC) #samSecT column labels
qgtLm[1 + seq_along(TARFST), 1] <- max(qgtLm) + seq_along(TARFST) #tarFst row labels
qgtLm[1 + seq_along(TARFST), 2] <- max(qgtLm) + 1 #y-axis label
qgtLm[nrow(qgtLm), 2 + seq_along(SAMSEC)] <- max(qgtLm) + 1 #x-axis label
qgtLm[nrow(qgtLm), c(1,2)] <- max(qgtLm) + 1 #invScenario line type legend
qgtLm[1 + seq_along(TARFST), 2 + seq_along(SAMSEC)] <- max(qgtLm) + #data
  matrix(1:(length(SAMSEC)*length(TARFST)), nrow = length(TARFST), ncol = length(SAMSEC), byrow = F)
#panel dimensions
qgtW <- c(1.5, 2, rep(6, length(SAMSEC)))
qgtH <- c(1, rep(3.5, length(TARFST)), 1.5)
}

#make plot
{
png(filename = paste0('13a.plots/s3.SR_pan.png'), 
    width = sum(qgtW)*100, height = sum(qgtH)*100)
layout(qgtLm, widths = qgtW, heights = qgtH)
par(mar = c(0.2, 0, 0.2, 0), cex = 2.5, xpd = NA)

#labels etc.
{
for(samSecT in SAMSEC){ #samSecT column labels
  plot.new()
  text(0.5, 0.6, cex = 1.2, paste(samSecT, 'gens.')); lines(c(0.05, 0.95), c(0.2, 0.2), lwd = 3)
}
for(i in seq_along(TARFST)){ #tarFst row labels
  plot.new(); text(0.7, 0.5, cex = 1.2, srt = 90, TARFST[i])
  lines(c(0.95, 0.95), c(0.05, 0.95), lwd = 3)
  if(i == ceiling(length(TARFST) / 2)){
    text(0.2, ifelse(length(TARFST)%%2==1, 0.5, 0.0), 
         'Target Fst level', srt = 90, font = 2, cex = 1.2)
  }
}
plot.new() #y-axis label
text(0.35, 0.5, cex = 1.2, srt = 90, 'Proportion of replicates successful')
plot.new() #x-axis label
text(0.5, 0.4, cex = 1.2, 'Q-threshold')
plot.new() #invScenario line type legend
legend('center', bty = 'n', legend = prettyNames[names(invLty)], 
       lty = invLty, lwd = 3, col = 'black'
)
}

par(mar = c(0.8, 0.6, 0.6, 0.2))
strucScenario <- 'pan'
#plot data
for(samSecT in SAMSEC){
  for(tarFst in TARFST){
    plot(NA, xlim = c(0, 0.6), ylim = c(0, 1.2), xaxt = 'n', yaxt = 'n', bty = 'n', xlab = '', ylab = '')
    axis(2, at = c(0, 0.5, 1), las = 1, labels = (samSecT == SAMSEC[1]), lwd = 2)
    axis(1, at = seq(0, 0.5, by = 0.1), labels = (tarFst == max(TARFST)), lwd = 2)

    for(invScenario in c('Allo', 'Sec')){
      pRow <- which(pRes$pKey == paste(strucScenario, invScenario, tarFst, samSecT, sep = '_'))
      
      #plot curve
      ec <- pRes$sr_ecdf[[pRow]]
      if(invScenario == 'Allo'){
        curve(ec(x), from = 0, to = 0.5, lty = invLty[invScenario], lwd = 3.5, add = T)
        sr_Allo <- ec(seq(0, 0.5, by = 0.01))
      } else {
        curve(1 - ec(x), from = 0, to = 0.5, lty = invLty[invScenario], lwd = 3.5, add = T)
        sr_Sec <- 1 - ec(seq(0, 0.5, by = 0.01))
      }
      
      #display number of simReps that got K=2
      text(0.55, ifelse(invScenario == 'Allo', 1, 0.05), length(pRes$simRepsK2[[pRow]]))
    }
    
    #across both invScenarios, find the qT-range where overall mean success rate is maximized
    comboSR <- (sr_Allo + sr_Sec) / 2
    rangeMaxComboSR <- seq(0, 0.5, by = 0.01)[range(which(abs(max(comboSR) - comboSR) < 0.00001))]
    
    #display max mean success rate, and the optimal qT-range
    text(mean(rangeMaxComboSR), 1.2, round(max(comboSR), digits = 2))
    points(mean(rangeMaxComboSR), max(comboSR), pch = 3, cex = 2, lwd = 4)
    arrows(x0 = rangeMaxComboSR[1], y0 = 1.1, 
           x1 = ifelse(rangeMaxComboSR[1] == rangeMaxComboSR[2], rangeMaxComboSR[1] + 0.001, rangeMaxComboSR[2]), 
           angle = 90, code = 3, lwd = 2, length = 0.15)
    text(rangeMaxComboSR[1] - 0.01, 1.1, 
         rangeMaxComboSR[1], adj = 1)
    if(rangeMaxComboSR[1] != rangeMaxComboSR[2]){
      text(rangeMaxComboSR[2] + 0.01, 1.1, 
           rangeMaxComboSR[2], adj = 0)
    }
  }
}
dev.off()
}
}

##Fig 2: 2-panel allo/sec success rate over Q-threshold range, pan only
{
#plotting stuff
{
invLty <- c(Allo = 2, Sec = 1)
prettyNames <- c(Allo = 'No-invasion', Sec = 'Invasion')
#layout matrix
qgtPLm <- matrix(0, nrow = 2, ncol = 3)
qgtPLm[1, 1] <- max(qgtPLm) + 1 #y-axis label
qgtPLm[2, 2:3] <- max(qgtPLm) + 1 #x-axis label
qgtPLm[2, 1] <- max(qgtPLm) + 1 #line type legend
qgtPLm[1, 2:3] <- max(qgtPLm) + c(1, 2) #panels
#panel dimensions
qgtPW <- c(1, 2.5, 2.5)
qgtPH <- c(2, 0.6)
}

#make plot
{
png(filename = paste0('13a.plots/2.SR.png'), 
    width = sum(qgtPW)*250, height = sum(qgtPH)*250)
layout(qgtPLm, widths = qgtPW, heights = qgtPH)

par(mar = c(0.2, 0, 0.2, 0), cex = 2.5, xpd = NA)
#labels etc.
{
plot.new() #y-axis label
text(0.35, 0.5, cex = 1.2, srt = 90, 'Proportion of replicates\nsuccessful')
plot.new() #x-axis label
text(0.5, 0.4, cex = 1.2, 'Q-threshold')
plot.new() #invScenario line type and strucScenario colour legend
legend('bottomleft', bty = 'n', legend = prettyNames, lty = invLty, lwd = 3)
}

par(mar = c(0.8, 0.2, 0.6, 0.6))
#plot data
samSecT <- 10
strucScenario <- 'pan'
for(tarFst in c(0.01, 0.1)){
  plot(NA, xlim = c(0, 0.5), ylim = c(0, 1.2), xaxt = 'n', yaxt = 'n', bty = 'n', xlab = '', ylab = '')
  axis(2, at = c(0, 0.5, 1), las = 1, lwd = 2, labels = tarFst == 0.01)
  axis(1, at = seq(0, 0.5, by = 0.1), lwd = 2)
  
  for(invScenario in c('Allo', 'Sec')){
    pRow <- which(pRes$pKey == paste(strucScenario, invScenario, tarFst, samSecT, sep = '_'))
    
    #plot curve
    ec <- pRes$sr_ecdf[[pRow]]
    if(invScenario == 'Allo'){
      curve(ec(x), from = 0, to = 0.5, lty = invLty[invScenario], lwd = 3.5, add = T)
      sr_Allo <- ec(seq(0, 0.5, by = 0.01))
    } else {
      curve(1 - ec(x), from = 0, to = 0.5, lty = invLty[invScenario], lwd = 3.5, add = T)
      sr_Sec <- 1 - ec(seq(0, 0.5, by = 0.01))
    }
  }
  
  #across both invScenarios, find the qT-range where overall mean success rate is maximized
  comboSR <- (sr_Allo + sr_Sec) / 2
  rangeMaxComboSR <- seq(0, 0.5, by = 0.01)[range(which(abs(max(comboSR) - comboSR) < 0.00001))]
  
  #display max mean success rate, and the optimal qT-range
  text(mean(rangeMaxComboSR), 1.2, round(max(comboSR), digits = 2))
  points(mean(rangeMaxComboSR), max(comboSR), pch = 3, cex = 2, lwd = 4)
  arrows(x0 = rangeMaxComboSR[1], y0 = 1.1, 
         x1 = ifelse(rangeMaxComboSR[1] == rangeMaxComboSR[2], rangeMaxComboSR[1] + 0.001, rangeMaxComboSR[2]), 
         angle = 90, code = 3, lwd = 2, length = 0.15)
  text(rangeMaxComboSR[1] - 0.01, 1.1, 
       rangeMaxComboSR[1], adj = 1)
  if(rangeMaxComboSR[1] != rangeMaxComboSR[2]){
    text(rangeMaxComboSR[2] + 0.01, 1.1, 
         rangeMaxComboSR[2], adj = 0)
  }
  
  #add in-plot Fst label
  text(0.03, 0.8, adj = 0, paste0('Fst = ', tarFst))
}
dev.off()
}
}

##Tables S1/S2/S3 data:
#Table S1: pan-only optQT ranges, SR and cID at selected optQT, K2 rep sample sizes
#Table S2: same as ts1 but sub-only, still evaluated at pan-selected consensus optQT though
#Table S3: sub-only optQT ranges, SR at those selected optQT, K2 rep sample sizes
{
#part 1: K2 sample sizes, optQT ranges, and SR (evaluated at pan- or sub-selected optQT)
cat('\n\nOptimal Q-threshold ranges, maximal success rates, and K=2 sample sizes',
    '\nIn substructured scenario, using sub-selected optQT.')
for(strucScenario in c('pan', 'sub')){
  cat('\n\nstrucScenario:', strucScenario, '\n')
  
  for(tarFst in TARFST){
    for(samSecT in SAMSEC){
      
      k2sr <- list(); ecdf <- list()
      
      for(invScenario in c('Allo', 'Sec')){
        pRow <- which(pRes$pKey == paste(strucScenario, invScenario, tarFst, samSecT, sep = '_'))
        
        #which simReps got K=2?
        k2sr[[invScenario]] <- pRes$simRepsK2[[pRow]]
        
        #get the ecdf
        ecdf[[invScenario]] <- pRes$sr_ecdf[[pRow]]
      }
      
      #across both invScenarios, find the qT-range where overall mean success rate is maximized
      comboSR <- (ecdf$Allo(seq(0, 0.5, by = 0.01)) + 1 - ecdf$Sec(seq(0, 0.5, by = 0.01))) / 2
      optQT <- seq(0, 0.5, by = 0.01)[range(which(abs(max(comboSR) - comboSR) < 0.00001))]
      
      
      ##print stuff
      cat(
        '\ntarFst:', tarFst, ', samSecT:', samSecT,
        #optimal Q-threshold range
        '\noptQT range:', paste0(optQT[1], ifelse(optQT[1] == optQT[2], '', paste0('-', optQT[2]))),
        #SR and K2 sample size of each invScenario
        '\nAllo:', ecdf$Allo(optQT[1]), '(', length(k2sr$Allo), ')',
        '\nSec:', 1 - ecdf$Sec(optQT[1]), '(', length(k2sr$Sec), ')',
        '\nAv:', max(comboSR),
        '\n')
    }
  }
}

#pan-selected optQT (see paper for explanation)
sel_optQT <- c('10' = 0.35, '30' = 0.3, 'rest' = 0.1)

#part 2: substructure-scenario success rates when using pan-selected optQT
strucScenario <- 'sub'
cat('\n\nSubstructured scenario only: success rates at pan-selected optQT\n')
for(tarFst in TARFST){
  for(samSecT in SAMSEC){
    #get the selected optimal qT
    optQT <- ifelse(tarFst == 0.01, sel_optQT[paste(samSecT)], sel_optQT['rest'])
    
    apRow <- which(pRes$pKey == paste(strucScenario, 'Allo', tarFst, samSecT, sep = '_'))
    spRow <- which(pRes$pKey == paste(strucScenario, 'Sec', tarFst, samSecT, sep = '_'))
    
    cat('\ntarFst:', tarFst, ', samSecT:', samSecT, ', optQT:', optQT,
        '\nAllo:', pRes$sr_ecdf[[apRow]](optQT), '(', length(pRes$simRepsK2[[apRow]]), ')',
        '\nSec:', 1 - pRes$sr_ecdf[[spRow]](optQT), '(', length(pRes$simRepsK2[[spRow]]), ')',
        '\nAv:', (pRes$sr_ecdf[[apRow]](optQT) + 1 - pRes$sr_ecdf[[spRow]](optQT)) / 2, '\n'
        )
  }
}

#part 3: using pan-selected optQT: cID rates (both pan and sub)
cat('\n\nCorrect ID rates, using pan-selected optQT\n')
for(strucScenario in c('pan', 'sub')){
  cat('\nstrucScenario:', strucScenario, '\n')
  
  for(tarFst in TARFST){
    for(samSecT in SAMSEC){
      optQT <- ifelse(tarFst == 0.01, sel_optQT[paste(samSecT)], sel_optQT['rest'])
      
      cat('tarFst:', tarFst, ', samSecT:', samSecT, ', optQT:', optQT, '\n')
      
      for(invScenario in c('Allo', 'Sec')){
        #get cID from pRes
        cID <- pRes$cID[[which(pRes$pKey == paste(strucScenario, invScenario, tarFst, samSecT, sep = '_'))]]
        
        cat(strucScenario, invScenario, cID[as.integer(optQT*100)], '\n')
      }
    }
  }
}
}

##Fig 3: pan-only optimal qT ranges, showing overlap range and cases where sr=1 in bold
{
optQT_ranges <- list() #for storing the optimal qT range of each scenario
max_sr <- list() #for storing max success rate of each scenario
strucScenario <- 'pan'
#get the optQT ranges and max success rates
for(tarFst in TARFST){
  for(samSecT in SAMSEC){
    
    apRow <- which(pRes$pKey == paste(strucScenario, 'Allo', tarFst, samSecT, sep = '_'))
    spRow <- which(pRes$pKey == paste(strucScenario, 'Sec', tarFst, samSecT, sep = '_'))
    
    #across both invScenarios, find the qT-range where overall mean success rate is maximized
    comboSR <- (pRes$sr_ecdf[[apRow]](seq(0, 0.5, by = 0.01)) + 1 - pRes$sr_ecdf[[spRow]](seq(0, 0.5, by = 0.01))) / 2
    
    #store the optQT range and max SR
    optQT_ranges[[paste0(tarFst, '_', samSecT)]] <- seq(0, 0.5, by = 0.01)[range(which(abs(max(comboSR) - comboSR) < 0.00001))]
    max_sr[[paste(tarFst, samSecT, sep = '_')]] <- max(comboSR)
  }
}

#make plot
{
png(filename = paste0('13a.plots/3.optQT_ranges.png'), 
    width = 1200, height = 800)
par(mar = c(0.1, 0.1, 1.1, 0.1), xpd = NA, cex = 3)

plot(NA, xlim = c(-0.2, 0.5), ylim = c(-1.5, 8), bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
axis(1, at = seq(0, 0.5, by = 0.1), pos = 0)
text(0.25, -1.5, 'Q-threshold')
text(-0.17, 8.5, 'Fst')
text(-0.05, 8.5, 'Sampling\ntime (gens)')

idx <- 7.5
for(tarFst in TARFST){
  for(samSecT in SAMSEC){
    #add tarFst label
    if(samSecT == SAMSEC[1]){
      text(-0.17, idx - 0.5, tarFst)
      lines(x = rep(-0.11, 2), y = idx + c(0.2, -1.2), lwd = 4)
    }
    #add samSecT label
    text(-0.05, idx, samSecT)
    
    optQT_r <- optQT_ranges[[paste(tarFst, samSecT, sep = '_')]]
    arrows(x0 = optQT_r[1], 
           y0 = idx,
           x1 = ifelse(optQT_r[1] != optQT_r[2], optQT_r[2], optQT_r[1] + 0.0005),
           lwd = ifelse(max_sr[[paste0(tarFst, '_', samSecT)]] > 0.99, 6, 3),
           code = ifelse(optQT_r[1] != optQT_r[2], 3, 1), 
           length = 0.4, angle = 90)
    
    idx <- idx - 1
  }
}

#add shaded overlap area; see paper for details
rect(0.1, 0.2, 0.14, 5.8, col = adjustcolor('gray', alpha.f = 0.6), border = NA)

dev.off()
}
}

##Fig S4: allo/sec success rate over Q-threshold range, sub only
#using pan-selected optQTs
{
#plotting stuff
{
invLty <- c(Allo = 2, Sec = 1)
prettyNames <- c(Allo = 'No-invasion', Sec = 'Invasion')
#layout matrix
qgtLm <- matrix(0, nrow = 2 + length(TARFST), ncol = 2 + length(SAMSEC))
qgtLm[1, 2 + seq_along(SAMSEC)] <- max(qgtLm) + seq_along(SAMSEC) #samSecT column labels
qgtLm[1 + seq_along(TARFST), 1] <- max(qgtLm) + seq_along(TARFST) #tarFst row labels
qgtLm[1 + seq_along(TARFST), 2] <- max(qgtLm) + 1 #y-axis label
qgtLm[nrow(qgtLm), 2 + seq_along(SAMSEC)] <- max(qgtLm) + 1 #x-axis label
qgtLm[nrow(qgtLm), c(1,2)] <- max(qgtLm) + 1 #invScenario line type legend
qgtLm[1 + seq_along(TARFST), 2 + seq_along(SAMSEC)] <- max(qgtLm) + #data
  matrix(1:(length(SAMSEC)*length(TARFST)), nrow = length(TARFST), ncol = length(SAMSEC), byrow = F)
#panel dimensions
qgtW <- c(1.5, 2, rep(6, length(SAMSEC)))
qgtH <- c(1, rep(3.5, length(TARFST)), 1.5)
}

sel_optQT <- c('10' = 0.35, '30' = 0.3, 'rest' = 0.1)

#make plot
{
png(filename = paste0('13a.plots/s4.SR_sub.png'), 
    width = sum(qgtW)*100, height = sum(qgtH)*100)
layout(qgtLm, widths = qgtW, heights = qgtH)
par(mar = c(0.2, 0, 0.2, 0), cex = 2.5, xpd = NA)

#labels etc.
{
for(samSecT in SAMSEC){ #samSecT column labels
  plot.new()
  text(0.5, 0.6, cex = 1.2, paste(samSecT, 'gens.')); lines(c(0.05, 0.95), c(0.2, 0.2), lwd = 3)
}
for(i in seq_along(TARFST)){ #tarFst row labels
  plot.new(); text(0.7, 0.5, cex = 1.2, srt = 90, TARFST[i])
  lines(c(0.95, 0.95), c(0.05, 0.95), lwd = 3)
  if(i == ceiling(length(TARFST) / 2)){
    text(0.2, ifelse(length(TARFST)%%2==1, 0.5, 0.0), 
         'Target Fst level', srt = 90, font = 2, cex = 1.2)
  }
}
plot.new() #y-axis label
text(0.35, 0.5, cex = 1.2, srt = 90, 'Proportion of replicates successful')
plot.new() #x-axis label
text(0.5, 0.4, cex = 1.2, 'Q-threshold')
plot.new() #invScenario line type legend
legend('center', bty = 'n', legend = prettyNames[names(invLty)], 
       lty = invLty, lwd = 3, col = 'black'
)
}

par(mar = c(0.8, 0.6, 0.6, 0.2))
#plot data
strucScenario <- 'sub'
for(samSecT in SAMSEC){
  for(tarFst in TARFST){
    plot(NA, xlim = c(0, 0.6), ylim = c(0, 1.2), xaxt = 'n', yaxt = 'n', bty = 'n', xlab = '', ylab = '')
    axis(2, at = c(0, 0.5, 1), las = 1, labels = (samSecT == SAMSEC[1]), lwd = 2)
    axis(1, at = seq(0, 0.5, by = 0.1), labels = (tarFst == max(TARFST)), lwd = 2)
    
    for(invScenario in c('Allo', 'Sec')){
      pRow <- which(pRes$pKey == paste(strucScenario, invScenario, tarFst, samSecT, sep = '_'))
      
      #plot curve
      ec <- pRes$sr_ecdf[[pRow]]
      if(invScenario == 'Allo'){
        curve(ec(x), from = 0, to = 0.5, lty = invLty[invScenario], lwd = 3.5, add = T)
        sr_Allo <- ec(seq(0.01, 0.49, by = 0.01))
      } else {
        curve(1 - ec(x), from = 0, to = 0.5, lty = invLty[invScenario], lwd = 3.5, add = T)
        sr_Sec <- 1 - ec(seq(0.01, 0.49, by = 0.01))
      }
      
      #display number of simReps that got K=2
      text(0.55, ifelse(invScenario == 'Allo', 1, 0.05), length(pRes$simRepsK2[[pRow]]))
    }
    
    #get selected optimal Q-threshold for this scenario
    optQT <- ifelse(tarFst == 0.01, sel_optQT[paste(samSecT)], sel_optQT['rest'])
    
    #get average SR across both invScenarios at optQT
    comboSR <- (sr_Allo + sr_Sec) / 2
    avSR <- comboSR[as.integer(optQT*100)]
    
    #display max mean success rate, and the optimal qT-range
    text(optQT, 1.2, round(avSR, digits = 2))
    points(optQT, avSR, pch = 3, cex = 2, lwd = 4)
    text(optQT - 0.01, 1.1, optQT, adj = 1)
  }
}
dev.off()
}
}

##Figs S5 (pan)/S6 (sub): correct individual-ID rates over Q-threshold range
{
#plotting stuff
{
invLty <- c(Allo = 2, Sec = 1)
prettyNames <- c(Allo = 'No-invasion', Sec = 'Invasion')
#layout matrix
qgtLm <- matrix(0, nrow = 2 + length(TARFST), ncol = 2 + length(SAMSEC))
qgtLm[1, 2 + seq_along(SAMSEC)] <- max(qgtLm) + seq_along(SAMSEC) #samSecT column labels
qgtLm[1 + seq_along(TARFST), 1] <- max(qgtLm) + seq_along(TARFST) #tarFst row labels
qgtLm[1 + seq_along(TARFST), 2] <- max(qgtLm) + 1 #y-axis label
qgtLm[nrow(qgtLm), 2 + seq_along(SAMSEC)] <- max(qgtLm) + 1 #x-axis label
qgtLm[nrow(qgtLm), c(1,2)] <- max(qgtLm) + 1 #invScenario line type legend
qgtLm[1 + seq_along(TARFST), 2 + seq_along(SAMSEC)] <- max(qgtLm) + #data
  matrix(1:(length(SAMSEC)*length(TARFST)), nrow = length(TARFST), ncol = length(SAMSEC), byrow = F)
#panel dimensions
qgtW <- c(1.5, 2, rep(6, length(SAMSEC)))
qgtH <- c(1, rep(3.5, length(TARFST)), 2)
}

#pan-selected optQT:
sel_optQT <- c('10' = 0.35, '30' = 0.3, 'rest' = 0.1)

#make plots (one for each strucScenario)
for(strucScenario in c('pan', 'sub')){
  png(filename = paste0('13a.plots/s', ifelse(strucScenario == 'pan', 5, 6), '.cID_', strucScenario, '.png'), 
      width = sum(qgtW)*100, height = sum(qgtH)*100)
  layout(qgtLm, widths = qgtW, heights = qgtH)
  par(mar = c(0.2, 0, 0.2, 0), cex = 2.5, xpd = NA)
  
  #labels etc.
  {
  for(samSecT in SAMSEC){ #samSecT column labels
    plot.new()
    text(0.5, 0.6, cex = 1.2, paste(samSecT, 'gens.')); lines(c(0.05, 0.95), c(0.2, 0.2), lwd = 3)
  }
  for(i in seq_along(TARFST)){ #tarFst row labels
    plot.new(); text(0.7, 0.5, cex = 1.2, srt = 90, TARFST[i])
    lines(c(0.95, 0.95), c(0.05, 0.95), lwd = 3)
    if(i == ceiling(length(TARFST) / 2)){
      text(0.2, ifelse(length(TARFST)%%2==1, 0.5, 0.0), 
           'Target Fst level', srt = 90, font = 2, cex = 1.2)
    }
  }
  plot.new() #y-axis label
  text(0.35, 0.5, cex = 1.2, srt = 90, 'Mean proportion of individuals correctly identified')
  plot.new() #x-axis label
  text(0.5, 0.5, cex = 1.2, 'Q-threshold')
  plot.new() #legend: invScenario line type, strucScenario colour, optimal qT linetype
  legend('bottomleft', bty = 'n', ncol = 2, text.width = 0.65,
         legend = c(prettyNames, '', 'Optimal Q-threshold'), 
         lty = c(invLty, NA, 3), lwd = c(3, 3, NA, 3)
  )
  }
  
  par(mar = c(0.8, 0.6, 0.6, 0.6))
  #plot data
  for(samSecT in SAMSEC){
    for(tarFst in TARFST){
      plot(NA, xlim = c(0, 0.5), ylim = c(0, 1.2), xaxt = 'n', yaxt = 'n', bty = 'n', xlab = '', ylab = '')
      axis(2, at = c(0, 0.5, 1), las = 1, labels = (samSecT == SAMSEC[1]), lwd = 2)
      axis(1, at = seq(0, 0.5, by = 0.1), labels = (tarFst == max(TARFST)), lwd = 2)
      
      for(invScenario in c('Allo', 'Sec')){
        
        pKey <- paste(strucScenario, invScenario, tarFst, samSecT, sep = '_')
        pRow <- which(pRes$pKey == pKey)
        
        lines(x = seq(0.01, 0.49, by = 0.01), y = pRes$cID[[pRow]],
              lty = invLty[invScenario], lwd = 3)
      }
      
      #display the selected optimal qT
      optQT <- ifelse(tarFst == 0.01, sel_optQT[as.character(samSecT)], sel_optQT['rest'])
      lines(x = rep(optQT, 2), y = c(0, 1), lty = 3, lwd = 2)
      text(optQT, 1.1, optQT)
    }
  }
  dev.off()
}
}

logMsg('good job 3.main results section')
rm(list = ls()[!ls() %in% notDelete])
}

######### 4.plot heatmap using optimal Q-thresholds ############################

#plot fig S7: individual ID heatmap
{
logMsg('start 4.plot heatmap using optimal Q-thresholds')
#plot only displays data from simReps that got K=2 in deltaK, though this first
# part (making the counts table) gets data from all simReps.
# simReps that didn't get K=2 are filtered out when making the plot.

#pan-selected optQT, selected in section 3:
sel_optQT <- c('10' = 0.35, '30' = 0.3, 'rest' = 0.1)

#get data
{
#make counts table with one row per uKey, containing counts of individuals of
# each possible combination of true ID and inferred label

tRows <- which(tRes$K == 2 & tRes$isMaxLK == 1)
counts <- tRes[tRows, c('simRep', 'strucScenario', 'invScenario', 'tarFst', 'samSecT', 'pKey', 'uKey')]

#add columns for number of individuals with each true-inferred label combination
allIDs <- c('rec', 'hyb', 'inv')
allIDcombos <- with(expand.grid(allIDs, allIDs), paste(Var1, Var2, sep = '_'))
counts[, allIDcombos] <- NA

#get counts of individuals of each true-inferred label combo
for(row in 1:nrow(counts)){
  uKey <- counts$uKey[row]
  iData <- indData[[uKey]]
  
  #determine which cluster probably should correspond to inv
  invCl <- which.max(c(cor(iData$qAnc1, iData$tAncInv), cor(iData$qAnc2, iData$tAncInv)))
  
  #get the optimal Q-threshold of this tarFst/samSecT scenario
  optQT <- ifelse(counts$tarFst[row] == 0.01, sel_optQT[as.character(counts$samSecT[row])], sel_optQT['rest'])
  
  #get each individual's STRUCTURE-inferred ID with this optQT (rec, inv, hyb)
  qID <- ifelse(iData[[paste0('qAnc', invCl)]] > (1 - optQT), 'inv', 
         ifelse(iData[[paste0('qAnc', invCl)]] < optQT, 'rec', 'hyb'))
  
  #make combined true_inferred labels (trueID_qID)
  tqID <- paste(iData$trueID, qID, sep = '_')
  
  for(comboID in allIDcombos){
    counts[row, comboID] <- sum(tqID == comboID)
  }
}
}

#plotting stuff
{
#laymat
hpLm <- matrix(0, nrow = 3 + length(TARFST), ncol = 3 + (4*length(SAMSEC)))
hpLm[1, 3:(ncol(hpLm) - 1)] <- max(hpLm) + rep(c(1,2), each = 2*length(SAMSEC)) #strucScenario column labels
hpLm[2, 3:(ncol(hpLm) - 1)] <- max(hpLm) + rep(c(1:(2*length(SAMSEC))), each = 2) #samSecT column labels
hpLm[2 + seq_along(TARFST), 1] <- max(hpLm) + seq_along(TARFST) #tarFst row labels
hpLm[3:(nrow(hpLm) - 1), 2] <- max(hpLm) + 1 #y-axis label
hpLm[nrow(hpLm), 3:(ncol(hpLm) - 1)] <- max(hpLm) + 1 #x-axis label
hpLm[nrow(hpLm), 1:2] <- max(hpLm) + 1 #invScenario lty legend
hpLm[1:2, ncol(hpLm)] <- max(hpLm) + 1 #optimal Q-thresholds column label
hpLm[3:(nrow(hpLm) - 1), ncol(hpLm)] <- max(hpLm) + 1 #display optimal Q-thresholds used
hpLm[3:(nrow(hpLm) - 1), 3:(ncol(hpLm) - 1)] <- max(hpLm) + 
  matrix(1:(length(SAMSEC)*length(TARFST)*4), nrow = length(TARFST), ncol = length(SAMSEC)*4, byrow = F)

#dims
hpW <- c(1.25, 1.25, rep(c(3, 4), 2*length(SAMSEC)), 4)
hpH <- c(1, 1, rep(4, length(TARFST)), 1.5)

#colours and names
strucCols <- c(pan = '#5B8FA8', sub = '#B1746F')
invLty <- c(Allo = 2, Sec = 1)
prettyNames <- c(Allo = 'No-invasion', Sec = 'Invasion', pan = 'Panmictic', sub = 'Substructured',
                 rec = 'Recipient', hyb = 'Hybrid', inv = 'Invader')
allIDs <- c('rec', 'hyb', 'inv')
#replaced recipient with native in text, so using nat for labels, but data still uses rec
labIDs <- c('nat', 'hyb', 'inv')
}

#make plot
{
png(filename = paste0('13a.plots/s7.heatmap.png'), 
    width = sum(hpW)*200, height = sum(hpH)*200)
layout(hpLm, widths = hpW, heights = hpH)
par(cex = 7, mar = c(0.1, 0.1, 0.1, 0.1), mgp = c(3, 0, 0), xpd = NA)

#labels etc.
{
for(strucScenario in c('pan', 'sub')){
  plot.new()
  text(0.5, 0.7, cex = 1.2, prettyNames[strucScenario]); lines(c(0.05, 0.95), c(0.3, 0.3), lwd = 8)
}
for(samSecT in rep(SAMSEC, 2)){ #samSecT column labels
  plot.new()
  text(0.5, 0.7, cex = 1.2, paste(samSecT, 'gens.')); lines(c(0.05, 0.95), c(0.3, 0.3), lwd = 8)
}
for(i in seq_along(TARFST)){ #tarFst row labels
  plot.new(); text(0.75, 0.5, srt = 90, cex = 1.2, TARFST[i])
  if(i == ceiling(length(TARFST) / 2)){
    text(0.25, ifelse(length(TARFST)%%2==1, 0.5, 0.0), 
         'Target Fst level', srt = 90, font = 2, cex = 1.2)
  }
}
plot.new() #y-axis label
text(0.5, 0.5, srt = 90, cex = 1.2, 'STRUCTURE-inferred identity category')
plot.new() #x-axis label
text(0.5, 0.3, cex = 1.2, 'True identity category')
plot.new() #invScenario lty legend
legend('bottomleft', bty = 'n', legend = prettyNames[names(invLty)], 
       lty = invLty, lwd = 8, col = 'black', ncol = 2, cex = 1.2
)
plot.new() #column label for Q-thresholds used
text(0.5, 0.5, 'Q-thresholds used:')
par(mar = c(0, 0, 0, 0))
plot.new() #display Q-thresholds used
#Fst = 0.01 optQTs
text(0.5, 0.875, '10 gens: 0.35\n30 gens: 0.30')
lines(c(0.01, 0.01), c(0.775, 0.975), lwd = 8)
#Fst > 0.01 optQTs
text(0.5, 0.375, '0.10')
lines(c(0.01, 0.01), c(0.075, 0.675), lwd = 8)
}

#heatmap panels
for(strucScenario in c('pan', 'sub')){
  colsPropID <- colorRampPalette(c('white', strucCols[strucScenario]))(100)
  for(samSecT in SAMSEC){
    for(invScenario in c('Allo', 'Sec')){
      par(mar = c(0.2, ifelse(invScenario == 'Allo', 0.6, 0.2), 0.2, ifelse(invScenario == 'Allo', 0.2, 0.6)))
      for(tarFst in TARFST){
        if(invScenario == 'Allo'){
          plot(NA, type = 'n', xlim = c(0.5, 2.5), ylim = c(0.5, 3.5), 
               xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')
          if(tarFst==max(TARFST)){axis(1, at = 1:2, labIDs[c(1, 3)], tick = F)}
          if(samSecT==min(SAMSEC) & strucScenario=='pan'){
            axis(2, at = 1:3, labIDs, tick = F)
          }
          #true_inf ID categories
          for(tID in c(1, 3)){
            tID_pos <- ifelse(tID == 3, 2, tID)
            for(qID in 1:3){
              props <- c()
              for(simRep in pRes$simRepsK2[[which(pRes$pKey == paste(strucScenario, 'Allo', tarFst, samSecT, sep = '_'))]]){
                cRow <- which(counts$uKey == paste(simRep, strucScenario, 'Allo', tarFst, samSecT, sep = '_'))
                props <- c(props, counts[cRow, paste(allIDs[c(tID,qID)], collapse = '_')]/totIndSampSize)
              }
              col_idx <- ceiling(mean(props)*99)+1
              rect(xleft = tID_pos - 0.5, ybottom = qID - 0.5,
                   xright = tID_pos + 0.5, ytop = qID + 0.5,
                   col = colsPropID[col_idx], border = 'grey40', lwd = 8, lty = invLty[invScenario])
              if(mean(props)>0){text(tID_pos, qID, round(mean(props), 2))}
            }
          }
          #highlight true squares
          rect(xleft = 0.5, ybottom = 0.5, xright = 1.5, ytop = 1.5, border = 'black', lwd = 10)
          rect(xleft = 1.5, ybottom = 2.5, xright = 2.5, ytop = 3.5, border = 'black', lwd = 10)
          
        } else if(invScenario == 'Sec'){
          plot(NA, type = 'n', xlim = c(0.5, 3.5), ylim = c(0.5, 3.5), 
               xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')
          if(tarFst==max(TARFST)){axis(1, at = 1:3, labIDs, tick = F)}
          #true_inf ID categories
          for(tID in 1:3){
            for(qID in 1:3){
              props <- c()
              for(simRep in pRes$simRepsK2[[which(pRes$pKey == paste(strucScenario, 'Sec', tarFst, samSecT, sep = '_'))]]){
                cRow <- which(counts$uKey == paste(simRep, strucScenario, 'Sec', tarFst, samSecT, sep = '_'))
                props <- c(props, counts[cRow, paste(allIDs[c(tID,qID)], collapse = '_')]/totIndSampSize)
              }
              col_idx <- ceiling(mean(props)*99)+1
              rect(xleft = tID - 0.5, ybottom = qID - 0.5,
                   xright = tID + 0.5, ytop = qID + 0.5,
                   col = colsPropID[col_idx], border = 'grey40', lwd = 8, lty = invLty[invScenario])
              if(mean(props)>0){text(tID, qID, round(mean(props), 2))}
            }
          }
          #highlight true squares
          for(i in 1:3){
            rect(xleft = i-0.5, ybottom = i-0.5, xright = i+0.5, ytop = i+0.5, border = 'black', lwd = 10)
          }
        }
      }
    }
  }
}

dev.off()
}

rm(list = ls()[!ls() %in% notDelete])
logMsg('good job 4.plot heatmap using optimal Q-thresholds')
}

######### 5.plot Q-real vs Q-STRUCTURE #########################################

#plots Figs S8-S15: STRUCTURE-inferred Q-values vs true invader-ancestry
#proportions. Data plotted for all individuals sampled from simReps that got K=2
{
logMsg('start 5.plot QvQ')
dir.create('13b.QvQplots', showWarnings = FALSE)

#plotting Q-STRUCTURE vs Q-real for the invader pop for each pKey
#heatmap of all individuals from all simReps (that got K=2) on each plot
for(tarFst in TARFST){
  for(samSecT in SAMSEC){
    
    jpeg(filename = paste0('13b.QvQplots/QvQ_K2_', sprintf('%.2f', tarFst), '_', samSecT, '.jpeg'), 
         quality = 100, width = 2000, height = 2000)
    layout(
      matrix(c(0, 1, 2, 0, 5, 8, 9, 3, 5, 10, 11, 4, 7, 6, 6, 0), nrow = 4, byrow = T), 
      widths = c(1, 6, 6, 1), heights = c(1, 6, 6, 1))
    par(cex = 4, xpd = NA, mar = c(0.1, 0.1, 0.1, 0.1))
    
    #labels
    {
    #column labels: population structure scenario
    plot.new()
    text(0.5, 0.5, cex = 1.2, 'Panmictic')
    plot.new()
    text(0.5, 0.5, cex = 1.2, 'Substructured')
    #row labels: invasion scenario
    plot.new()
    text(0.5, 0.5, cex = 1.2, srt = 270, 'No-invasion')
    plot.new()
    text(0.5, 0.5, cex = 1.2, srt = 270, 'Invasion')
    
    plot.new() #y-axis label
    text(0.5, 0.5, cex = 1.2, srt = 90, 'Inferred invader ancestry')
    plot.new() #x-axis label
    text(0.5, 0.5, cex = 1.2, 'True invader ancestry')
    
    plot.new() #show optimal Q-threshold
    #decided against this
    #text(0.2, 0.5, adj = 0, paste0('Optimal Q-threshold: ', ifelse(tarFst == 0.01, ifelse(samSecT == 10, '0.35', '0.30'), 0.10)))
    }
    
    par(cex = 3, mar = c(2.1, 2.1, 0.6, 0.6))
    for(invScenario in c('Allo', 'Sec')){
      for(strucScenario in c('pan', 'sub')){
        
        #initialize vectors in which to store data from all simReps
        #(faster to do this, and combine them into a dataframe after)
        tInv_v <- numeric(); qInv_v <- numeric(); pop_v <- integer()
        
        #which simReps got K=2
        k2SimReps <- pRes$simRepsK2[[which(pRes$pKey == paste(strucScenario, invScenario, tarFst, samSecT, sep = '_'))]]
        #get data from all simReps that got K=2
        for(simRep in k2SimReps){
          udat <- indData[[paste(simRep, strucScenario, invScenario, tarFst, samSecT, sep = '_')]]
          
          #which STRUCTURE cluster best corresponds to the real invader pop
          invClus <- which.max(c(cor(udat$tAncInv, udat$qAnc1), cor(udat$tAncInv, udat$qAnc2)))
          
          #add data to vectors
          tInv_v <- c(tInv_v, udat$tAncInv)
          qInv_v <- c(qInv_v, udat[[paste0('qAnc', invClus)]])
          pop_v <- c(pop_v, udat$pop)
        }
        
        #combine data from all simReps for this pKey
        pdat <- data.frame(
          tInv = tInv_v, qInv = qInv_v, pop = pop_v
        )
        
        #separate bins for tInv=0 or 1
        nbins <- 50
        int <- 1/nbins
        pdat$tInvMod <- pdat$tInv
        pdat$tInvMod[pdat$tInvMod == 0] <- -int/2
        pdat$tInvMod[pdat$tInvMod == 1] <- 1 + (int/2)
        xbreaks <- seq(-int, 1 + (int), length.out = nbins + 3)
        ybreaks <- seq(0, 1, length.out = nbins + 1)
        cts <- table(
          cut(pdat$tInvMod, xbreaks, include.lowest = T),
          cut(pdat$qInv, ybreaks, include.lowest = T)
        )
        log_cts <- log10(cts)
        log_cts[!is.finite(log_cts)] <- NA
        image(
          x = xbreaks, y = ybreaks, z = log_cts, xlab = '', ylab = '',
          col = colorRampPalette(c('#F5F8EA', '#74D3C0', '#2F327D'))(100)
        )
        
        #display the pan-selected optimal qT
        optQT <- ifelse(tarFst == 0.01, ifelse(samSecT == 10, 0.35, 0.3), 0.1)
        lines(x = range(xbreaks), y = rep(optQT, 2), lty = 3, lwd = 2)
        lines(x = range(xbreaks), y = 1 - rep(optQT, 2), lty = 3, lwd = 2)
        
        #display the K=2 sample size in bottom-right corner
        text(0.9, 0.1, cex = 1.5, length(k2SimReps))
      }
    }
    dev.off()
  }
}

rm(list = ls()[!ls() %in% notDelete])
logMsg('good job 5.plot QvQ')
}
