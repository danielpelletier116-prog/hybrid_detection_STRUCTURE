import sys, os, atexit, warnings, tskit, msprime, pyslim
import numpy as np
import pandas as pd
from datetime import datetime
from itertools import combinations
import random

startTime = datetime.now()

def make_logFun(writeOption):
    def logMsg(*messages, time_since=startTime):
        elapsed_sec = (datetime.now() - time_since).total_seconds()
        if elapsed_sec < 60:
            elapsed = f'{round(elapsed_sec)}s'
        else:
            elapsed = f'{round(elapsed_sec/60, 1)}m'
        writeOption(''.join([f"{elapsed}: {msg}\n" for msg in messages]))
    return logMsg
    
simRep = sys.argv[1] #this is a string, padded with 0s for 2 characters

logf = open(f'10.SLURM_outputs/1_simRep{simRep}.txt', 'a', buffering=1)
atexit.register(logf.close)
logMsg = make_logFun(logf.write)
logMsg(f'START {startTime.strftime('%Y-%m-%d %H:%M:%S')}', f'READING run parameters for simRep{simRep}')

runParamDict = {}
with open('00b.runParameters.txt') as f:
    exec(f.read(), {}, runParamDict)
N = int(runParamDict['N'])
r = float(runParamDict['R'])
u = float(runParamDict['U'])
L = int(runParamDict['L'])
tBur = int(runParamDict['TBUR'])
tAllo = int(runParamDict['TALLO'])
samAllo = int(runParamDict['SAMALLO'])
samSec = [int(item) for item in runParamDict['SAMSEC'].split(',')]
tarFst = [float(item) for item in runParamDict['TARFST'].split(',')]
recIndSampSize = int(runParamDict['RECINDSAMPSIZE'])
invIndSampSizeRatio = float(runParamDict['INVINDSAMPSIZERATIO'])
locSampSize = int(runParamDict['LOCSAMPSIZE'])

random_seed = np.random.randint(1, 2**31)
logMsg(f'GENERATED random seed for msprime mutations: {random_seed}')

ts = tskit.load(f'11a.treeSeqs/treeSeq_rep{simRep}.trees')
logMsg(f'TREESEQ contains {ts.num_trees} trees')
if any(t.num_roots > 1 for t in ts.trees()):
    logMsg('RECAPITATING')
    warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)
    rts = pyslim.recapitate(ts, recombination_rate=r, ancestral_Ne=N)
    logMsg(f'RECAP added {rts.num_trees - ts.num_trees} trees',
        f'RECAP removed {max(t.num_roots for t in ts.trees())-max(t.num_roots for t in rts.trees())} max roots',
        f'RECAP added {round(rts.max_root_time - ts.max_root_time)} root time',
        f'RECAP tree max age: {round(rts.max_root_time)}')
else:
    logMsg('SKIPPING recapitation')
    rts = ts
    logMsg(f'NO-RECAP tree max age: {round(rts.max_root_time)}')
mrts = msprime.sim_mutations(rts, rate=u, model=msprime.BinaryMutationModel(state_independent=False), random_seed = random_seed)
logMsg(f'MUT added {mrts.num_mutations} neutral mutations',
    f'MUT at {mrts.num_sites} sites'
)
del [ts, rts]

mrts.dump(f'11a.treeSeqs/treeSeq_rep{simRep}.trees')
logMsg('SAVED mrts, overwrote ts')

#get allopatry sampling/2nd contact initialization times
slimSamAlloT = np.arange(tBur, tBur+tAllo+samAllo, samAllo)
pySamAlloT = pyslim.slim_time(mrts, time = slimSamAlloT)

#set allopatry popIDs
invader = [1]
alloRecPan = [2]
alloRecSub = [3,4,5]
popD = {(['inv', 'pan', 'sub'] + [f'sub{i}' for i in alloRecSub])[i] : 
    (invader+alloRecPan+[alloRecSub]+alloRecSub)[i] for i in range(6)}
allAllo = [1,2,3,4,5]

#dictionary for Fst data:
#hFst=Hudson's Fst (equivalent to WC84 in this scenario), between the invader pop and either the pan or sub recipient pop
#demeHFst=mean pairwise Fst between subpops
fData = {
    'hFst_pan' : [],
    'hFst_sub' : [],
    'demeHFst_sub' : []
}

for counter, (slimT, pyT) in enumerate(zip(slimSamAlloT, pySamAlloT), start=1):
    logMsg(f'GET Fst data for slimT {slimT}, {counter}/{len(slimSamAlloT)}')
    
    #get simplified tree for all sample nodes at pyT, and node IDs for each population/grouping
    fts = mrts.simplify(
        samples = np.where((mrts.nodes_time==pyT) & (mrts.nodes_flags==1) & np.isin(mrts.nodes_population, allAllo))[0], 
        filter_populations = False)
    nodes = {l : np.where((fts.nodes_time==pyT) & (fts.nodes_flags==1) & np.isin(fts.nodes_population, p))[0]
         for l, p in popD.items()}
    
    #get Fst between invader pop and either the pan or sub recipient pop
    for l in ['pan', 'sub']:
        n = nodes[l]
        tHw = (fts.diversity(nodes['inv'], span_normalise=False) + fts.diversity(n, span_normalise=False)) / 2
        tHb = fts.divergence([nodes['inv'], n], span_normalise=False)
        fData[f'hFst_{l}'].append(1 - (tHw/tHb))
    
    #get mean pairwise Fst between subpops
    mean_tHw = sum([
        (fts.diversity(nodes[f'sub{i}'], span_normalise=False) + 
        fts.diversity(nodes[f'sub{j}'], span_normalise=False)) / 2
        for i,j in combinations(alloRecSub, 2)]) / 3
    mean_tHb = sum([fts.divergence([nodes[f'sub{i}'], nodes[f'sub{j}']], span_normalise=False)
        for i,j in combinations(alloRecSub, 2)]) / 3
    fData['demeHFst_sub'].append(1 - (mean_tHw/mean_tHb))

logMsg('GOT Fst data')
with open(f'21b.Fst.txt', 'a', buffering = 1) as fstf:
    for metric, data in fData.items():
        strData = " ".join(f'{v:.10f}' for v in data)
        fstf.write(f'{simRep}_{metric} {strData}\n')
logMsg('WROTE Fst data')
for var in ['slimT', 'pyT', 'counter', 'fts', 'nodes', 'l', 'n', 'tHw', 'tHb', 'mean_tHw', 'mean_tHb']:
    locals().pop(var, None)

#read in 21a.popIDs.txt table which contains: tick (slim time when 2nd contact starts), popIds of pan/subSec
popIDs = pd.read_csv('21a.popIDs.txt', sep='\t', dtype={0: str})
popIDs['sub'] = popIDs[['sub3', 'sub4', 'sub5']].values.tolist()

#output true ancestry files and SNP files
for strucScenario, alloRec in zip(['pan', 'sub'], [alloRecPan, alloRecSub]):
    Fst_series = np.array(fData[f'hFst_{strucScenario}'], dtype = float)
    for tarF in tarFst: #for each target Fst, find the allopatry sample time with Fst closest to target
        selT, pySelT = slimSamAlloT[np.abs(Fst_series - tarF).argmin()], pySamAlloT[np.abs(Fst_series - tarF).argmin()]
        logMsg(f'START {strucScenario}, tarF={tarF}, selT={selT}, Fst={round(Fst_series[np.abs(Fst_series - tarF).argmin()],3)}')
        #get nodes of pre-2nd contact ancestral invader population
        invAncInds = pyslim.individuals_alive_at(mrts, time = pySelT, population = invader)
        invAncNodes = [int(node) for ind in invAncInds for node in mrts.individual(ind).nodes]
        #logMsg(f'GOT ancestral invader nodes for selT {selT}')
        #get the 2nd contact recipient popIDs
        secRec = np.array(popIDs.loc[popIDs.tick == str(selT), strucScenario].values[0])
        for samSecT in samSec:
            slimT, pyT = selT + samSecT, pySelT - samSecT
            invAlives = pyslim.individuals_alive_at(mrts, pyT, population = invader) #get IDs of invader pop inds at this time
            #split for allo vs sec file production
            for invScenario, recPop in zip(['Allo', 'Sec'], [alloRec, secRec]):
                recAlives = pyslim.individuals_alive_at(mrts, pyT, population = recPop) #get IDs of recipient pop inds at this time
                
                #sample individuals from invader and rec alives, and get their node IDs
                invSampInds = np.random.choice(invAlives, size=max(round(invIndSampSizeRatio*recIndSampSize), 1), replace=False)
                recSampInds = np.random.choice(recAlives, size=recIndSampSize, replace=False)
                sampInds = np.concatenate([invSampInds, recSampInds])
                sampNodes = [int(node) for ind in sampInds for node in mrts.individual(ind).nodes]
                totIndSampSize = len(sampInds)
                logMsg(f'SAMPLED {totIndSampSize} total inds')
                
                #make simplified tree with only the relevant individuals
                simpNodes = np.concatenate([invAncNodes, sampNodes])
                sts, nodeMap = mrts.simplify(samples = simpNodes, map_nodes = True, filter_populations = False)
                newInvAncNodes = nodeMap[invAncNodes]
                newSampNodes = nodeMap[sampNodes]
                newSampInds = np.unique([sts.node(n).individual for n in newSampNodes])
                logMsg('MADE simplified tree')
                
                #table for ancestry of invader pop in all sampled nodes
                la = sts.tables.link_ancestors(samples = newSampNodes, ancestors = newInvAncNodes)
                logMsg('GOT invader-ancestry table')
                #get true invader ancestry of sampInds and write to file
                #(even though it's obvious for any inv/allo inds)
                tAncFilePath = (f'11b.true_ancestry/tAnc_{strucScenario}{invScenario}_simRep{simRep}'
                    f'_tarF{tarF}_samT{str(samSecT).zfill(2)}.txt')
                with open(tAncFilePath, 'w') as tancf:
                    for ind in newSampInds:
                        label = f"p{sts.individual(ind).population}i{ind}"
                        true_inv_anc = (sum(sum(la.right[la.child==node]-la.left[la.child==node])
                            for node in sts.individual(ind).nodes) / (2*L))
                        true_nat_anc = 1.0 - true_inv_anc
                        tancf.write(f"{label} {true_inv_anc} {true_nat_anc}\n")
                logMsg(f'WROTE true ancestry to {tAncFilePath}')
                
                #sample loci
                numSampSites = round(locSampSize * sts.num_sites / sts.sequence_length)
                sampLoci = np.sort(np.random.choice(np.arange(sts.num_sites), numSampSites, replace=False))
                polySites = sampLoci[np.where(sts.segregating_sites(sample_sets=newSampNodes, windows="sites", span_normalise=False)[sampLoci])]
                logMsg(f'GOT sites: {sts.num_sites} total, {numSampSites} sampled, {len(polySites)} poly in indivSample')
                
                #write SNP file with polymorphic sampled loci
                snpFilePath = (f'11c.SNP_files/SNP_{strucScenario}{invScenario}_simRep{simRep}_tarF{tarF}_'
                    f'samT{str(samSecT).zfill(2)}.txt')
                #get the table of variants for the sampled nodes
                vars = list(sts.variants(samples = newSampNodes))
                logMsg('GOT vars')
                varsPos = {var.site.position: var for var in vars}
                logMsg('GOT varsPos')
                
                #apply minor allele count: only keep sites from polySites where minor allele count is 3 or more
                macPolyPos = [] #storing the actual genome position, rather than site ID
                for polySite in polySites:
                    varGeno = varsPos.get(sts.site(polySite).position).genotypes
                    mac = min(np.sum(varGeno), len(varGeno) - np.sum(varGeno))
                    if mac >= 3:
                        macPolyPos.append(sts.site(polySite).position)
                logMsg(f'GOT positions of {len(macPolyPos)} poly sites with MAC >=3, rem {len(polySites)-len(macPolyPos)} sing/doub')
                
                with open(snpFilePath, 'w') as snpf:
                    for nodeIndex in np.arange(len(newSampNodes)):
                        node = newSampNodes[nodeIndex]
                        label = f"p{sts.node(node).population}i{sts.node(node).individual}"
                        genos = []
                        for macPolyP in macPolyPos:
                            var = varsPos.get(macPolyP)
                            genos.append(str(var.genotypes[nodeIndex]))
                        geno = " ".join(genos)
                        snpf.write(f"{label} {geno}\n")
                logMsg(f'WROTE SNP file to {snpFilePath}')
                del [sts, vars, varsPos]

logMsg('FINISHED writing all tAnc and SNP files')