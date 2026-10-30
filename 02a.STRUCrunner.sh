#!/bin/bash
logS() { echo "$(date '+%F %T') $*" >> 20a.mainLog.txt; }
logS "2.$arrayID started: jobID=$SLURM_JOB_ID"
start_time=$(date +%s)

#get parameters from 20b.jobMap.txt	
IFS=$'\t' read -r arrayID status jobID parentJobID simRep invScenario \
    strucScenario tarFst samSecT K strucRep seed tAncFP snpFP runTime \
    < <(awk -v id="$arrayID" -F $'\t' 'NR>1 && $1==id {print; exit}' 20b.jobMap.txt)

outFile="12.STRUC_outputs/arrayID${arrayID}"
numPolySites=$(awk 'NR==1 {print NF - 1; exit}' "$snpFP")
totIndSampSize=$(( $(wc -l < "$snpFP") / 2 ))
echo "$(date '+%F %T') Running STRUCTURE: jobID=$SLURM_JOB_ID | arrayID=$arrayID |\
$strucScenario$invScenario | simRep=$simRep | K=$K | strucRep=$strucRep | seed=$seed |\
numPolySites=$numPolySites | totIndSampSize=$totIndSampSize"

module load StdEnv/2020 gcc/9.3.0 structure/2.3.4 && \
structure -m "02b.mainParams.txt" \
	-e "02c.extraParams.txt" \
	-K "$K" \
	-N "$totIndSampSize" \
	-L "$numPolySites" \
	-i "$snpFP" \
	-o "$outFile" \
	-D "$seed"
if [[ $? -eq 0 ]]; then
	echo "$(date '+%F %T') STRUCTURE succeeded"
	elapsed_time=$(( $(date +%s) - start_time ))
	#set jobMap status to comp, and log elapsed runtime (in seconds)
	( flock -x 200
		awk -F'\t' -v OFS='\t' -v ai="$arrayID" -v et="$elapsed_time" '
			$1==ai { $2="comp"; $15=et } {print}
			' 20b.jobMap.txt > 20b.jobMap.tmp && mv 20b.jobMap.tmp 20b.jobMap.txt
	) 200>20b.jobMap.lock
	#rename STRUCTURE full output file to have .txt at the end
	mv "${outFile}_f" "${outFile}_f.txt"
else
	echo "$(date '+%F %T') STRUCTURE failed"
	logS "2.$arrayID failed: jobID=$SLURM_JOB_ID"
	#set jobMap status to fail
	( flock -x 200
		awk -F'\t' -v OFS='\t' -v ai="$arrayID" '$1==ai { $2="fail" } {print}' \
			20b.jobMap.txt > 20b.jobMap.tmp && mv 20b.jobMap.tmp 20b.jobMap.txt
	) 200>20b.jobMap.lock
	exit 1
fi