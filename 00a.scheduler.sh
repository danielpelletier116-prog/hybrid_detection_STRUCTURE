#!/bin/bash
#SBATCH --mem-per-cpu=1000
#SBATCH --time=01:00:00
#SBATCH --job-name=0.1
#SBATCH --output=20a.mainLog.txt
#SBATCH --open-mode=append

logS() { echo "$(date '+%F %T') $*" >> 20a.mainLog.txt; }
source 00b.runParameters.txt

NJOBSSUB=0 #count of jobs scheduled
SCIDS="" #IDs of all jobs scheduled

if [[ -z $NSCHED ]]; then #if this is the first scheduler
	find . -type f -print0 | xargs -0 dos2unix
	mkdir -p {10.SLURM_outputs,11a.treeSeqs,11b.true_ancestry,11c.SNP_files,12.STRUC_outputs}
	NSCHED=1 #start count of schedulers
	
	#create jobMap, update mainParams if needed
	module load StdEnv/2023 r/4.4.0
	Rscript 00c.mapper.R
	logS "Made 20b.jobMap.txt"
	logS "$NUMSIMREPS simulations plus $(tail -n +2 20b.jobMap.txt | wc -l) STRUCTURE jobs total"
	
	#schedule simulation jobs
	for ((i=1; i<=NUMSIMREPS; i++)); do
		simRep=$(printf "%02d" "$i")
		
		SIMSCRIPT=$(mktemp)
cat > "$SIMSCRIPT" <<'EOF'
#!/bin/bash
logS() { echo "$(date '+%F %T') $*" >> 20a.mainLog.txt; }
logS "1.$simRep started: jobID=$SLURM_JOB_ID"
echo "Running simulation: simRep=$simRep | jobID=$SLURM_JOB_ID"
source 00b.runParameters.txt
module load StdEnv/2023  gcc/12.3  slim/4.3
slim -l -m \
-d simRep='"'"$simRep"'"' \
-d N=$N -d L=$L -d R=$R -d M=$M \
-d TBUR=$TBUR -d TALLO=$TALLO \
-d NUMINV=$NUMINV \
-d SAMALLO=$SAMALLO -d SAMSEC='"'"$SAMSEC"'"' \
01a.simulation.slim
echo "Done simulation, running python job"
module load python/3.12.4  scipy-stack/2025a
python 01b.treeseq_processing.py $simRep
echo "Done python job, updating jobMap parentJobID to done"
( flock -x 200
awk -F'\t' -v OFS='\t' -v sr="$simRep" '
NR==1 {print; next} { if ($5==sr) $4="done"; print }
' 20b.jobMap.txt > 20b.jobMap.tmp && mv 20b.jobMap.tmp 20b.jobMap.txt
) 200>20b.jobMap.lock
echo "Done updating jobMap parentJobID to done, good job"
logS "1.$simRep ended: jobID=$SLURM_JOB_ID"
EOF
		SIMJOBID=$(sbatch --parsable \
		--dependency=afterok:$SLURM_JOB_ID \
		--mem-per-cpu=64000 --time=23:58:00 \
		--job-name=1.${simRep}.${NUMSIMREPS} \
		--output=10.SLURM_outputs/1_simRep${simRep}.txt \
		--export=simRep=$simRep \
		"$SIMSCRIPT")
		rm "$SIMSCRIPT"
		#update parentJobID of struc rows for the same simRep
		( flock -x 200
			awk -F'\t' -v OFS='\t' -v sr="$simRep" -v sjid="$SIMJOBID" '
				NR==1 {print; next}   # keep header
				{ if ($5==sr) $4=sjid; print }
				' 20b.jobMap.txt > 20b.jobMap.tmp && mv 20b.jobMap.tmp 20b.jobMap.txt
		) 200>20b.jobMap.lock
		logS "1.$simRep scheduled: jobID=$SIMJOBID"
		
		((NJOBSSUB++)) #increase count of scheduled jobs by 1
		#append scheduled job ID to list
		if [[ -z "$SCIDS" ]]; then SCIDS=$SIMJOBID; else SCIDS="$SCIDS:$SIMJOBID"; fi
	done
fi

NSTRUCJOBS=$(tail -n +2 20b.jobMap.txt | wc -l)
logS "0.$NSCHED started: jobID=$SLURM_JOB_ID | STRUCTURE job progress: \
$(awk -F'\t' '$2 == "comp"' 20b.jobMap.txt | wc -l)/$NSTRUCJOBS"

#IF not all status in jobMap are comp
if [ $(awk -F'\t' '$2 != "comp" && NR > 1' 20b.jobMap.txt | wc -l) -gt 0 ]; then

#run until NJOBSSUB=997 or run out of jobMap rows
while IFS=$'\t' read -r arrayID status jobID parentJobID simRep _; do
	if [ "$status" != "comp" ]; then
		
		if [ "$status" != "notSched" ]; then #applies to sched and fail statuses
			#get slurm job status if it exists
			slurm_stat=$(sacct -j "$jobID" --format=State,ExitCode --noheader | head -n 1 | xargs)
			logS "2.$arrayID job error, rescheduling: status=$status | jobStatus=$slurm_stat |" \
				"simRep=$simRep | original jobID=$jobID | parentJobID=$parentJobID"
			mkdir -p 10b.STRUC_fail_logs
			src="10.SLURM_outputs/2_arrayID${arrayID}.txt"
			if [ -f "$src" ]; then mv "$src" "10b.STRUC_fail_logs/2_arrayID${arrayID}.txt"; fi
		fi
		
		#schedule struc job with parent jobID afterok dependency if it's not yet done, pass arrayID
		if [ "$parentJobID" == "done" ]; then DEPJOBID="$SLURM_JOB_ID"
		else DEPJOBID="$SLURM_JOB_ID:$parentJobID"; fi
		SCID=$(sbatch --parsable --dependency=afterok:"$DEPJOBID" \
		--mem-per-cpu=1000 --time=05:58:00 \
		--job-name=2.${arrayID}.${NSTRUCJOBS} \
		--output=10.SLURM_outputs/2_arrayID${arrayID}.txt \
		--export=arrayID=$arrayID \
		02a.STRUCrunner.sh)
		logS "2.$arrayID scheduled: simRep=$simRep | jobID=$SCID"
		
		#update status and jobID of row in jobMap
		( flock -x 200
			awk -F'\t' -v OFS='\t' -v ai="$arrayID" -v jid="$SCID" '
				$1==ai { $2="sched"; $3=jid } {print}
				' 20b.jobMap.txt > 20b.jobMap.tmp && mv 20b.jobMap.tmp 20b.jobMap.txt
		) 200>20b.jobMap.lock
		
		((NJOBSSUB++)) #increase count of scheduled jobs by 1
		#append scheduled job ID to list
		if [[ -z "$SCIDS" ]]; then SCIDS=$SCID; else SCIDS="$SCIDS:$SCID"; fi
		#stop scheduling jobs if at max
		if [ $NJOBSSUB -ge 997 ]; then break; fi
		
	fi #end of if jobStatus != comp
done < <(tail -n +2 20b.jobMap.txt)

#scheduler reschedules itself (still within IF all status != comp)
((NSCHED++))
NEXTSCHEDNAME="0.${NSCHED}_$(awk -F'\t' '$2 == "comp"' 20b.jobMap.txt | wc -l).$NSTRUCJOBS"
SCHEDJOBID=$(sbatch --parsable --dependency=afterany:$SCIDS \
--export=NSCHED=$NSCHED --job-name=$NEXTSCHEDNAME 00a.scheduler.sh)
logS "0.$NSCHED scheduled: jobID=$SCHEDJOBID | numJobsScheduled=$NJOBSSUB"


#ELSE all status in jobMap are comp: schedule R job
else

rm 20b.jobMap.lock

RSCRIPT=$(mktemp)
cat > "$RSCRIPT" <<'EOF'
#!/bin/bash
echo "$(date '+%F %T') 3 started, running R job"
module load StdEnv/2023 r/4.4.0
Rscript 03.STRUCTUREAnalysis.R
echo "$(date '+%F %T') good job R, compressing only 11a.treeSeqs"
tar -c --use-compress-program="pigz -p $SLURM_CPUS_PER_TASK" -f 11a.treeSeqs.tar.gz 11a.treeSeqs
echo "$(date '+%F %T') compressed 11a.treeSeqs, good job 3"
EOF

RJOBID=$(sbatch --parsable \
--mem=8000 --cpus-per-task=4 --time=23:58:00 \
--job-name=3.R --output=20a.mainLog.txt --open-mode=append \
"$RSCRIPT")
logS "3 scheduled: jobID=$RJOBID"

fi