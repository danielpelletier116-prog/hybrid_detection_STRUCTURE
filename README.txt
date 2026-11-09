This README.txt file was generated on 2026-11-09 by Daniel Pelletier

--------------------
GENERAL INFORMATION
--------------------

Title of Dataset:
Scripts and results for the working paper "The program STRUCTURE is good at detecting contemporary hybridization" (working title)

Author Contact Information:
	Name: Daniel Pelletier
	Institution:  Simon Fraser University
	Email: danielpelletier116@gmail.com

Date range of data simulation and analysis:
2025-11-27 to 2026-10-30

Information about funding sources that supported the collection of the data: 
The full pipeline was run on the Fir cluster computer of the Digital Research Alliance of Canada.
Work was supported by funding from NSERC (CGS-M), SFU, and the Maud-Menten Institute.

---------------------------
SHARING/ACCESS INFORMATION
---------------------------

Licenses/restrictions placed on the data:
These data are available under a CC BY 4.0 license <https://creativecommons.org/licenses/by/4.0/> 

Links to publications that cite or use the data:
These scripts and results were generated for the working paper "The program STRUCTURE is good at detecting contemporary hybridization",
to be submitted for review.

Links/relationships to ancillary data sets or software packages:
Pipeline is designed to be run with the SLURM workload manager: https://slurm.schedmd.com/overview.html
Simulations and processing of simulated data use the following software and packages:
	-SLiM v4.3: https://messerlab.org/slim/
	-python v3.12.4: https://www.python.org/downloads/release/python-3124/
	-tskit v1.0.2, msprime v1.4.1, and pyslim v1.1.2 python libraries: https://tskit.dev/
	-STRUCTURE v2.3.4: https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html
	-R v4.4.0: https://cran.r-project.org/bin/windows/base/old/4.4.0/
	-hierfstat R package: https://cran.r-project.org/web/packages/hierfstat/index.html

Recommended citation for this dataset: 
To be updated upon dataset archiving in a public repository.

---------------------
DATA & FILE OVERVIEW
---------------------

This dataset includes 10 scripts and parameter text files, and 1 results file.
The results file includes all seeds necessary to reproduce the full run.

The purpose of the scripts is to simulate genetic data under a variety of evolutionary parameters,
and test the performance of the software STRUCTURE at detecting hybridization from that genetic data.
This is designed to be run by submitting jobs to a cluster computer through the SLURM scheduler.
Briefly, the main scheduling script reads parameters from the runParameters file, and schedules the
appropriate number of replicate simulation jobs. It then schedules STRUCTURE jobs to analyse the genetic
data output from each simulation. Due to a maximum of 1000 schedule-able jobs at a time, it only schedules
up to 997 jobs, then re-schedules itself to be run once all those jobs are done, repeating until all
required jobs have been run successfully. Afterwards, it schedules an R job which extracts and
processes outputs of all the STRUCTURE jobs. The scheduler doesn't run the final R analysis job, this one
is meant to be run locally by the user.

Note that terminology used in the paper and in scripts is different in two cases (see paper for details):
	1: The populations and individuals referred to as "native" in the paper are referred to as "recipient" in
		scripts, abbreviated to "rec".
	2: The invasion scenarios referred to as "Invasion" and "No-invasion" in the paper are referred to as "Allo"
		(short for Allopatric) and "Sec" (short for Secondary Contact), respectively, in the scripts.

File and folder naming:
All files and folders (including ones produced by running the full pipeline, which are not included here)
start with a 2-digit number, and optionally a letter, before the dot.
The 1st digit indicates the type of object:
	0: a script or parameter text file, required in order to run the pipeline
	1: a folder containing results files
	2: a standalone results file
The 2nd digit indicates which part of the pipeline the object is relevant to:
	0: the initialization or top-level organization of the pipeline
	1: the evolutionary simulation part of the pipeline
	2: running the STRUCTURE analysis software
	3: processing STRUCTURE outputs, turning them into analysable results
	4: evaluating the performance of the STRUCTURE software
The optional letter is just to differentiate multiple files or folders that have the same 2-digit code.

Script files:
00a.scheduler.sh:
	The central scheduler script. To start the full pipeline, run: sbatch 00a.scheduler.sh, in
	a directory containing all the other files.
00b.runParameters.txt:
	Contains the parameter values used for simulations and analysis. The values in this file are the ones used in
	the full run that produced the results (below); they can be modified to different values before starting the full run.
00c.mapper.R:
	00a starts by running this, to produce a table (jobMap) containing one row for each future
	STRUCTURE job to be run.
01a.simulation.slim:
	The evolutionary simulation script, written for SLiM v.4.3 (not compatible with newer SLiM releases).
	This is run for each simulation replicate, and outputs tree-sequence files.
01b.treeseq_processing.py:
	This is run right after every 01a run. Uses the tskit group of python packages to
	process the tree-sequences produced by 01a, simulate neutral mutations on them, and produce SNP files for
	sampled simulated individuals for each scenario. Also writes Fst results from this simulation replicate
	to the 21b.Fst.txt file.
02a.STRUCrunner.sh:
	This script is run once for each STRUCTURE analysis run. 00a gives it a specific ID number,
	which this script uses to source the specific parameters for one STRUCTURE run.
02b.mainParams.txt:
	One of the two files of general parameters that STRUCTURE needs.
02c.extraParams.txt:
	The second file of general parameters that STRUCTURE needs.
03.STRUCTUREAnalysis.R:
	This is run once all STRUCTURE jobs are done. It extracts and processes the outputs from all STRUCTURE jobs,
	saving them in a large file named 23.fullResults_<date>.Rdata, with the date when the file was written in
	YYYY-MM-DD format.
04.resultsAnalysis.R:
	This final analysis script is not scheduled in the main pipeline; instead, it's meant to be run locally
	on the user's computer, only requiring the 23.fullResults_<date>.Rdata and 00b.runParameters.txt files.
	It produces all plots in folders 14a.plots/ and 14b.QvQplots/, and prints out results for tables.

Results file:
23.fullResults_2026-10-30.Rdata:
	This file contains all results needed to produce the plots and tables in the paper. When read into R, it contains 8 objects:
	
	empFst: Contains the empirical-POV Fst measures, taken from the SNP files. Table with 1 row for each of the 32
		unique combination of simulation parameters, and 1 column per simulation replicate (plus 5 columns
		for the simulation parameters of each row).
	fData: Contains the pre-invasion measures of Fst, taken directly from the tree-sequences. Table with 3
		rows per simulation replicate; for each replicate, 1 row each for: Fst in panmictic-scenario, Fst in
		substructured-scenario, and Fst between the 3 subpopulations in substructured-scenario. 1 column per
		sampling time.
	indData: Large list of tables, with 1 table for each replicate simulation of each unique combination of evolutionary
		parameters. Each table has 1 row per sampled simulated individual, and contains information on their true (simulated)
		and STRUCTURE-inferred population ancestry.
	pRes: Table with 1 row for each of the 32 unique simulation parameter combinations, summarizing some of the results
		across replicate simulations and STRUCTURE runs.
	selT: Small table with the selected times (post-divergence) at which invasion was initiated (based on reaching target
		Fst levels), for each simulation replicate and each unique combination of parameters.
	tRes: Large table with 1 row for each of the 96000 STRUCTURE runs performed, containing information on the
		specific simulation parameters, STRUCTURE run parameters, and results for each.
	pySeeds: Vector containing the seeds for the python simulation step, ordered by simulation replicate number (1-50).
	simSeeds: Vector containing the seeds for SLiM simulations, ordered by simulation replicate number (1-50).

Reproducibility the full pipeline run that generated the results presented in the paper:
Seeds for SLiM simulations are in the "simSeeds" vector object, ordered by simulation replicate number (1-50), which can
	be read in R by loading the	23.fullResults_2026-10-30.Rdata file. 
Seeds for python simulations are in the "pySeeds" vector object, ordered by simulation replicate number (1-50), which can
	be read in R by loading the	23.fullResults_2026-10-30.Rdata file. 
Seeds for each STRUCTURE run are in the "seed" column of the tRes table, which can be read in R by loading the
	23.fullResults_2026-10-30.Rdata file.
