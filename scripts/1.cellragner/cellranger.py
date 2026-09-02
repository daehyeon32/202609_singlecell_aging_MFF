# python3 crc.cellranger.01.run.py

import sys
import argparse
import os
import glob
import re
import subprocess
import matplotlib
import pandas as pd
matplotlib.use('Agg')

## macrogen / nicem / DNAlink / rokit / macrogen_multiomeseq
# (macrogen) /BiO2/Store/UNIST-CRC-YUHS-2023-11/1_Fresh_tissue/01_scRNAseq/2019_multiomics_scRNA_macrogen_3case
# (nicem) /BiO2/Store/UNIST-CRC-YUHS-2023-11/1_Fresh_tissue/01_scRNAseq/2019_multiomics_scRNA_nicem_47case
# (DNAlink) /BiO2/Store/UNIST-CRC-YUHS-2023-11/1_Fresh_tissue/01_scRNAseq/2021_multiomics_scRNA_DNAlink_5case
# (rokit) /BiO2/Store/UNIST-CRC-YUHS-2023-11/1_Fresh_tissue/01_scRNAseq/2021_multiomics_scRNA_rokit_30case
# (macrogen_multiomeseq) /BiO2/Store/UNIST-CRC-YUHS-2023-11/1_Fresh_tissue/05_Multiomeseq

dir_input = '/BiO/Store/UNIST-Aging-MFF-2026-09/Tabula_Muris_Senis_FASTQ/'
dir_output = '/BiO/Live/dleogus32/202609Aging_MFF/Analysis/1.cellranger/'

if not os.path.exists(dir_output):
	os.makedirs(dir_output)
if not os.path.exists(dir_output + "/sh"):
	os.makedirs(dir_output + "/sh")
if not os.path.exists(dir_output + "/stdeo"):
	os.makedirs(dir_output + "/stdeo")
		
#------------------------------------------------------------------------------------------------------------------------------------------

dir_cellranger = '/BiO/Share/Tools/cellranger-9.0.1/cellranger'
dir_ref = '/BiO/Share/Database/cellranger_reference/refdata-cellranger-mm10-3.0.0/'

#------------------------------------------------------------------------------------------------------------------------------------------

tissues = ["Heart", "Limb_Muscle"]

# Heart/*, Limb_Muscle/* 아래의 샘플 폴더 수집
fastqs_dirs = []

for tissue in tissues:
    tissue_dir = os.path.join(dir_input, tissue)
    sample_dirs = [ d for d in glob.glob(os.path.join(tissue_dir, "*")) if os.path.isdir(d)]
    fastqs_dirs.extend(sample_dirs)

fastqs_dirs.sort()
#fastqs_dirs = fastqs_dirs[0:1]

def get_fastq_id(fastqs_dir):
    r1_files = glob.glob(os.path.join(fastqs_dir, "*_R1_*.fastq.gz"))

    if not r1_files:
        raise RuntimeError(f"R1 FASTQ가 없습니다: {fastqs_dir}")

    fastq_ids = set()

    for r1_file in r1_files:
        filename = os.path.basename(r1_file)
        match = re.match(r"^(.+)_S\d+_L\d{3}_R1_\d{3}\.fastq\.gz$", filename)

        if match is None:
            raise RuntimeError(f"FASTQ 이름을 해석할 수 없습니다: {filename}")

        fastq_ids.add(match.group(1))

    if len(fastq_ids) != 1:
        raise RuntimeError(f"여러 fastq_id가 발견되었습니다: {fastqs_dir} -> {fastq_ids}")

    return fastq_ids.pop()


job_ids = []

for fastqs_dir in fastqs_dirs:
    tissue = os.path.basename(os.path.dirname(fastqs_dir))
    fastq_id = get_fastq_id(fastqs_dir)
    run_id = f"{tissue}__{fastq_id}"

    print(f"tissue={tissue}, fastq_id={fastq_id}, run_id={run_id}")

    sh_cellranger = os.path.join(dir_output, "sh", f"cellranger_{run_id}.sh")

    order = (
        f"{dir_cellranger} count "
        f"--id={run_id} "
        f"--transcriptome={dir_ref} "
        f"--fastqs={fastqs_dir} "
        f"--sample={fastq_id} "
        f"--create-bam=true "
        f"--localcores=8 "
        f"--localmem=100"
    )

    with open(sh_cellranger, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(order + "\n")

    dependency = ""
    if len(job_ids) >= 10:
        dependency = f"--dependency=afterok:{job_ids[-10]} "

    s_cellranger = (
        f"sbatch "
        f"{dependency}"
        f"--job-name=cellranger_{run_id} "
        f"--cpus-per-task=8 "
        f"--mem=100G "
        f"--error={dir_output}/stdeo/%x.e%A "
        f"--output={dir_output}/stdeo/%x.o%A "
        f"--chdir={dir_output} "
        f"{sh_cellranger}"
    )

    job_id = str(
        subprocess.check_output(s_cellranger, shell=True),
        "utf-8"
    ).split()[-1]

    job_ids.append(job_id)
