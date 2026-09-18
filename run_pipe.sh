#!/bin/bash

conda activate nextflow

export root_dir="/g"
export work_dir=$PWD
export executor="slurm"
export partition_fast_short="normal_prio"
export partition_slow_long="normal_prio_long"
export partition_slowest_unlimited="normal_prio_unlim"

mkdir -p log/

nextflow -log $PWD/log/nextflow.log run main.nf -resume
