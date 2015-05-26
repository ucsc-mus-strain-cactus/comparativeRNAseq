#!/usr/bin/bash
batchSystem="parasol"
jobTree=".jobTree"
maxThreads="15"
log="kallisto.log"

export PATH=./sonLib/bin:./jobTree/bin:${PATH}
export PYTHONPATH=./:${PYTHONPATH}

if [ -d ${jobTree} ]; then
    rm -rf ${jobTree}
fi

umask 0002

python bin/kallisto_from_fastq.py --jobTree ${jobTree} --batchSystem ${batchSystem} --logLevel DEBUG --maxThreads ${maxThreads} &> ${log}
