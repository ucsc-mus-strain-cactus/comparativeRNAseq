#!/usr/bin/bash
batchSystem="singleMachine"
jobTree=".STARjobTree"
maxThreads="10"
log="STAR.log"
defaultMemory="8589934592"

export PATH=./sonLib/bin:./jobTree/bin:${PATH}
export PYTHONPATH=./:${PYTHONPATH}

if [ -d ${jobTree} ]; then
    rm -rf ${jobTree}
fi

umask 0002

python bin/STAR_from_fastq.py --jobTree ${jobTree} --batchSystem ${batchSystem} --defaultMemory ${defaultMemory} \
--logLevel DEBUG --maxThreads ${maxThreads} &> ${log}
