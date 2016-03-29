import sys
import os
import errno
import gzip
import re
from jobTree.src.bioio import fastqRead


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def hardcoded_experiments(a):
    """
    Uses hard coded information about what experiments we expect to see to rename things.
    """
    if "Caroli" in a or "Pahari" in a:
        return a.split("_")[0]
    elif a.startswith("ERR"):
        return a.split("_")[0]
    elif a.startswith("SRR") and "_" in a:
        return a.split("_")[0]
    elif a.startswith("SRR"):
        return a.split(".")[0]
    elif re.match("^[BT]", a).pos == 0:
        return "_".join(a.split("_")[:3])
    # if we get here we need to write more rules
    raise RuntimeError("Error: need to hardcode new experiment names.")


def find_paired_fastqs(source_dir, base_path, files):
    experiment_map = {}
    for x in files:
        experiment = hardcoded_experiments(x)
        p = os.path.join(source_dir, base_path, x)
        if experiment in experiment_map:
            assert not type([]) == type(experiment_map[experiment]), (experiment_map[experiment], p)
            experiment_map[experiment] = [experiment_map[experiment], p]
        else:
            experiment_map[experiment] = p
    return experiment_map


def is_paired_sequencing(bamfile):
    # TODO: this is scary. Should check for unpaired being 0, and number paired == total number
    r = pysam.flagstat(bamfile)
    paired = int(r[5].split()[0])
    if paired != 0:
        return True
    else:
        return False


def fastq_read_size(fastq_path, num_reads=10000):
    sizes = []
    fastq_handle = fastqRead(gzip.open(fastq_path))
    for i in xrange(num_reads):
        name, seq, qual = fastq_handle.next()
        sizes.append(len(seq))
    return 1.0 * sum(sizes) / len(sizes)
