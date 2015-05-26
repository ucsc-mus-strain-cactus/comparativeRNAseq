import sys
import os
import errno
from jobTree.src.bioio import fastqRead

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def find_paired_fastqs(source_dir, base_path, files):
    fastq_files = {}
    base_filenames = {x.split(".")[0] for x in files}
    experiments = {x.split("_")[0] for x in base_filenames}
    for experiment in experiments:
        if experiment in base_filenames:
            # this experiment is unpaired
            path = os.path.join(source_dir, base_path, experiment + ".fastq.gz")
            assert os.path.exists(path)
            fastq_files[experiment] = path
        else:
            # this experiment is paired
            fwd = os.path.join(source_dir, base_path, experiment + "_1.fastq.gz")
            rev = os.path.join(source_dir, base_path, experiment + "_2.fastq.gz")
            assert all([os.path.exists(fwd), os.path.exists(rev)])
            fastq_files[experiment] = [fwd, rev]
    return fastq_files


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
    fastq_handle = fastqRead(fastq_path)
    for i in xrange(num_reads):
        name, seq, qual = fastq_handle.next()
        sizes.append(len(seq))
    return 1.0 * sum(sizes) / len(sizes)