#!/usr/bin/env python3
import os
import shlex
from argparse import ArgumentParser
from subprocess import PIPE, CalledProcessError, Popen
import shutil


def run_command(cmd):
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    out = stdout.decode("utf-8")
    err = stderr.decode("utf-8")
    if p.returncode != 0:
        print('STDERR from called program: {}'.format(stderr))
        print('STDOUT from called program: {}'.format(stdout))
        raise CalledProcessError(p.returncode, cmd)
    else:
        return [out, err]


def genome_size_bbtools(forward):
    # cmd = "kmercountexact.sh in={} in2={} peaks=peak.log".format(forward, reverse)
    cmd = "kmercountexact.sh in={} peaks=peak.log".format(forward)
    peak_log = "peak.log"
    cmd_to_run = shlex.split(cmd)
    # Remove peak_log if existed.
    if os.path.exists(peak_log):
        os.unlink(peak_log)
    p = Popen(cmd_to_run, stderr=PIPE, stdout=PIPE)
    _, stderr = p.communicate()
    gs = cov = 0
    if stderr != 0:
        with open(peak_log, "r") as fh:
            handle = fh.readlines()
            for line in handle:
                if "#haploid_genome_size" in line:
                    gs = line.split("\t")[1].strip()
                if "#haploid_fold_coverage" in line:
                    cov = line.split("\t")[1].strip()
                    break
        os.unlink(peak_log)
        # with open("{}.tsv".format(sample_id), "w") as oh:
        #     oh.writelines("{},{}".format(gs, cov))
    return [gs, float(cov)]


def genome_size_kmc(fastq, minkc=3):
    # kmc -m$ram_int -sm -n256 -ci$minkc -k25 -t$cpus \Q$R1\E kmc $dirname
    tmp = "kmc_tmp"
    if not os.path.exists(tmp):
        os.makedirs(tmp)

    cmd = "kmc -sm -n256 -ci{} -k25 -t4 {} kmc kmc_tmp 2>&1".format(
        minkc, fastq)
    cmd_to_run = shlex.split(cmd)
    stdout, _ = run_command(cmd_to_run)
    gs = 0
    lines = stdout.split("\n")
    for line in lines:
        if "unique counted" in line:
            gs = line.split(":")[1].strip()
    shutil.rmtree(tmp)
    os.unlink("kmc.kmc_pre")
    os.unlink("kmc.kmc_suf")
    return int(gs)


def genome_size_mash(fastq, mincov=25, minkc=3):
    cmd = "mash sketch -o tmp -k 32 -p 4 -m {} -c {} -r {}".format(
        minkc, mincov, fastq)
    cmd_to_run = shlex.split(cmd)
    _, stderr = run_command(cmd_to_run)
    lines = stderr.split("\n")
    genome_size = lines[1].replace("Estimated genome size: ", "")
    gs = int(float(genome_size))
    os.unlink("tmp.msh")
    return gs


def estimate_coverage_mash(fastq, genome_size):
    cmd = "seqkit stat {}".format(fastq)
    cmd_to_run = shlex.split(cmd)
    stdout, _ = run_command(cmd_to_run)
    lines = stdout.split("\n")
    bp = lines[1].split()[4].replace(",", "")
    est_cov = int(bp) * 2 / genome_size
    return est_cov


def find_refseq(forward, reverse, refseq, temp="tmp"):
    if forward.endswith("gz"):
        extr = "zcat"
    else:
        extr = "cat"
    # Sketching read
    sketch_cmd = "{} \"{}\" \"{}\" | mash sketch -o {} -m 2 -r -".format(
        extr, forward, reverse, temp)
    sketch_cmd_to_run = shlex.split(sketch_cmd)
    run_command(sketch_cmd_to_run)
    # Estimate distance between refseq and read
    dist_cmd = "mash dist {} {}.msh | sort -gk3 | head -n1".format(
        refseq, temp)
    dist_cmd_to_run = shlex.split(dist_cmd)
    out, _ = run_command(dist_cmd_to_run)
    lines = out.split("\n")
    refseq = lines.split()[0]
    return refseq


def parse_args():
    parser = ArgumentParser(
        description="Estimate genome size and closely related refseq")
    parser.add_argument(
        '--R1', help="Read (Forward or Reverse or SE)", required=True)
    parser.add_argument('--R2', help="Reverse read", required=False)
    parser.add_argument('--refseq', help="RefSeq Sketch Database (msh)", required=False,
                        default="/home/ubuntu/data2/mash/refseq.k21s1000.042018.msh")
    parser.add_argument(
        '--method', help="Estimate genome size and coverage methods [mash | bbtools | kmc]", required=False, default="mash")
    return parser.parse_args()


def main():
    args = parse_args()
    forward = args.R1
    # reverse = args.R2
    # refseq = args.refseq
    method = args.method

    if (method == "mash"):
        gs = genome_size_mash(forward)
        cv = estimate_coverage_mash(forward, gs)
    elif (method == "kmc"):
        gs = genome_size_kmc(forward)
        cv = estimate_coverage_mash(forward, gs)
    elif (method == "bbtools"):
        gs, cv = genome_size_bbtools(forward)
    else:
        print("Hmm! Please choose a method for estimating the genome size and depth coverage, either mash or bbtoolsor kmc")

    print("{}\t{:.1f}\t{}".format(gs, cv, method))


if __name__ == "__main__":
    main()
