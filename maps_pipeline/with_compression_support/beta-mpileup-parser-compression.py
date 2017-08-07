#! /usr/bin/env python2.6

import os, sys, math, datetime, gc, time
import threading, multiprocessing
import argparse
import subprocess
import gzip
import csv
import re

from collections import Counter

INSERT_RE = re.compile('\+[0-9]+[ACGTNacgtn]+')

#Edited 08-06-2017 by Christian Schudoma to make usable

#Edited 11-24-2015 by Tyson Howell to output to stdout so that it can be compressed or piped. Also added an option to provide an input file containing the output of `wc -l` run on the input (recommended to generate this when writing out the mpileup in the previous step) to reduce I/O. Intermediate files are also compressed.

#Comai Lab, Ucdavis Genome Center
#Meric Lieberman, 2012
# This work is the property of UC Davis Genome Center - Comai Lab

# Use at your own risk.
# We cannot provide support.
# All information obtained/inferred with this script is without any
# implied warranty of fitness for any purpose or use whatsoever.
#------------------------------------------------------------------------------

#Part 2: mpileup-parser.py
#
#This program parses a mpileup file to a simplified format to be used with our MAPS
#mutation and genotyping package.
#
#INPUT:
#This program takes an mpileup.txt file as input
#
#OUTPUT:
#This program outputs a parsed mpileup.txt file.
#
#NOTE:
#This program reads the entire file into memory before parsing, so it is recommended not to be run on systems with limited memory. It typically requires 1.5 times the size of the mpileup file in RAM to run. Many machines will not have this, and in this case it is recommended to break the mpileup file into smaller chunks to be processed separately. (i. e. by chromosome or scaffold)
#If being used with the MAPS package: the first step in MAPS is also threaded so it may be best to leave the chunks separate and process them individually in MAPS as well. Results can be combined at the end without compromising the results.
#This program is threaded, so it can be used with the -t flag to specify the number of cores to be used
#
#PARAMETERS:
#1. REQUIRED, default value in []:
#-f or --mpileup_file, The input mpileup file. [required]
#2. OPTIONAL:
#-t or --thread, Number of cores to be used while processing. [1]

#split file into chunks
def splitter(l, n):
    i = 0
    chunk = l[:n]
    while chunk:
        yield chunk
        i += n
        chunk = l[i:i+n]

# this is translated from original script,
# in which
#  - deletion handling is missing
#  - insertions are only counted once
def scanChanges(changes, validChars):
    c_changes, c_inserts = Counter(), Counter()
    p = 0
    last = None
    while p < len(changes):
        if changes[p] in validChars:
            c_changes[changes[p].upper()] += 1
            last = changes[p].upper()
        elif changes[p] == "^":
            p += 1
        elif changes[p] == '+' or changes[p] == '-':
            c_changes[last] -= 1
            inslen = re.match('[0-9]+', changes[p + 1:]).group()
            p += int(inslen) + len(inslen)
            if changes[p] == '+':
                c_inserts['.{}'.format(changes[p + 1:p + 1 + int(inslen) + len(inslen)]).upper()] += 1
        p += 1
    return c_changes, c_inserts



#parse one pileup line
#based on parsing script from Joeseph Fass

class MyThread (multiprocessing.Process):

    def __init__ (self, filen, res, startpos, endpos):
        self.filen = filen
        self.res = res
        self.startpos = startpos
        self.endpos = endpos
        multiprocessing.Process.__init__ (self)

    def run (self):
        sys.stderr.write("{} {}\n".format(self.startpos, self.endpos))
        counter, ctgood = 1, 0
        all_pos = self.endpos - self.startpos
        with gzip.open('temp-parse-{}.txt.Z'.format(self.res), 'wt') as gzip_out, gzip.open(self.filen, 'rt') as gzip_in:
            next(gzip_in)
            for lcount, row in enumerate(gzip_in, start=1):
                if lcount < self.startpos:
                    continue
                break
            mpreader = csv.reader(gzip_in, delimiter='\t')
            for lcount, row in enumerate(mpreader, start=lcount):
                if lcount > self.endpos:
                    break
                if ctgood % 100008 == 0:
                    sys.stderr.write("{} {}/{}\n".format(self.res, ctgood, all_pos))
                ctgood += 1
                result, refbase = row[:3], row[2]
                samples = list(splitter(row[3:], 3))
                for depth, changes, qualities in samples:
                    try:
                        # TODO: make sure that quals are phred33!
                        mean_SQ = sum(map(lambda x: ord(x)-33, list(qualities))) / len(qualities)
                    except ZeroDivisionError:
                        mean_SQ = 0.0
                    if mean_SQ < 20.0:
                        result.extend(['.', '.', '.', '.'])
                    else:
                        inserts, quals = Counter(), Counter()
                        valid = set(['a','A','c','C','g','G','t','T','.',',','*'])
                        quals, inserts = scanChanges(changes, valid)
                        total_HQ = sum(quals.values())
                        Aa_HQ, Tt_HQ, Cc_HQ, Gg_HQ = quals['A'], quals['T'], quals['C'], quals['G']
                        match_HQ = quals['.'] + quals[',']
                        dels = quals['*']
                        scan = list()

                        if total_HQ == 0 and not inserts:
                            result.extend(['.', '.', '.', '.'])
                        else:
                            if inserts:
                                inname, incount = sorted(inserts.items(), key=lambda x:x[1], reverse=True)[0]
                                inper = incount / float(incount + total_HQ) * 100
                                delper = dels / float(incount + total_HQ) * 100
                                oin = sum(inserts.values())
                                total_HQ += oin
                                scan.append([inname, inper])
                            else:
                                inname, inper, oin = '.', '.', '.'
                                delper = dels / float(total_HQ) * 100

                            aPer = Aa_HQ / float(total_HQ) * 100
                            tPer = Tt_HQ / float(total_HQ) * 100
                            cPer = Cc_HQ / float(total_HQ) * 100
                            gPer = Gg_HQ / float(total_HQ) * 100
                            matchPer = match_HQ / float(total_HQ) * 100
                            scan.extend([['A',aPer],['T',tPer],['C',cPer],['G',gPer],['*',delper]])

                            for w in scan:
                               if w[0] == refbase:
                                  w[1] += matchPer
                                  break

                            for w in sorted(scan, key=lambda x:x[1], reverse=True)[:3]:
                               if w[1] == 0.0:
                                  result.append('.')
                               else:
                                  result.append('{}_{}'.format(w[0], round(w[1], 2)))
                            result.append(str(total_HQ))

            gzip_out.write('{}\n'.format('\t'.join(result)))


# Uses wc to get te number of lines in the file
def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])


if __name__ == '__main__':

    start = time.time()

    ap = argparse.ArgumentParser()
    ap.add_argument('--thread', '-t', dest='nthreads', type=int, default=1, help='How many threads to use during processing (DEFAULT == 1)')
    ap.add_argument('--mpileup_file', '-f', help='Input mpileup file')
    ap.add_argument('--length', '-l', dest='length_file', help='File containing length in lines of input mpileup')
    args = ap.parse_args()

    numThreads = args.nthreads

    t1 =  datetime.datetime.now()

    try:
       filename = args.mpileup_file
       f = gzip.open(filename, 'rt')
    except:
       sys.stderr.write("Could not open input. Please check your command line parameters with -h or --help.\n")
       sys.exit(1)

    flen = 0
    if args.length_file:
        try:
            with open(args.length_file) as line_file:
                flen = int(next(line_file).strip().split()[0])
        except:
            sys.stderr.write("Could not open linecount file. Falling back to counting lines myself.\n")

    if flen == 0:
        flen = file_len(filename)

    cutnum = (flen-1)/numThreads+1
    cutset = [[i * cutnum + 1, min((i + 1) * cutnum, flen)] for i in range(numThreads)]

    mpreader = csv.reader(f, delimiter='\t')
    try:
        header = next(mpreader)
    except:
        sys.exit("Empty pileup file")


    h2 = list(splitter(header[3:],3))
    newhead = header[:3]
    libs = []
    for lab in h2:
        lname = ('-'.join(lab[0].split('-')[1:]))
        libs.append(lname)
        newhead += ["SNP1-{}".format(lname),"SNP2-{}".format(lname),"SNP3-{}".format(lname),"Cov-{}".format(lname)]

    sys.stdout.write('\t'.join(newhead)+'\n')

    threads = list()
    cat = ["zcat"]
    sys.stderr.write("{}\n".format(cutset))
    for counter, x in enumerate(cutset, start=1):
        sys.stderr.write("{}\n".format(counter))
        threads.append(MyThread(filename, counter, x[0], x[1]))
        cat.append("temp-parse-{}.txt.Z".format(counter))
        threads[-1].start()

    for counter, thread in enumerate(threads, start=1):
        thread.join()
        sys.stderr.write("{} joined\n".format(counter))
    sys.stderr.write('{}\n'.format(cat))
    os.system(' '.join(cat))
    os.system("rm -f temp-parse-*")

    now = time.time()-start
    sys.stderr.write("{} {}".format(int(now/60), int(now%60)))
