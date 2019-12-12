from __future__ import print_function
import os
import re
import gzip
import itertools
import argparse
import time
import logging

__author__ = 'Martin Aryee'

logger = logging.getLogger('root')

#args = {'out_dir':'/PHShome/ma695/tmp', 'min_reads':10}
#base = '/data/joung/sequencing_bcl/131007_M01326_0075_000000000-A6B33/Data/Intensities/BaseCalls'
#args['read1'] = os.path.join(base, 'Undetermined_S0_L001_R1_001.fastq.gz')
#args['read2'] = os.path.join(base, 'Undetermined_S0_L001_R2_001.fastq.gz')
#args['index1'] = os.path.join(base, 'Undetermined_S0_L001_I1_001.fastq.gz')
#args['index2'] = os.path.join(base, 'Undetermined_S0_L001_I2_001.fastq.gz')

def fq(file):
    if re.search('.gz$', file):
        fastq = gzip.open(file, 'rb')
    else:
        fastq = open(file, 'r')
    with fastq as f:
        while True:
            l1 = f.readline()
            if not l1:
                break
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            yield [l1, l2, l3, l4]


def get_sample_id(i1, i2, sample_names, max_hamming):
    seq1 = i1[1]
    seq2 = i2[1]

    #	Why originally 1:8? Seems that would ignore the first character?
    #sample_barcode = seq1[1:8] + seq2[1:8]
    #if sample_barcode in sample_names:
    #    return sample_names[sample_barcode]
    #else:
    #    return sample_barcode

    sample_barcode = seq1[0:8] + seq2[0:8]

    hamming_dist = {}
    for barcode in sample_names.keys():
        i5dist = 0
        i7dist = 0
        #	Note ...
        #	>>> range(8)
        #	[0, 1, 2, 3, 4, 5, 6, 7]
        for character in range(8):
            if barcode[character] != sample_barcode[character]:
                i7dist +=1
            if barcode[character + 7] != sample_barcode[character + 7]:
                i5dist +=1
        if i5dist <= max_hamming and i7dist <= max_hamming:
            hamming_dist[barcode] = i5dist + i7dist
    if len(hamming_dist) != 0:
        m = min(hamming_dist.values())
        if hamming_dist.values().count(m) == 1:
            min_barcode = [barcode for barcode, distance in hamming_dist.items() if hamming_dist[barcode] == m]
            return sample_names[min_barcode[0]]
        else:
            return 'multiplematches'
    else:
        return 'hdist_toofar'

def demultiplex(read1, read2, index1, index2, sample_barcodes, out_dir, min_reads, max_hamming):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if type(sample_barcodes) != dict:
        sample_names = {}
        if not sample_barcodes==None:
            for line in open(sample_barcodes, 'r'):
                fields = line.strip().split('\t')
                if len(fields)==3:
                    sampleid, i7barcode, i5barcode = fields
                    sample_names[i7barcode + i5barcode] = sampleid
    else:
        sample_names = sample_barcodes

    outfiles_r1 = {}
    outfiles_r2 = {}
    outfiles_i1 = {}
    outfiles_i2 = {}

    total_count = 0
    count = {}
    buffer_r1 = {}
    buffer_r2 = {}
    buffer_i1 = {}
    buffer_i2 = {}

    #it = itertools.izip(fq(args['read1']), fq(args['read2']), fq(args['index1']), fq(args['index2']))
    #for r1,r2,i1,i2 in itertools.islice(it, 0, 100):
    start = time.time()
    for r1,r2,i1,i2 in itertools.izip(fq(read1), fq(read2), fq(index1), fq(index2)):
        total_count += 1
        if total_count % 1000000 == 0:
            logger.info("Processed %d reads in %.1f minutes.", total_count, (time.time()-start)/60)
        sample_id = get_sample_id(i1, i2, sample_names, max_hamming)

        # Increment read count and create output buffers if this is a new sample barcode
        if not count.has_key(sample_id):
            count[sample_id] = 0
            buffer_r1[sample_id] = []
            buffer_r2[sample_id] = []
            buffer_i1[sample_id] = []
            buffer_i2[sample_id] = []
        count[sample_id] += 1

        # Write the reads to file or to sample-specific output buffers if we haven't
        # reached min_reads
        # The buffers are used to avoid writing fastqs for sample barcodes with very few reads.
        if count[sample_id] < min_reads:
            buffer_r1[sample_id].append(r1)
            buffer_r2[sample_id].append(r2)
            buffer_i1[sample_id].append(i1)
            buffer_i2[sample_id].append(i2)
        elif count[sample_id] == min_reads:
            outfiles_r1[sample_id] = open(os.path.join(out_dir, '%s.r1.fastq' % sample_id), 'w')
            outfiles_r2[sample_id] = open(os.path.join(out_dir, '%s.r2.fastq' % sample_id), 'w')
            outfiles_i1[sample_id] = open(os.path.join(out_dir, '%s.i1.fastq' % sample_id), 'w')
            outfiles_i2[sample_id] = open(os.path.join(out_dir, '%s.i2.fastq' % sample_id), 'w')
            # Spill the buffers to sample-specific fastqs
            for record in buffer_r1[sample_id] + r1:
                outfiles_r1[sample_id].write(''.join(record))
            for record in buffer_r2[sample_id] + r2:
                outfiles_r2[sample_id].write(''.join(record))
            for record in buffer_i1[sample_id] + i1:
                outfiles_i1[sample_id].write(''.join(record))
            for record in buffer_i2[sample_id] + i2:
                outfiles_i2[sample_id].write(''.join(record))
            del buffer_r1[sample_id]
            del buffer_r2[sample_id]
            del buffer_i1[sample_id]
            del buffer_i2[sample_id]
        else:
            for line in r1:
                print (line, file=outfiles_r1[sample_id], end="")
            for line in r2:
                print (line, file=outfiles_r2[sample_id], end="")
            for line in i1:
                print (line, file=outfiles_i1[sample_id], end="")
            for line in i2:
                print (line, file=outfiles_i2[sample_id], end="")

    # Write remaining buffered reads to a single fastq.
    # (These reads correspond to barcodes that were seen less than min_reads times)
    undetermined_r1 = open(os.path.join(out_dir, 'undetermined.r1.fastq'), 'w')
    undetermined_r2 = open(os.path.join(out_dir, 'undetermined.r2.fastq'), 'w')
    undetermined_i1 = open(os.path.join(out_dir, 'undetermined.i1.fastq'), 'w')
    undetermined_i2 = open(os.path.join(out_dir, 'undetermined.i2.fastq'), 'w')
    for sample_id in buffer_r1.keys():
        for record in buffer_r1[sample_id]:
            undetermined_r1.write(''.join(record))
        for record in buffer_r2[sample_id]:
            undetermined_r2.write(''.join(record))
        for record in buffer_i1[sample_id]:
            undetermined_i1.write(''.join(record))
        for record in buffer_i2[sample_id]:
            undetermined_i2.write(''.join(record))

    seq_counts = open(os.path.join(out_dir, 'seqcounts.txt'), 'w')
    samples = []
    strings = []
    for sample in count.keys():
        if sample.isdigit():
            samples.append(int(sample))
        else:
            strings.append(sample)
    samples.sort()
    samples.extend(strings)
    samples.append('undetermined')
    # Close files
    for sample_id in outfiles_r1:
        outfiles_r1[sample_id].close()
        outfiles_r2[sample_id].close()
        outfiles_i1[sample_id].close()
        outfiles_i2[sample_id].close()
    undetermined_r1.close()
    undetermined_r2.close()
    undetermined_i1.close()
    undetermined_i2.close()

    total_reads = 0
    for sample in samples:
        if os.path.exists("%s.r1.fastq" % str(sample)):
            reads = len(open("%s.r1.fastq" % str(sample)).readlines(  ))/4
            seq_counts.write("%s had %d reads \n" % (sample, reads))
            total_reads += reads
    seq_counts.write("Total Reads:%.0f" % int(total_reads))
    num_fastqs = len([v for k,v in count.iteritems() if v>=min_reads])
    logger.info('Wrote FASTQs for the %d sample barcodes out of %d with at least %d reads.', num_fastqs, len(count), min_reads)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read1', required=True)
    parser.add_argument('--read2', required=True)
    parser.add_argument('--index1', required=True)
    parser.add_argument('--index2', required=True)
    parser.add_argument('--min_reads', type=int, default=10000)
    parser.add_argument('--max_hamming', type=int, default=0)
    parser.add_argument('--sample_barcodes')
    parser.add_argument('--out_dir', default='.')
    args = vars(parser.parse_args())

    demultiplex(args['read1'], args['read2'], args['index1'], args['index2'], args['sample_barcodes'], args['out_dir'], min_reads=args['min_reads'], max_hamming=args['max_hamming'])

if __name__ == '__main__':
    main()
