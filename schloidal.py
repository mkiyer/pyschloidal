"""
schloidal
"""
import logging
import argparse

import pysam
import numpy as np
import h5py


def get_refs_from_bam(bamfile):
    f = pysam.Samfile(bamfile, "rb")
    refs = list(zip(f.references, f.lengths))
    f.close()
    return refs


def get_refs_from_sam(samfile):
    f = pysam.Samfile(samfile, "r")
    refs = list(zip(f.references, f.lengths))
    f.close()
    return refs


class GenomeTrack:
    # constants
    DTYPE = np.uint64
    COMPRESSION = 'lzf'

    # attribute names
    REFS = 'refs'
    NUM_FRAGS = 'num_frags'

    def __init__(self, filename, refs=None, mode='a'):
        if h5py.is_hdf5(filename):
            # open existing file
            self.h5f = h5py.File(filename, mode)
        else:
            # create new file
            self.h5f = h5py.File(filename, mode)
            # set attributes
            self.h5f.attrs[GenomeTrack.REFS] = refs
            self.h5f.attrs[GenomeTrack.NUM_FRAGS] = 0
            # create datasets
            for ref, length in refs:
                self.h5f.create_dataset(ref, (2, length,),
                                        dtype=GenomeTrack.DTYPE,
                                        chunks=True,
                                        compression=GenomeTrack.COMPRESSION,
                                        shuffle=True,
                                        fillvalue=0)

    def close(self):
        self.h5f.flush()
        self.h5f.close()


    def add_alignments(self, pysamiter):
        for a in pysamiter:
            print(a.reference_start, a.reference_end, a.cigartuples)


def add_vector_track(parser, options):
    tf = open_trackfactory(parser, options)
    t = tf.create_track(options.name, VectorTrack,
                        pe=options.pe,
                        strand=options.strand,
                        allele=options.allele)
    logging.info("added %s '%s' to trackfactory %s" %
                 (VectorTrack.__name__, options.name,
                  options.file))
    datafile = check_datafile(parser, options)
    if datafile is not None:
        logging.info("inserting data file %s (type=%s)" %
                     (options.data_file, options.file_type))
        if options.file_type == "bam":
            # add BAM file
            bamfh = pysam.Samfile(options.data_file, "rb")
            intervalcoviter = BamCoverageIterator(bamfh,
                                                  norm_rlen=options.norm_rlen,
                                                  num_hits_tag=options.bam_nh,
                                                  hit_prob_tag=options.bam_prob,
                                                  max_multimaps=options.max_multihits,
                                                  keep_dup=options.keep_dup,
                                                  keep_qcfail=False,
                                                  flip_read2_strand=options.fr)
            t.fromintervals(intervalcoviter)
            # store coverage statistics to allow normalization calculations
            stats = intervalcoviter.stats
            logging.debug("\tProcessed '%d' valid reads" % (stats.num_reads))
            logging.debug("\tTotal coverage '%f'" % (stats.total_cov))
            bamfh.close()
    tf.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bam_file', help="input bam file")
    parser.add_argument('h5_file', help="output HDF5 file")
    args = parser.parse_args()

    refs = get_refs_from_bam(args.bam_file)
    gtracks = GenomeTrack(args.h5_file, refs)
    print(dict(gtracks.h5f.attrs['refs']))

    pysamfh = pysam.AlignmentFile(args.bam_file)
    gtracks.add_alignments(pysamfh)



def fromtophat(self, accepted_hits_bam, junctions_bed,
               max_multimaps=None,
               flip_read2_strand=True):
    # insert splice junction track
    track_name = self.hdf_group._v_name
    rnames = set(self.get_rnames())
    logging.info("[RnaseqTrack] adding junctions")
    junc_iter = tophat_bed_to_juncs(track_name, open(junctions_bed))
    for junc in junc_iter:
        if junc[REF_COL_NAME] not in rnames:
            logging.debug('Skipping junc %s' % str(junc))
        else:
            junc['id'] = self.junc_track.num_intervals
            self.junc_track.add(junc)
    self.junc_track.index(persist=True)
    # insert coverage track
    logging.info("creating coverage track")
    bamfh = pysam.Samfile(accepted_hits_bam, "rb")
    #cmdline = bamfh.header["PG"][0]["CL"]
    #re.search(r'--max-multihits(?:\s+|=)(\d+)', cmdline)
    intervalcoviter = BamCoverageIterator(bamfh,
                                          norm_rlen=True,
                                          num_hits_tag="NH",
                                          hit_prob_tag=None,
                                          max_multimaps=None,
                                          keep_dup=True,
                                          keep_qcfail=False,
                                          flip_read2_strand=True)
    self.cov_track.fromintervals(intervalcoviter)
    # store coverage statistics to allow normalization calculations
    stats = intervalcoviter.stats
    logging.debug("\tProcessed '%d' valid reads" % (stats.num_reads))
    logging.debug("\tTotal coverage '%f'" % (stats.total_cov))
    bamfh.close()


CIGAR_M = 0 #match  Alignment match (can be a sequence match or mismatch)
CIGAR_I = 1 #insertion  Insertion to the reference
CIGAR_D = 2 #deletion  Deletion from the reference
CIGAR_N = 3 #skip  Skipped region from the reference
CIGAR_S = 4 #softclip  Soft clip on the read (clipped sequence present in <seq>)
CIGAR_H = 5 #hardclip  Hard clip on the read (clipped sequence NOT present in <seq>)
CIGAR_P = 6 #padding  Padding (silent deletion from the padded reference sequence)

def get_genomic_intervals(read):
    intervals = []
    rseq = read.seq
    qseq = array.array('c')
    qstart = 0
    astart = read.pos
    aend = astart
    for op,length in read.cigar:
        if (op == CIGAR_D):
            aend += length
        elif (op == CIGAR_I) or (op == CIGAR_S):
            qstart += length
        elif (op == CIGAR_M):
            qseq.fromstring(rseq[qstart:qstart + length])
            qstart += length
            aend += length
        elif (op == CIGAR_N):
            if aend > astart:
                if len(qseq) != (aend - astart):
                    logging.error("Read %s has aend != astart" % (str(read)))
                else:
                    intervals.append((astart, aend, qseq))
            astart = aend + length
            aend = astart
            qseq = array.array('c')
    if aend > astart:
        if len(qseq) != (aend - astart):
            logging.error("Read %s has aend != astart" % (str(read)))
        else:
            intervals.append((astart, aend, qseq))
    if aend != read.aend:
        logging.error("Read %s has aend != read.aend" % (str(read)))
    return intervals



if __name__ == "__main__":
    # execute only if run as a script
    main()
