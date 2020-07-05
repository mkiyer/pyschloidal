"""
schloidal
"""
import os
import sys
import re
import argparse
import collections

import pysam
import numpy as np
import h5py

from schloidal.lib.hts import parse_alignments, get_refs_from_hts_file

STRAND_POS = 0
STRAND_NEG = 1
STRAND_NONE = 2
STRAND_STR_TO_INT = {'+': 0, '-': 1, '.': 2}
STRAND_INT_TO_STR = ['+', '-', '.']

BEDGRAPH_FORMAT_STRING = "{}\t{}\t{}\t{}"

_interval_string_re = re.compile(r'^(?P<ref>\w+)(?:\[(?P<strand>[+-.])\])?(?::(?P<start>[\d,]+)(?:-(?P<end>[\d,]+))?)?')


def _parse_interval_string(interval_string):
    """ref[strand]:start-end

    `strand` can be '+', '-', or '.'
    """
    m = _interval_string_re.match(interval_string)
    if m is None:
        ref, start, end, strand = None, None, None, STRAND_NONE
        # if interval string is just a single character,
        # try to extract something useful
        if len(interval_string) == 1:
            strand = STRAND_STR_TO_INT[interval_string[0]]
    else:
        ref = m.group('ref')
        start = m.group('start')
        end = m.group('end')
        strand = m.group('strand')
        if start is not None:
            start = int(start.replace(",", ""))
        if end is not None:
            end = int(end.replace(",", ""))
        if strand is None:
            strand = STRAND_NONE
        else:
            strand = STRAND_STR_TO_INT[strand]
    return ref, start, end, strand


def _parse_interval_tuple(interval):
    if interval is None:
        return None, None, None, STRAND_NONE
    if len(interval) == 1:
        return interval[0], None, None, STRAND_NONE
    elif len(interval) == 2:
        return interval[0], interval[1], interval[1] + 1, STRAND_NONE
    elif len(interval) == 3:
        return interval[0], interval[1], interval[2], STRAND_NONE
    else:
        if isinstance(interval[3], str):
            strand = STRAND_STR_TO_INT(interval[3])
        else:
            strand = interval[3]
        return interval[0], interval[1], interval[2], strand


def parse_interval(interval):
    """parses a genomic interval specifier in either the string
    format `<ref[strand]:start-end>` or tuple (ref, start, end, strand)

    :param interval: interval specifier
    :type interval: str or tuple

    >>> parse_interval("chr1:100-200")
    ("chr1", 100, 200, ".")
    >>> parse_interval("chr1:1,000-20,000")
    ("chr1", 1000, 20000)
    >>> parse_interval("chr1[-]:1,000-20,000")
    ("chr1", 1000, 20000, "-")
    >>> parse_interval("chrX")
    ("chrX", None, None, ".")
    >>> parse_interval("chrX:5")
    ("chrX", 5, 6, ".")
    >>> parse_interval(("chr2", 100, 300))
    ("chr2", 100, 300, ".")
    """
    if isinstance(interval, str):
        interval = _parse_interval_string(interval)
    interval = _parse_interval_tuple(interval)
    return interval





def iter_pe_alignments(aligniter):
    """
    generator function to return a tuple of lists of paired-end reads
    iterator must be sorted by query name

    :param aligniter: pysam AlignedSegment iterator
    :return:
    """
    pe_alignments = ([], [])
    num_reads = 0
    prev_qname = None
    for a in aligniter:
        if not a.is_proper_pair:
            continue
        if a.is_unmapped or a.is_qcfail or a.is_secondary:
            continue
        # get read attributes
        qname = a.query_name
        readnum = 1 if a.is_read2 else 0
        # if query name changes we have finishes with this fragment
        if num_reads > 0 and qname != prev_qname:
            yield pe_alignments
            pe_alignments = ([], [])
            num_reads = 0
        pe_alignments[readnum].append(a)
        prev_qname = qname
        num_reads += 1
    if num_reads > 0:
        yield pe_alignments


def get_concordant_mate_pairs(pe_alignments):
    """
    :param pe_alignments: tuple of lists ([read1 alignments], [read2 alignments])
    :return: list of mate pairs [[read1, read2], [read1, read2]]
    """
    if any((len(a) == 0) for a in pe_alignments):
        return []
    # index read1 by mate reference name and position
    mate_dict = collections.defaultdict(lambda: collections.deque())
    for a1 in pe_alignments[0]:
        mate_dict[(a1.next_reference_id, a1.next_reference_start)].append(a1)
    # iterate through read2 and get mate pairs
    mate_pairs = []
    for a2 in pe_alignments[1]:
        a1 = mate_dict[(a2.reference_id, a2.reference_start)].popleft()
        mate_pairs.append((a1, a2))
    return mate_pairs


def _array_to_bedgraph(fileh, ref, start, arr, factor):
    end = start
    val = 0
    for i in range(0, arr.shape[0]):
        n = arr[i]
        if val != n:
            if (start != end) and (val != 0):
                print(BEDGRAPH_FORMAT_STRING.format(ref, start, end, factor * val), file=fileh)
            start = end
            val = n
        end += 1
    if start != end:
        print(BEDGRAPH_FORMAT_STRING.format(ref, start, end, factor * val), file=fileh)


def array_to_bedgraph(fileh, ref, start, end, strand, arr, factor, chunk_size):
    if end < start:
        end = arr.shape[1]
    if chunk_size > (end - start):
        chunk_size = (end - start)
    if strand == STRAND_NONE:
        strands = (STRAND_POS, STRAND_NEG)
    else:
        strands = (strand,)
    while start < (end - chunk_size):
        # extract chunk from array
        chunk_arr = arr[strands, start:start + chunk_size].sum(axis=0)
        # write chunk
        _array_to_bedgraph(fileh, ref, start, chunk_arr, factor)
        start += chunk_size
    if start < end:
        chunk_arr = arr[strands, start:end].sum(axis=0)
        _array_to_bedgraph(fileh, ref, start, chunk_arr, factor)


class GenomeTracks:
    # constants to control chunking behavior
    H5_RDCC_NBYTES = 1 << 24
    H5_RDCC_W0 = 0
    H5_RDCC_NSLOTS = 10007
    H5_FILE_KWARGS = {'rdcc_nbytes': H5_RDCC_NBYTES,
                      'rdcc_w0': H5_RDCC_W0,
                      'rdcc_nslots': H5_RDCC_NSLOTS}

    DTYPE = np.uint64
    COMPRESSION = 'lzf'
    CHUNK_SIZE = 1<<20

    # attribute names
    REFS = 'refs'
    NUM_FRAGS = 'num_frags'

    def _parse_interval(self, key):
        ref, start, end, strand = parse_interval(key)
        if ref is None:
            return None, None, None, strand
        assert ref in self._refdict
        if start is None:
            start = 0
        if end is None:
            end = int(self._refdict[ref])
        return ref, start, end, strand

    def _create(self, filename, refs, mode):
        # create new file
        self.h5f = h5py.File(filename, mode, **GenomeTracks.H5_FILE_KWARGS)
        # set attributes
        self.h5f.attrs[GenomeTracks.REFS] = refs
        self.h5f.attrs[GenomeTracks.NUM_FRAGS] = 0
        # create datasets
        for ref, length in refs:
            chunk_size = min(length, GenomeTracks.CHUNK_SIZE)
            self.h5f.create_dataset(ref, (2, length,),
                                    dtype=GenomeTracks.DTYPE,
                                    chunks=(1, chunk_size),
                                    compression=GenomeTracks.COMPRESSION,
                                    shuffle=True,
                                    fillvalue=0)

    def __init__(self, filename, refs=None, mode='a'):
        if h5py.is_hdf5(filename):
            # open existing file
            self.h5f = h5py.File(filename, mode, **GenomeTracks.H5_FILE_KWARGS)
        else:
            # create new file
            self._create(filename, refs, mode)
        # store dictionary of refs for easy lookup
        self._refdict = collections.OrderedDict(self.h5f.attrs[GenomeTracks.REFS])

    def close(self):
        self.h5f.flush()
        self.h5f.close()

    def get_references(self):
        for ref, length in self.h5f.attrs[GenomeTracks.REFS]:
            print(f"{ref}\t{length}")

    def add_alignments(self, pysamiter):
        num_frags = 0
        for pe_alignments in iter_pe_alignments(pysamiter):
            pairs_list = get_concordant_mate_pairs(pe_alignments)
            assert len(pairs_list) == 1
            a1, a2 = pairs_list[0]
            reference_name = a1.reference_name
            strand = a1.get_tag('XS')
            assert strand == a2.get_tag('XS')
            strand = STRAND_STR_TO_INT[strand]
            if a1.reference_start < a2.reference_start:
                start = a1.reference_start
                end = a2.reference_end
            else:
                start = a2.reference_start
                end = a1.reference_end
            assert start < end
            assert reference_name == a2.reference_name
            dset = self.h5f[reference_name]
            dset[np.s_[strand, start:end]] += 1
            num_frags += 1
            if num_frags % 1000 == 0:
                print('num_frags', num_frags)
        print('num_frags', num_frags)
        self.h5f.attrs[GenomeTracks.NUM_FRAGS] += num_frags

    def tobedgraph(self, fileh, interval, norm=True, multiplier=1.0e6):
        ref, start, end, strand = self._parse_interval(interval)
        if ref is None:
            refs = list(self._refdict.keys())
        else:
            refs = [ref]
        if start is None: start = 0
        if end is None: end = -1
        if norm:
            num_frags = self.h5f.attrs[GenomeTracks.NUM_FRAGS]
            factor = 1.0 * multiplier / num_frags
        else:
            factor = 1.0
        for ref in refs:
            array_to_bedgraph(fileh, ref, start, end, strand, self.h5f[ref], factor, self.CHUNK_SIZE)


def main_new(parser, args):
    h5_file = args.h5_file
    bam_file = args.bam_file
    if os.path.exists(h5_file):
        parser.error(f"HDF5 file {h5_file} already exists")
    refs = get_refs_from_hts_file(bam_file)
    gtracks = GenomeTracks(args.h5_file, refs)
    pysamfh = pysam.AlignmentFile(bam_file)
    pysamfh.close()
    gtracks.close()


def main_add(parser, args):
    gtracks = GenomeTracks(args.h5_file)
    pysamfh = pysam.AlignmentFile(args.bam_file)
    parse_alignments(pysamfh)
    #gtracks.add_alignments(pysamfh)
    pysamfh.close()
    gtracks.close()


def main_view(parser, args):
    h5_file = args.h5_file
    region = parse_interval(args.region)
    gtracks = GenomeTracks(h5_file)
    gtracks.tobedgraph(sys.stdout, region)
    gtracks.close()


def main_refs(parser, args):
    gtracks = GenomeTracks(args.h5_file)
    gtracks.get_references()
    gtracks.close()


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="commands",
                                       description="valid subcommands",
                                       help="")
    # create the parser for the "new" command
    parser_new = subparsers.add_parser("new", help="create a new genome tracks file")
    parser_new.add_argument("h5_file", help="genome tracks HDF5 file to create")
    parser_new.add_argument('bam_file', help="BAM file")
    parser_new.set_defaults(func=main_new)

    # create the parser for the "add" command
    parser_add = subparsers.add_parser("add", help="add a BAM file")
    parser_add.add_argument("h5_file", help="genome tracks HDF5 file")
    parser_add.add_argument('bam_file', help="BAM file to add")
    parser_add.set_defaults(func=main_add)

    # create parser for "view" command
    parser_view = subparsers.add_parser("view", help="view data in a track")
    parser_view.add_argument('h5_file', help='genome tracks HDF5 file')
    parser_view.add_argument('region', nargs="?", default=None,
                             help="genomic region")
    parser_view.set_defaults(func=main_view)

    # create parser for "refs" command
    parser_view = subparsers.add_parser("refs", help="get references (chrom.sizes)")
    parser_view.add_argument('h5_file', help='genome tracks HDF5 file')
    parser_view.set_defaults(func=main_refs)

    # parse args
    args = parser.parse_args()
    args.func(parser, args)


if __name__ == "__main__":
    # execute only if run as a script
    main()
