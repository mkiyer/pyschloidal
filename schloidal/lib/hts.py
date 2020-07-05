import collections
import pysam
import numpy as np
from enum import IntEnum

class Config:
    MAXSIZE = 1000000
    DTYPE = np.uint32


class Strand:
    POS = 0
    NEG = 1
    NONE = 2
    _STR_TO_INT = {'+': POS, '-': NEG, '.': NONE}
    _INT_TO_STR = ['+', '-', '.']

    @staticmethod
    def from_str(s):
        return Strand._STR_TO_INT[s]

    @staticmethod
    def to_str(s):
        return Strand._INT_TO_STR[s]


class Cigar(IntEnum):
    """
    CIGAR operations
    CIGAR: M = match or mismatch
    CIGAR: I = insertion to the reference
    CIGAR: D = deletion from the reference
    CIGAR: N = skip on the reference (e.g. spliced alignment)
    CIGAR: S = clip on the read with clipped sequence
    CIGAR: H = clip on the read with clipped sequence trimmed off
    CIGAR: P = padding
    CIGAR: equals = match
    CIGAR: X = mismatch
    """
    MATCH = 0
    INS = 1
    DEL = 2
    REF_SKIP = 3
    SOFT_CLIP = 4
    HARD_CLIP = 5
    PAD = 6
    EQUAL = 7
    DIFF = 8


def get_refs_from_hts_file(htsfile):
    f = pysam.AlignmentFile(htsfile)
    refs = list(zip(f.references, map(int, f.lengths)))
    f.close()
    return refs


def get_exons_and_introns(a):
    start = a.reference_start
    pos = start
    exons = []
    introns = []
    for op, nbases in a.cigar:
        if op == Cigar.MATCH or op == Cigar.EQUAL or op == Cigar.DIFF or op == Cigar.DEL:
            pos += nbases
        elif op == Cigar.REF_SKIP:
            if pos > start:
                exons.append((start, pos))
                introns.append((pos, pos + nbases))
            start = pos + nbases
            pos = start
    if pos > start:
        exons.append((start, pos))
    return exons, introns


class PileupBuffer(object):
    def __init__(self, flush_callback_func, maxsize, dtype):
        self.flush_callback_func = flush_callback_func
        self.maxsize = maxsize
        self.ref = None
        self.start = 0
        self.end = 0
        self.cov = np.zeros((maxsize, 3), dtype=dtype)
        self.cov_spliced = np.zeros((maxsize, 2), dtype=dtype)
        self.recur = np.zeros((maxsize, 3), dtype=dtype)
        self.recur_spliced = np.zeros((maxsize, 2), dtype=dtype)

    def clear(self, ref=None, start=0, end=None):
        self.ref = ref
        self.start = start
        end = start if end is None else end
        self.end = end
        self.cov.fill(0)
        self.cov_spliced.fill(0)
        self.recur.fill(0)
        self.recur_spliced.fill(0)

    def flush(self, pos):
        if pos > self.start:
            self.flush_callback_func(self, self.start, pos)

    def _write(self, ref, start, end, strand, is_spliced, value):
        astart = start - self.start
        aend = end - self.start
        self.cov[astart:aend, strand] += value
        if is_spliced:
            self.cov_spliced[astart:aend, strand] += value

    def add(self, ref, strand, intervals, value=1):
        tstart, tend = intervals[0][0], intervals[-1][-1]
        tsize = tend - tstart
        is_spliced = len(intervals) > 1

        # reset buffer when reference changes
        if ref != self.ref:
            self.flush(self.ref, self.start, self.end)
            self.clear(ref, tstart)

        # flush buffer when full
        if tend > self.start + self.maxsize:
            self.flush(self.ref, self.start, self.end)
            self.clear(ref, tstart)




        # reset buffer when empty or when reference changes
        if ref != self.ref:
            self.flush(self.ref, self.start, self.end)
            self.clear(ref, tstart)

        for start, end in intervals:
            # check if need to flush buffer
            if end > self.maxend:
                self.flush(ref, start)
                # need to shift the data to the left

            while end > self.maxend:
                # store part of this interval before flushing
                if start < self.maxend:
                    self.end = self.maxend
                    self._write(ref, start, self.maxend, strand, is_spliced, value)
                    start = self.maxend
                # reset buffer
                self.flush(ref, start)
                self.clear(ref, start)
            # keep track of part of buffer that is being used            
            if end > self.end:
                self.end = end
            # write the rest of the interval
            self._write(ref, start, end, strand, is_spliced, value)


def parse_alignments(aligniter):

    def flusher(buf):
        print("flushing", buf.ref, buf.start, buf.end)

    readbuf = PileupBuffer(flusher, Config.MAXSIZE, Config.DTYPE)
    introns = collections.OrderedDict()

    
    for a in aligniter:
        if a.is_unmapped or a.is_qcfail or a.is_secondary:
            continue
        ref = a.reference_id
        strand = a.get_tag('XS')
        strand = Strand.from_str(strand)
        exons, introns = get_exons_and_introns(a)
        readbuf.add(ref, strand, exons)
    
    readbuf.flush()

