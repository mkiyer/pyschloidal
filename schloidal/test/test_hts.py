import numpy as np

from schloidal.lib.hts import PileupBuffer, Strand


def test_pileup():
    def flush(b):
        print(b.ref, b.start, b.end)
    buf = PileupBuffer(flush, 5, np.uint32)
    # add intervals
    buf.add(0, Strand.POS, [(0, 1), (3, 4)])
    buf.add(0, Strand.POS, [(0, 5)])
    assert np.array_equal(buf.cov[0:5, Strand.POS], np.array([2, 1, 1, 2, 1]))
    assert np.array_equal(buf.cov_spliced[0:5, Strand.POS], np.array([1, 0, 0, 1, 0]))


def test_flush():
    x = [(3, 13), (13, 23)]
    y = []
    a = np.zeros((29,3), np.uint32)

    def flush(b):
        y.append((b.start, b.end))
        a[b.start:b.end, :] = b.cov
        print(b.ref, b.start, b.end, b.cov)

    buf = PileupBuffer(flush, 10, np.uint32)
    buf.add(0, Strand.POS, [(3, 28)])
    buf.add(0, Strand.POS, [(4, 29)])
    print(a)

    #assert all(x[i] == y[i] for i in range(len(x)))
    #assert np.array_equal()

    #assert np.array_equal(buf.cov[0:5, Strand.POS], np.array([2, 1, 1, 2, 1]))
    #assert np.array_equal(buf.cov_spliced[0:5, Strand.POS], np.array([1, 0, 0, 1, 0]))
