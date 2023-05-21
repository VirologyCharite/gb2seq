from functools import lru_cache

try:
    from importlib.resources import open_text
except ImportError:
    from importlib_resources import open_text

from dark.fasta import FastaReads


@lru_cache()
def getSequence(filename, id_=None):
    """
    Load a FASTA file.

    @param filename: The C{str} FASTA file name.
    @param id_: The C{str} id of the sequence wanted, or C{None} if to retrieve
        the first sequence.
    """
    with open_text("gb2seq.data", filename) as fp:
        for read in FastaReads(fp):
            if id_ is None or read.id.split()[0] == id_:
                return read

    raise ValueError(f"Sequence {id_} not found in {filename!r}.")
