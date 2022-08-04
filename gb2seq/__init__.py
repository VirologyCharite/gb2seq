from pathlib import Path


class Gb2SeqError(Exception):
    "A gb2seq library error occurred"


__version__ = "0.2.0"

DATA_DIR = Path(__file__).parent.parent / "data"
