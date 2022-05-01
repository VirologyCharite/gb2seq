from pathlib import Path


class Sars2SeqError(Exception):
    "A sars2seq library error occurred"


__version__ = "0.1.2"

DATA_DIR = Path(__file__).parent.parent / "data"
