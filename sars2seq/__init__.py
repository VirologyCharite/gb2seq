from pathlib import Path

__version__ = '0.0.3'


class Sars2SeqError(Exception):
    'A sars2seq library error'


DATA_DIR = Path(__file__).parent.parent / 'data'
