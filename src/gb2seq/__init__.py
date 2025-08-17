from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("gb2seq")
except PackageNotFoundError:
    # Package is not installed.
    __version__ = "unknown"


class Gb2SeqError(Exception):
    "A gb2seq library error occurred."
