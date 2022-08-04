from warnings import warn

from .alignment import (
    SARS2Alignment as SARS2Genome,
    addAlignerOption,
    getGappedOffsets,
    alignmentEnd,
    offsetInfoMultipleGenomes,
    AlignmentError,
)

# Use the above so flake8 doesn't complain.
_ = (
    SARS2Genome,
    addAlignerOption,
    getGappedOffsets,
    alignmentEnd,
    offsetInfoMultipleGenomes,
    AlignmentError,
)

warn(
    "The SARS2Genome class has been renamed to SARS2Alignment. "
    "Please use that instead. Support for SARS2Genome will be removed "
    "in a future version."
)
