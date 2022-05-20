from typing import Optional, Tuple


def splitChange(s: str) -> Tuple[Optional[str], int, Optional[str]]:
    """
    Split a change string.

    @param s: The change C{str} to split. This has the format RNG, where R is
        a reference base, N is a 1-based site (i.e., genome location), and G
        is a genome sequence base. So, e.g., 'L28S P1003Q' indicates that we
        expected a change from 'L' to 'S' at site 28 and from 'P' in the
        reference to 'Q' in the genome we're examining at site 1003. The
        reference or genome base (but not both) may be absent.
    @return: A 3-C{tuple} with a C{str} reference base or C{None}, 0-based
        C{int} offset, and a C{str} genome base or C{None}.
    """
    if not s:
        # Passing any false value (e.g., 0 or None), not just the empty string,
        # will trigger this exception.
        raise ValueError(f"Could not parse change argument {s!r}.")

    try:
        site = int(s)
    except ValueError:
        pass
    else:
        raise ValueError(
            f"Change string {s!r} does not include a " f"reference or genome base."
        )

    try:
        site = int(s[1:])  # E.g., D400
    except ValueError:
        pass
    else:
        return s[0], site - 1, None

    # Don't let a leading gap indicator be taken as negative site.
    if s[0] != "-":
        try:
            site = int(s[:-1])  # E.g., 400D
        except ValueError:
            pass
        else:
            return None, site - 1, s[-1]

    try:
        site = int(s[1:-1])  # E.g., C400D
    except ValueError:
        pass
    else:
        return s[0], site - 1, s[-1]

    raise ValueError(f"Could not parse change string {s!r}.")
