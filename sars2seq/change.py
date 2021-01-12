def splitChange(s):
    """
    Split a change string.

    @param s: The change C{str} to split. This has the format RNG, where R is
        a reference base, N is a 1-based location, and G is a genome sequence
        base. So, e.g., 'L28S P1003Q' indicates that we expected a change from
        'L' to 'S' at offset 28 and from 'P' in the reference to 'Q' in the
        genome we're examining at offset 1003. The reference or genome base
        (but not both) may be absent.
    @return: A 3-C{tuple} with reference base or C{None}, 0-based C{int}
        offset, genome base or C{None}.
    """
    try:
        location = int(s)
    except ValueError:
        pass
    else:
        raise ValueError(f'Change string {s!r} does not include a '
                         f'reference or genome base.')

    try:
        location = int(s[1:])  # E.g., D400
    except ValueError:
        pass
    else:
        return s[0], location - 1, None

    # Don't let a leading gap indicator be taken as negative location.
    if s[0] != '-':
        try:
            location = int(s[:-1])  # E.g., 400D
        except ValueError:
            pass
        else:
            return None, location - 1, s[-1]

    try:
        location = int(s[1:-1])  # E.g., C400D
    except ValueError:
        pass
    else:
        return s[0], location - 1, s[-1]

    raise ValueError(f'Could not parse change string {s!r}.')
