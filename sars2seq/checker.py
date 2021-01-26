import copy

class Checker:
    def __init__(self, name, changes, nt):
        """
        Check genome changes occurred.

        @param name: The C{str} name of the feature to check (e.g., 'nsp2').
        @param changes: A C{str} specification in the form of space-separated
            RNS strings, where R is a reference base, N is an integer offset,
            and S is a sequence base. So, e.g., 'L28S P1003Q' indicates that
            we expected a change from 'L' to 'S' at offset 28 and from 'P' to
            'Q' at offset 1003.
        @param nt: If C{True} check nucleotide sequences. Else protein.
        """
        self.description = (f'Check {name!r} for {"nt" if nt else "aa"} '
                            f'changes {changes!r}')

        def check(genome):
            _, errorCount, _ = genome.feature(name).check(changes, nt)
            return errorCount == 0

        self._func = check

    def __str__(self):
        return self.description

    def __call__(self, genome):
        return self._func(genome)

    def __and__(self, other):

        def check(genome):
            return self(genome) and other(genome)

        newChecker = copy.copy(self)
        newChecker._func = check
        newChecker.description = f'({self.description} AND {other.description})'
        return newChecker

    def __or__(self, other):

        def check(genome):
            return self(genome) or other(genome)

        newChecker = copy.copy(self)
        newChecker._func = check
        newChecker.description = f'({self.description} OR {other.description})'
        return newChecker


class AAChecker(Checker):
    def __init__(self, name, changes):
        super().__init__(name, changes, False)


class NTChecker(Checker):
    def __init__(self, name, changes):
        super().__init__(name, changes, True)
