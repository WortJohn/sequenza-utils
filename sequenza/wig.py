class Wiggle:
    """
    Read/write wiggle files as iterable objects.
    """

    def __init__(self, wig):
        self._wig = wig
        self._chromosome = None
        self._span = None
        self._start = None
        self._end = None

    _sentinel = object()

    def __next__(self):
        return self.next()

    def next(self):
        try:
            wig_line = next(self._wig)
        except StopIteration:
            raise StopIteration
        try:
            pos, self._val = wig_line.strip().split("\t", 1)
            self._start = int(pos)
            self._end = self._start + self._span
            return ((self._chromosome, self._start, self._end), (self._val,))
        except ValueError:
            line, chrom, span = wig_line.strip().split()
            self._chromosome = chrom.split("chrom=")[1]
            self._span = int(span.split("span=")[1])
            pos, self._val = next(self._wig).strip().split("\t", 1)
            self._start = int(pos)
            self._end = self._start + self._span
            return ((self._chromosome, self._start, self._end), (self._val,))
        except StopIteration:
            raise StopIteration
        except ValueError:
            if self._chromosome:
                raise StopIteration
            else:
                raise WiggleError(
                    (
                        "Error: the file is not in the wiggle format: "
                        "https://genome.ucsc.edu/goldenpath/help/wiggle.html"
                    )
                )

    def write(self, item):
        try:
            try:
                chromosome = item[0]
            except AttributeError:
                chromosome = item[0]
            span = int(item[2]) - int(item[1])
            if span < 0:
                raise WiggleError("Negaive span value")
            if chromosome != self._chromosome or span != self._span:
                self._chromosome = chromosome
                self._span = span
                self._wig.write(
                    "".join(
                        map(
                            str,
                            [
                                "variableStep chrom=",
                                chromosome,
                                " span=",
                                int(span),
                                "\n",
                            ],
                        )
                    )
                )
                self._wig.write("\t".join(map(str, [item[1], item[3], "\n"])))
            else:
                self._wig.write("\t".join(map(str, [item[1], item[3], "\n"])))
        except WiggleError:
            raise

    def close(self):
        self._wig.close()

    def __iter__(self):
        return iter(self.next, self._sentinel)


class WiggleError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)
