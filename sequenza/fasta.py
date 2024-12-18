class Fasta:
    """
    Creates an iterable with genomic coordinates from a fasta file
    """

    def __init__(self, file, n=60):
        self.fasta = file
        self.chromosome = None
        self.start = 1
        self.n = int(n)
        self.buffer = ""

    _sentinel = object()

    def next(self):
        if len(self.buffer) > self.n:
            fasta_line = ""
        else:
            if len(self.buffer) > 0:
                try:
                    fasta_line = next(self.fasta).strip()
                except StopIteration:
                    results = (
                        self.chromosome,
                        self.start,
                        self.start + len(self.buffer) - 1,
                        self.buffer,
                    )
                    self.buffer = ""
                    return results
            else:
                fasta_line = next(self.fasta).strip()
        if fasta_line.startswith(">"):
            if len(self.buffer) > 0:
                results = (
                    self.chromosome,
                    self.start,
                    self.start + len(self.buffer) - 1,
                    self.buffer,
                )
            else:
                results = None
            self.chromosome = fasta_line.lstrip(">").split()[0]
            self.buffer = ""
            self.start = 1
        else:
            seq = self.buffer + fasta_line.upper()
            if len(seq) >= self.n:
                results = (
                    self.chromosome,
                    self.start,
                    self.start + self.n - 1,
                    seq[0 : self.n],
                )
                self.buffer = seq[self.n : len(seq)]
                self.start = self.start + self.n
            else:
                self.buffer = seq
                results = None
        return results

    def close(self):
        self.fasta.close()

    def __iter__(self):
        return iter(self.next, self._sentinel)
