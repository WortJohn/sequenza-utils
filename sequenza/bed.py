


class SimpleBed:
    """
    Read BED files as iterable objects.
    """

    def __init__(self, bed, offset_left=0, offset_right=0):
        self._bed = bed
        self._chromosome = None
        self._start = None
        self._end = None
        self._offset_left = offset_left
        self._offset_right = offset_right

    _sentinel = object()

    def __next__(self):
        return self.next()

    def next(self):
        try:
            bed_line = next(self._bed)
        except StopIteration:
            raise StopIteration
        try:
            self._chromosome, self._start, self._end, self._val = bed_line.strip().split("\t")
            start = int(self._start) + self._offset_left
            end = int(self._end) + self._offset_left
            return ((self._chromosome, start, end), (self._val,))
        except StopIteration:
            raise StopIteration
        except ValueError:
            raise StopIteration

    def close(self):
        self._wig.close()

    def __iter__(self):
        return iter(self.next, self._sentinel)



def seqz_from_beds(beds):
    seqz_template = "%s\t%i\tN\t%i\t%i\t%s\t1.0\t0\thom\t%s\t%s\tN\t.\t0\n"
    for line in beds:
        depth_0 = float(line[1][0])
        depth_1 = float(line[1][1])
        if depth_0 == 0:
            continue
        try:
            gc = int(line[1][2])
        except IndexError:
            gc = 50
        yield seqz_template % (
            line[0][0],
            int(line[0][2]),
            int(depth_0),
            int(depth_1),
            round(depth_1/float(depth_0), 3),
            gc, depth_1
        )
