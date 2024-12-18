from heapq import merge


class zip_coordinates:
    """
    Merge two object that have coordinate chromosome/position.
    The format of the objects must be a tuple with (coordinates, data)
    where coordinate is a tuple with chromosome,position_start, position_end
    and data is a tuple with the data. The data of the two object will be
    merged for matching lines.
    For the first object only the start coordinate is taken into account.
    """

    def __init__(self, item1, item2):
        self.c2 = item2
        try:
            coordinates, self._last_data = next(self.c2)
        except StopIteration:
            coordinates, self._last_data = ((None, 0, 0), (None,))
        self._chromosome, self._last_window_s, self._last_window_e = coordinates
        self.c1 = item1
        self._last_chromosome = None

    _sentinel = object()

    def __next__(self):
        return self.next()

    def next(self):
        self.c1_line = next(self.c1)
        going_on = True
        while going_on:
            if self._chromosome == self.c1_line[0][0]:
                self._last_chromosome = self._chromosome
                if (
                    self.c1_line[0][1] >= self._last_window_s
                    and self.c1_line[0][1] < self._last_window_e
                ):
                    data = self.c1_line[1] + self._last_data
                    return (self.c1_line[0], data)
                    going_on = False
                elif self.c1_line[0][1] < self._last_window_s:
                    self.c1_line = next(self.c1)
                elif self.c1_line[0][1] >= self._last_window_e:
                    coordinates, self._last_data = next(self.c2)
                    (
                        self._chromosome,
                        self._last_window_s,
                        self._last_window_e,
                    ) = coordinates
            else:
                if (
                    self._last_chromosome != self._chromosome
                    and self._last_chromosome is not None
                ):
                    self.c1_line = next(self.c1)
                else:
                    coordinates, self._last_data = next(self.c2)
                    (
                        self._chromosome,
                        self._last_window_s,
                        self._last_window_e,
                    ) = coordinates

    def close(self):
        self.c1.close()
        self.c2.close()

    def __iter__(self):
        return iter(self.next, self._sentinel)


def chrompos_keyfunc(line, i):
    chrom = line[0]
    chrom = chrom.replace("chr", "")
    try:
        chrom = int(chrom)
    except ValueError:
        chrom = int(chrom, 36)

    return (chrom, line[1], i)


def decorated_item(f, i):
    for line in f:
        yield (chrompos_keyfunc(line[0], i), line)


def merge_items(item1, item2):
    files = [item1, item2]

    for line in merge(*[decorated_item(files[i], i + 1) for i in [0, 1]]):
        yield line[1]


def zip_fast(item1, item2):
    """
    Use the native implementation of the heapq algorithm to sort
    and merge files chromosome-coordinate ordered.
    It assumes that the two files are position ordered and both
    files have the same chromosome order.
    It differs from zip_coordinates by the fact that this return
    all the position present in both files, group together the lines
    present in both
    """

    merged = merge_items(item1, item2)

    store_line = next(merged)
    for line in merged:
        if store_line[0][0] == line[0][0] and store_line[0][2] > line[0][1]:
            yield (store_line[0], (store_line[1][0], line[1][0]))
            if store_line[0][2] == line[0][2] and store_line[0][1] == line[0][1]:
                store_line = next(merged)
        else:
            yield (store_line[0], (store_line[1][0], None))
            store_line = line
    if store_line == line:
        yield (store_line[0], (store_line[1][0], None))
