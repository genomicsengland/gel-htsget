import math
import struct
import gzip


def deserialize(fmt, f):
    return struct.unpack(fmt, f.read(struct.calcsize(fmt)))[0]


def reg2bin(beg, end):
    """
    Beg inclusive.
    End exclusive.
    """
    if beg >> 14 == (end - 1) >> 14: return int(
        ((1 << 15) - 1) / 7 + (beg >> 14))
    if beg >> 17 == (end - 1) >> 17: return int(
        ((1 << 12) - 1) / 7 + (beg >> 17))
    if beg >> 20 == (end - 1) >> 20: return int(
        ((1 << 9) - 1) / 7 + (beg >> 20))
    if beg >> 23 == (end - 1) >> 23: return int(
        ((1 << 6) - 1) / 7 + (beg >> 23))
    if beg >> 26 == (end - 1) >> 26: return int(
        ((1 << 3) - 1) / 7 + (beg >> 26))

    return 0


def reg2bins(beg, end):
    """
    Generates bin ids which overlap the specified region.

    Beg inclusive.
    End exclusive.
    """
    # Based off the algorithm presented in:
    # https://samtools.github.io/hts-specs/SAMv1.pdf

    # Bin calculation constants.
    BIN_ID_STARTS = [0, 1, 9, 73, 585, 4681]
    BIN_ID_SHIFTS = [29, 26, 23, 20, 17, 14]

    # Maximum range supported by specifications.
    MAX_RNG = (2 ** 29) - 1

    for start, shift in zip(BIN_ID_STARTS, BIN_ID_SHIFTS):
        i = beg >> shift if beg > 0 else 0
        j = end >> shift if end < MAX_RNG else MAX_RNG >> shift

        for bin_id_offset in range(i, j + 1):
            yield start + bin_id_offset


class Index(object):
    """
    Logical representation of a BAM Index.
    """

    def __init__(self):
        self.references = list()
        self.unmapped = None

    @staticmethod
    def file_consume_magic_number(f):
        if f.read(4) != b'BAI\x01':
            raise ValueError('BAI magic number not found')

    @staticmethod
    def file_produce_magic_number(f):
        f.write(b'BAI\x01')

    @staticmethod
    def file_consume_reference_num(f):
        return deserialize('<l', f)

    @staticmethod
    def file_consume_unmapped(f):
        return deserialize('<Q', f)

    def header_length(self):
        """
        Returns the byte length of the header of the corresponding bam file.
        """
        # Header is everything up to first alignment record.
        return self.references[0].intervals[0] >> 16

    def region_offset(self, ref, beg, end):
        """
        Returns the byte offset of the first chunk containing alignments within
        the specified region, and the offset within the chunk of the first
        record within the specified region.
        """
        return self.references[ref].region_offset(beg, end)

    def region_offset_iter(self, ref, beg, end):
        """
        Returns the byte offset of the first chunk containing alignments within
        the specified region, and the offset within the chunk of the first
        record within the specified region.
        """
        return self.references[ref].region_offset_iter(beg, end)

    @classmethod
    def from_file(cls, filename=None, mode='rb', fileobj=None):
        """
        Generates a Index instance from a file object.
        The file should be opened in binary read mode.
        """

        assert any([
            filename is not None,
            fileobj is not None,
        ]), 'Must specify filename or fileobj.'

        assert None in [
            filename,
            fileobj,
        ], 'Filename and fileobj are mutually exclusive.'

        ifs = fileobj if fileobj is not None else open(filename, mode)

        # Check that this file is tagged as a BAI file.
        cls.file_consume_magic_number(ifs)

        self = cls()

        # Generate reference objects contained by the BAI.
        for i in range(cls.file_consume_reference_num(ifs)):
            self.references.append(Reference.from_file(fileobj=ifs))

        self.unmapped = cls.file_consume_unmapped(ifs)

        return self

    def to_dict(self):
        """
        Returns a dict representation of the BAM Index.
        """
        return {
            'references': list(map(lambda x: x.to_dict(), self.references)),
            'unmapped': self.unmapped,
        }

    @classmethod
    def from_dict(cls, d):
        """
        Generates a Index instance from a dict object.
        """
        self = cls()

        dict_to_reference = lambda x: Reference.from_dict(x)
        self.references = list(map(dict_to_reference, d['references']))

        self.unmapped = d['unmapped']

        return self


class Reference(object):
    """
    Logical representation of a BAM Index Reference.
    """

    # Smallest bp coverage increment.
    TILE_SIZE = 2 ** 14

    # Maximum range supported by specifications.
    MAX_RNG = (2 ** 29) - 1

    # Bin calculation constants.
    BIN_ID_STARTS = [0, 1, 9, 73, 585, 4681]
    BIN_ID_SHIFTS = [29, 26, 23, 20, 17, 14]

    # Metadata bin id.
    METADATA_ID = 37450

    def __init__(self):
        self.bins = dict()
        self.intervals = list()

        # Optional metadata
        self.unmapped_beg = None
        self.unmapped_end = None
        self.num_mapped = None
        self.num_unmapped = None

    @staticmethod
    def file_consume_num_bins(f):
        return deserialize('<i', f)

    @staticmethod
    def file_consume_num_intervals(f):
        return deserialize('<l', f)

    @staticmethod
    def file_consume_offset(f):
        return deserialize('<Q', f)

    @classmethod
    def region_to_bins(cls, beg, end):
        """
        Generates bin ids which overlap the specified region.
        """
        # Based off the algorithm presented in:
        # https://samtools.github.io/hts-specs/SAMv1.pdf
        for start, shift in zip(cls.BIN_ID_STARTS, cls.BIN_ID_SHIFTS):
            # TODO fix up handling of end
            end = end if end is not None else cls.MAX_RNG
            i = beg >> shift if beg > 0 else 0
            j = end >> shift if end <= cls.MAX_RNG else cls.MAX_RNG >> shift

            for bin_id_offset in range(i, j + 1):
                yield start + bin_id_offset

    def region_offset(self, beg, end):
        """
        Returns the byte offset of the first chunk containing alignments within
        the specified region, and the offset within the chunk of the first
        record within the specified region.
        """
        # Calculate the linear offset.
        linear = self.linear_offset(beg)
        offset = linear

        # Determine the bins overlapping the region.
        bins = type(self).region_to_bins(beg, end)
        for b in bins:

            # Skip any bins that weren't indexed - nothing in them.
            if b not in self.bins:
                continue

            # Find the earliest chunk that overlaps the linear offset.
            for lower, upper in self.bins[b].chunks:
                if upper < linear:
                    continue
                if lower < offset:
                    offset = lower

        # Convert and return the non-virtual offsets.
        return offset >> 16, offset & 0x0000FFFF

    def region_offset_iter(self, beg, end):
        """
        Returns the byte offset of the first chunk containing alignments within
        the specified region, and the offset within the chunk of the first
        record within the specified region.
        """
        # Calculate the linear offset.
        linear = self.linear_offset(beg)

        # Determine the bins overlapping the region.
        offsets = {}
        bins = type(self).region_to_bins(beg, end)
        for b in bins:

            # Skip any bins that weren't indexed - nothing in them.
            if b not in self.bins:
                continue

            # Find the earliest chunk that overlaps the linear offset.
            for lower, upper in self.bins[b].chunks:
                if upper < linear:
                    continue
                b = lower >> 16
                offsets[b] = lower

        for k in sorted(offsets.keys()):
            yield offsets[k] >> 16, offsets[k] & 0x0000FFFF

    def linear_offset(self, beg):
        """
        Returns the lower limit (virtual) file offset.
        """
        # Divides the genome up into TILE_SIZE bp chunks.
        # Calculates the tile a bp corresponds to, looks up the offset.
        try:
            return self.intervals[int(math.floor(beg / type(self).TILE_SIZE))]
        except IndexError:
            return -1

    @classmethod
    def from_file(cls, filename=None, mode='rb', fileobj=None):
        """
        Generates a Reference instance from a file object.
        The file should be opened in binary read mode.
        This function will consume the Reference subtree.
        """

        assert any([
            filename is not None,
            fileobj is not None,
        ]), 'Must specify filename or fileobj.'

        assert None in [
            filename,
            fileobj,
        ], 'Filename and fileobj are mutually exclusive.'

        ifs = fileobj if fileobj is not None else open(filename, mode)

        self = cls()

        for i in range(cls.file_consume_num_bins(ifs)):
            b = Bin.from_file(fileobj=ifs)

            # Process metadata bins separately.
            if b.bin_id == cls.METADATA_ID:
                self.unmapped_beg = b.chunks[0][0]
                self.unmapped_end = b.chunks[0][1]
                self.num_mapped = b.chunks[1][0]
                self.num_unmapped = b.chunks[1][1]
                continue

            self.bins[b.bin_id] = b

        for i in range(cls.file_consume_num_intervals(ifs)):
            self.intervals.append(cls.file_consume_offset(ifs))

        return self

    def to_dict(self):
        """
        Returns a dict representation of the BMA Index Reference.
        """
        bins = dict()
        for b in self.bins:
            bins[b] = self.bins[b].to_dict()

        return {
            'bins': bins,
            'intervals': self.intervals,
            'unmapped_beg': self.unmapped_beg,
            'unmapped_end': self.unmapped_end,
            'num_mapped': self.num_mapped,
            'num_unmapped': self.num_unmapped,
        }

    @classmethod
    def from_dict(cls, d):
        """
        Generates a Reference instance from a dict object.
        """
        self = cls()

        self.bins = dict()
        for b in d['bins']:
            self.bins[b] = Bin.from_dict(d['bins'][b])

        self.intervals = d['intervals']
        self.unmapped_beg = d['unmapped_beg']
        self.unmapped_end = d['unmapped_end']
        self.num_mapped = d['num_mapped']
        self.num_unmapped = d['num_unmapped']

        return self


class Bin(object):
    """
    Logical representation of a BAM Index Bin.
    """

    def __init__(self):
        self.bin_id = None
        self.chunks = list()

    @staticmethod
    def file_consume_bin_id(f):
        return deserialize('<I', f)

    @staticmethod
    def file_consume_num_chunks(f):
        return deserialize('<i', f)

    @staticmethod
    def file_consume_beg_off(f):
        return deserialize('<Q', f)

    @staticmethod
    def file_consume_end_off(f):
        return deserialize('<Q', f)

    @classmethod
    def from_file(cls, filename=None, mode='rb', fileobj=None):
        """
        Generates a Bin instance from a file object.
        The file should be opened in binary read mode.
        This function will consume the Bin subtree.
        """

        assert any([
            filename is not None,
            fileobj is not None,
        ]), 'Must specify filename or fileobj.'

        assert None in [
            filename,
            fileobj,
        ], 'Filename and fileobj are mutually exclusive.'

        ifs = fileobj if fileobj is not None else open(filename, mode)

        self = cls()

        self.bin_id = cls.file_consume_bin_id(ifs)
        for i in range(cls.file_consume_num_chunks(ifs)):
            self.chunks.append((
                cls.file_consume_beg_off(ifs),
                cls.file_consume_end_off(ifs),
            ))

        return self

    def to_dict(self):
        """
        Returns a dict representation of the BAM Index Bin.
        """
        return {
            'bin_id': self.bin_id,
            'chunks': self.chunks,
        }

    @classmethod
    def from_dict(cls, d):
        """
        Generates a Bin instance from a dict object.
        """
        self = cls()

        self.bin_id = d['bin_id']
        self.chunks = d['chunks']

        return self


class Header(object):
    """
    BAM Header
    https://github.com/LabAdvComp/python-bam/blob/master/bam/bam.py
    """

    def __init__(self):
        self.refs = list()
        self.lens = list()

    @classmethod
    def from_file(cls, filename=None, mode='rb', fileobj=None):

        assert any([
            filename is not None,
            fileobj is not None,
        ]), 'Must specify filename or fileobj.'

        assert None in [
            filename,
            fileobj,
        ], 'Filename and fileobj are mutually exclusive.'

        ifs = fileobj if fileobj is not None else gzip.GzipFile(filename, mode)

        self = cls()

        try:

            if ifs.read(4) != 'BAM\x01':
                raise ValueError('BAM magic number not found')

            hlen = deserialize('<i', ifs)
            self.head = ifs.read(hlen)

            if len(self.head) < hlen:
                raise ValueError('BAM header is truncated')

            # TODO once the SAM header is parsed, this can be changed to a map
            rlen = deserialize('<i', ifs)
            for i in range(rlen):
                nlen = deserialize('<i', ifs)
                name = deserialize('%dsx' % (nlen - 1), ifs)
                size = deserialize('<i', ifs)

                self.refs.append(name.decode('utf-8'))
                self.lens.append(size)

        except struct.error:

            raise ValueError('BAM header is truncated')

        return self

    def to_file(self, filename=None, mode='wb', fileobj=None):

        assert any([
            filename is not None,
            fileobj is not None,
        ]), 'Must specify filename or fileobj.'

        assert None in [
            filename,
            fileobj,
        ], 'Filename and fileobj are mutually exclusive.'

        ofs = fileobj if fileobj is not None else open(filename, mode)

        ofs.write(self.to_bytes())

    @classmethod
    def from_bytes(cls, b):
        raise NotImplementedError('TODO')

    def to_bytes(self):

        header = bytes()

        header += b'BAM\x01'
        header += struct.pack('<i', len(self.head))
        header += self.head
        header += struct.pack('<i', len(self.refs))
        for r, l in zip(self.refs, self.lens):
            header += struct.pack('<i', len(r) + 1)
            header += struct.pack('<%ds' % len(r), r.encode('utf-8'))
            header += b'\x00'
            header += struct.pack('<i', l)

        return header