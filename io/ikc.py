"""
Indexed k-mer count (IKC) utility classes and functions.
"""

import mmap
import os


class IKCHeaderV1:
    """
    Create a new version 1 IKC header.
    """

    def __init__(self, k_size, k_min_size, k_min_mask, offset_data, offset_index, offset_meta, offset_eof, id_string):
        self.k_size = k_size
        self.k_min_size = k_min_size
        self.k_min_mask = k_min_mask
        self.offset_data = offset_data
        self.offset_index = offset_index
        self.offset_meta = offset_meta
        self.offset_eof = offset_eof
        self.id_string = id_string

        self.kmer_util = util.kmer.KmerUtil(k_size, k_min_size, k_min_mask)


def get_header(mmap_file):
    """
    Get a header from an open file.

    :param mmap_file: Memory-mapped (mmap) file.

    :return: Header object.
    """

    # Check
    if mmap_file is None:
        raise RuntimeError('Cannot read header from file channel: None')

    # Read buffer
    mmap_file.seek(0)

    # Check magic
    header_magic = mmap_file.read(15)\

    if header_magic != b'Idx_Kmer_Count\0':
        raise RuntimeError('Bad file magic: The first 15 bytes of a valid IKC file should be "Idx_Kmer_Count\\0"')

    # Get version
    file_version = get_int_field(mmap_file, 1)

    if file_version != 1:
        raise RuntimeError('File version is incorrect: Currently supported version is 1: Found {}'.format(file_version))

    return get_ikc_header_v1(mmap_file)


def get_ikc_header_v1(mmap_file):
    """
    Get a version 1 IKC header. `mmap_file` is positioned after file magic and version.

    :param mmap_file: Memory-mapped file.

    :return: Header object
    """

    # Check reserved fields (should be 0)
    for index in range(7):
        reserved_byte = mmap_file.read_byte()

        if reserved_byte != 0:
            raise RuntimeError('IKC reserved field {} is not null: 0x{:x}'.format(index + 1, reserved_byte))

    # Get minimizer size
    k_min_size = get_int_field(mmap_file, 1)

    if not (0 < k_min_size < 16):
        raise RuntimeError(
            'Minimizer size must be between 0 and 15 (two-complement byte field): {:d}'.format(k_min_size)
        )

    # Get k-mer size
    k_size = get_int_field(mmap_file, 4)

    if k_size < 1:
        raise RuntimeError('IKC file has bad k-mer size: {}'.format(k_size))

    # Get k-mer mask
    k_min_mask = get_int_field(mmap_file, 4)

    # Get offsets
    offset_index = get_int_field(mmap_file, 8)
    offset_meta = get_int_field(mmap_file, 8)

    # Get optional ID
    id_string = mmap_file.read(32).split(b'\x00')[0]

    # Get header length (also the data offset)
    offset_data = mmap_file.tell()

    # Get file length
    mmap_file.seek(0, os.SEEK_END)
    offset_eof = mmap_file.tell()

    # Return header
    return IKCHeaderV1(k_size, k_min_size, k_min_mask, offset_data, offset_index, offset_meta, offset_eof, id_string)


def get_int_field(mmap_file, n_bytes, signed=True):

    return int.from_bytes(
        mmap_file.read(n_bytes),
        'big',
        signed=signed
    )


class IKCReader:
    """
    Read from an IKC (indexed k-mer count) file.
    """

    def __init__(self, ikc_file_name):

        self.__ikc_file = None
        self.__mmap_file = None

        # Init
        self.ikc_file_name = ikc_file_name

        self.__ikc_file = open(self.ikc_file_name, 'r')

        self.__mmap_file = mmap.mmap(
            self.__ikc_file.fileno(), 0,
            mmap.MAP_PRIVATE, mmap.PROT_READ
        )

        # Get header
        self.ikc_header = get_header(self.__mmap_file)

        self.__offset_data = self.ikc_header.offset_data
        self.__offset_index = self.ikc_header.offset_index
        self.__offset_meta = self.ikc_header.offset_meta
        self.__offset_eof = self.ikc_header.offset_eof

        self.__kmer_util = self.ikc_header.kmer_util
        self.__min_kmer_util = self.__kmer_util.min_kmer_util

        self.__kmer_bytes = self.__kmer_util.word_size * 4
        self.__count_bytes = 4

        # Declare variables set by __read_index()
        self.__last_min = None
        self.__last_min_offset = None
        self.__last_min_n = None

        # Read index: minimizer-keyed index to tuples of (offset, length)
        self.__read_index()

    def __del__(self):
        """
        Close files.
        """

        try:
            if self.__mmap_file is not None:
                self.__mmap_file.close()
                self.__mmap_file = None

        finally:
            if self.__ikc_file is not None:
                self.__ikc_file.close()
                self.__ikc_file = None

    def __read_index(self):
        """
        Read the IKC file index.
        """

        # Init
        mmap_file = self.__mmap_file

        self.__min_index = dict()

        end_location = self.__offset_meta

        kmer_bytes = self.__kmer_util.word_size_bytes  # Number of bytes per k-mer
        record_size = kmer_bytes + 4  # Number of bytes per data record (k-mer and a 4-byte count)

        # Check for 0 records
        if self.__offset_index == end_location:
            raise RuntimeError('No index in IKC file')

        # Read first record
        mmap_file.seek(self.__offset_index)

        minimizer = get_int_field(mmap_file, 4)
        offset = get_int_field(mmap_file, 8)

        min_count = 1

        if offset != self.__offset_data:
            raise RuntimeError(
                'Minimizer #{} in index (0x{:x}) does not point to the data offset (0x{:x}): 0x{:x}'.format(
                    min_count, minimizer, self.__offset_data, offset
                )
            )

        # Seed last minimizer - First search starts here
        self.__last_min = minimizer

        # Read remaining records
        while mmap_file.tell() < end_location:
            min_count += 1

            # Read next minimizer
            next_minimizer = get_int_field(mmap_file, 4, signed=False)
            next_offset = get_int_field(mmap_file, 8, signed=False)

            # Check bounds
            if not (offset < next_offset < end_location):
                raise RuntimeError((
                    'Minimizer #{} in index (0x{:x}) must be between the last offset (0x{:x}) and the end of the data '
                    'section (0x{:x}): 0x{:x}').format(
                    min_count, next_minimizer, offset, end_location, next_offset
                ))

            # Write record for the last minimizer (needed offset from this one to compute length)
            group_len = next_offset - offset

            if group_len % record_size != 0:
                raise RuntimeError((
                    'Size of minimizer group #{} in index (0x{:x}) is not a modulus of the k-mer word size (32-bit '
                    'words) and count (32-bit count): {} !% {}').format(
                    min_count, next_minimizer, group_len, record_size
                ))

            # Add to minimizer dict
            if minimizer in self.__min_index:
                raise RuntimeError((
                    'Found multiple entries for minimizer (0x{:x}): First duplicate at record #{} in index').format(
                    next_minimizer, min_count
                ))

            self.__min_index[minimizer] = (offset, group_len // record_size)

            # Setup for next record
            minimizer = next_minimizer
            offset = next_offset

        # Finish seeding the last minimizer (now that the length is known)
        self.__last_min_offset, self.__last_min_n = self.__min_index[self.__last_min]

    def get(self, kmer):
        """
        Get the count of `kmer`.

        :param kmer: Integer k-mer.

        :return: Count of `kmer` or `0` if the k-mer is not in the IKC file.
        """

        # Get count from last minimizer set (avoids recomputing minimizer)
        count = self._search(kmer, self.__last_min_offset, self.__last_min_n)

        if count > 0:
            return count

        # Get minimizer for this k-mer
        minimizer = self.__kmer_util.minimizer(kmer)

        if minimizer == self.__last_min or minimizer not in self.__min_index:
            return 0

        self.__last_min = minimizer
        self.__last_min_offset, self.__last_min_n = self.__min_index[minimizer]

        return self._search(kmer, self.__last_min_offset, self.__last_min_n)

    def _search(self, kmer, offset, n):
        """
        Search the IKC file for a k-mer.

        :param kmer: K-mer.
        :param offset: Minimizer group offset.
        :param n: Number of k-mers in the minimizer group.

        :return: K-mer count or `0` if the k-mer was not found.
        """

        print('Searching: 0x{:x} (offset = 0x{:d}, n = {:d})'.format(kmer, offset, n))

        # Init binary search.
        first = 0
        last = n

        mmap_file = self.__mmap_file

        kmer_bytes = self.__kmer_bytes
        count_bytes = self.__count_bytes
        record_bytes = kmer_bytes + count_bytes

        # Run binary search
        while first <= last:

            # Get k-mer
            mid = (first + last) // 2

            mmap_file.seek(offset + (mid * record_bytes))

            next_kmer = get_int_field(mmap_file, kmer_bytes, signed=False)

            print('\t* 0x{:x} (mid = {:d})'.format(next_kmer, mid))

            if next_kmer == kmer:
                return get_int_field(mmap_file, count_bytes)

            # Advance search
            if kmer > next_kmer:
                first = mid + 1
            else:
                last = mid - 1

        # No match found
        return 0

    def index_items(self):
        """
        Get a list of items in the index. Each is a tuple of "(minimizer, (offset, length))", where length is the number
        of records and offset is the file offset in bytes to the first record.

        :return:
        """
        return self.__min_index.items()

    def iter_min_order(self):
        """
        Iterate over elements in the IKC file in minimizer order then k-mer order.

        :return: An iterator of tuples (k-mer, count) in minimizer/k-mer sort order.
        """

        mmap_file = self.__mmap_file

        kmer_bytes = self.__kmer_bytes
        count_bytes = self.__count_bytes

        for minimizer in sorted(self.__min_index.keys()):
            # Init minimizer group
            offset, n_elements = self.__min_index[minimizer]

            mmap_file.seek(offset)

            for index in range(n_elements):
                yield (
                    get_int_field(mmap_file, kmer_bytes, signed=False),
                    get_int_field(mmap_file, count_bytes)
                )

    def iter_kmer_order(self):
        """
        Iterate records in k-mer sort order.

        :return: An iterator for records in k-mer sort order.
        """

        mmap_file = self.__mmap_file

        kmer_bytes = self.__kmer_bytes
        count_bytes = self.__count_bytes
        record_bytes = kmer_bytes + count_bytes

        # Track minimum k-mer per minimizer and next k-mer (minimum of all)
        group_dict = dict()  # Keyed by minimizer, tuple of (k-mer, count, offset, last_offset)

        min_kmer = None
        min_count = None
        min_minimizer = None

        for minimizer in self.__min_index.keys():
            offset, n_elements = self.__min_index[minimizer]
            last_offset = offset + (n_elements + 1) * record_bytes

            mmap_file.seek(offset)

            next_kmer = get_int_field(mmap_file, kmer_bytes, signed=False)
            next_count = get_int_field(mmap_file, count_bytes)

            group_dict[minimizer] = (
                next_kmer,
                next_count,
                offset,
                last_offset,
            )

            if min_kmer is None or next_kmer < min_kmer:
                min_kmer = next_kmer
                min_count = next_count
                min_minimizer = minimizer

        while min_kmer is not None:

            # Cache values for this iteration (returned at end)
            this_kmer = min_kmer
            this_count = min_count

            # Advance last minimizer group
            next_kmer, next_count, offset, last_offset = group_dict[min_minimizer]

            offset += record_bytes

            if offset < last_offset:
                mmap_file.seek(offset)

                next_kmer = get_int_field(mmap_file, kmer_bytes, signed=False)
                next_count = get_int_field(mmap_file, count_bytes)

                group_dict[min_minimizer] = (
                    next_kmer,
                    next_count,
                    offset,
                    last_offset
                )

            else:
                del(group_dict[min_minimizer])

            # Find the next minimum k-mer
            min_kmer = None
            min_count = None
            min_minimizer = None

            for minimizer in group_dict.keys():
                next_kmer, next_count, offset, last_offset = group_dict[minimizer]

                if min_kmer is None or next_kmer < min_kmer:
                    min_kmer = next_kmer
                    min_count = next_count
                    min_minimizer = minimizer

            # Yield minimum k-mer for this round
            yield this_kmer, this_count
