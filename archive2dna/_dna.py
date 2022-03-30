#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of archive2dna.
#
# archive2dna is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Foobar is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with archive2dna. If not, see <https://www.gnu.org/licenses/>
#
# Author : Jan Krause-Bilvin
# First release: 2022-02-02

from __future__ import annotations
from typing import SupportsIndex

# The order of bases in the following string is important as it determines the
# byte representation of each base. By defining them in this order, the
# base complement becomes a binary ones' complement.
#   A => 0b00
#   G => 0b01
#   C => 0b10
#   T => 0b11
_bases = "AGCT"

# Mapping between DNA sequences of length 4 and bytes.
_bytes2dna_dict = [
    f"{b1}{b2}{b3}{b4}"
    for b1 in _bases
    for b2 in _bases
    for b3 in _bases
    for b4 in _bases
]
_dna2bytes_dict = {dna: _bytes2dna_dict.index(dna) for dna in _bytes2dna_dict}


class Dna:
    _data: bytearray

    @classmethod
    def from_bytes(cls, dna: bytes) -> Dna:
        """Creates a DNA sequence from a binary representation where each pair
        of bits is converted using the following mapping:
          * 0b00 => A
          * 0b11 => T
          * 0b01 => G
          * 0b10 => C

        Args:
            dna (str): The DNA sequence binary representation.

        Returns:
            Dna: The newly created Dna class instance.
        """
        assert isinstance(dna, bytes), "dna expected to be of type bytes"
        return cls(dna)

    @classmethod
    def from_string(cls, dna: str) -> Dna:
        """Creates a DNA sequence from a textual representation with letters A,
        T, G, and C.

        Args:
            dna (str): The DNA sequence textual representation.

        Returns:
            Dna: The newly created Dna class instance.

        Remarks:
            This function expects a multiple of four bases.
        """
        assert isinstance(dna, str), "dna expected to be of type str"
        assert len(dna) % 4 == 0, "dna expected to be a multiple of 4 bases"
        try:
            return cls(
                bytes(
                    [
                        _dna2bytes_dict[dna[index : index + 4]]
                        for index in range(0, len(dna), 4)
                    ]
                )
            )
        except KeyError:
            raise Exception("base expected to be A, T, G or C")

    def __init__(self: Dna, dna: bytearray | bytes | str | None = None):
        """Creates a DNA sequence from a textual or binary representation, in
        which case each pair of bits is converted using the following mapping:
            * 0b00 => A
            * 0b11 => T
            * 0b01 => G
            * 0b10 => C

        Args:
            dna (str): The DNA sequence binary representation.

        Returns:
            Dna: The newly created Dna class instance.
        """
        if dna == None:
            self._data = bytearray()
        elif isinstance(dna, bytearray):
            self._data = dna
        elif isinstance(dna, bytes):
            self._data = bytearray(dna)
        elif isinstance(dna, str):
            self._data = Dna.from_string(dna)._data
        else:
            raise Exception("invalid dna sequence type")

    def copy(self: Dna) -> Dna:
        return Dna(self._data.copy())

    def to_bytes(self: Dna) -> bytes:
        """Returns a binary representation of the DNA sequence where each base
        is converted into 2 bits using the following mapping:
            * A => 0b00
            * T => 0b11
            * G => 0b01
            * C => 0b10

        Returns:
            bytes: The binary representation of the DNA sequence.
        """
        return bytes(self._data)

    def to_string(self: Dna) -> str:
        """Returns a textual representation of the DNA sequence with letters A,
        T, G, and C.

        Returns:
            bytes: The textual representation of the DNA sequence.
        """
        return "".join(map(lambda b: _bytes2dna_dict[b], self._data))

    def __repr__(self: Dna) -> str:
        return self.to_string()

    def __eq__(self: Dna, other: Dna) -> bool:
        if not isinstance(other, Dna):
            return NotImplemented
        return self._data == other._data

    def __len__(self: Dna):
        return len(self._data) * 4

    def __getitem__(self: Dna, __i: SupportsIndex) -> str:
        byte = self._data[__i // 4]
        byte_offset = 6 - 2 * (__i % 4)
        bits = byte >> byte_offset & 0b11
        return _bases[bits]

    def concat(self: Dna, other: Dna) -> Dna:
        """Returns the DNA sequence resulting from concatenating the current
        sequence with the one passed as argument.

        Args:
            other (Dna): The DNA sequence to concatenate to the end.

        Returns:
            Dna: The concatenation of the current sequence with the one passwed
            as argument.
        """
        return Dna(self._data + other._data)

    def append(self: Dna, other: Dna):
        self._data += other._data

    def __add__(self: Dna, other: Dna) -> Dna:
        return self.concat(other)

    def __radd__(self: Dna, other: Dna) -> Dna:
        return other.concat(self)

    def __iadd__(self: Dna, other: Dna) -> Dna:
        self.append(other)
        return self

    def __mul__(self: Dna, __n: SupportsIndex) -> Dna:
        return Dna(self._data * __n)

    def __rmul__(self: Dna, __n: SupportsIndex) -> Dna:
        return Dna(self._data * __n)

    def __imul__(self: Dna, __n: SupportsIndex) -> Dna:
        self._data *= __n
        return self

    def complement(self) -> Dna:
        """Returns the DNA sequence that is the complement of the current
        sequence: A becomes T, T becomes A, C becomes G, and G becomes C.

        Returns:
            Dna: The complement of the current sequence.
        """
        new_data = self._data.copy()
        for index in range(len(new_data)):
            new_data[index] ^= 0xFF
        return Dna(new_data)

    def __invert__(self: Dna) -> Dna:
        return self.complement()

    def addPrimer(self: Dna, primer: Dna) -> Dna:
        """Returns a copy of the current DNA sequence surrounded by the primer.
        Note that the primer is inverted at the opposite end.

        Args:
            primer (Dna): The primer to prepend and append to the current
            sequence.

        Returns:
            Dna: The current sequence surrounded by the primer.
        """
        new_dna = primer.copy()
        new_dna += self
        new_dna += ~primer
        return new_dna

    def removePrimer(self: Dna, primer: Dna) -> Dna:
        """Returns a copy of the current DNA sequence with its primer removed.

        Args:
            primer (Dna): The primer to remove from the start and end of the
            sequence.

        Returns:
            Dna: The current sequence with its primer remove.
        """
        assert self._data.startswith(primer._data), "start primer not detected"
        primer_comp = ~primer
        assert self._data.endswith(primer_comp._data), "end primer not detected"
        return Dna(self._data.lstrip(primer._data).rstrip(primer_comp._data))
