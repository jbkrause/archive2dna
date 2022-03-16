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

import hashlib


def sha256(filename):
    """Computes sha256 checksum of file."""
    h = hashlib.sha256()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(4096), b""):
            h.update(block)
    return h.hexdigest()


def split_bytes_in_four(b):
    """Split each byte in four bytes by pairs of bits, corresponding to a DNA base each"""
    b2 = bytes()
    for x in b:
        bits = bin(x).lstrip("0b").zfill(8)
        b2 += int(bits[0:2], 2).to_bytes(1, byteorder="big")
        b2 += int(bits[2:4], 2).to_bytes(1, byteorder="big")
        b2 += int(bits[4:6], 2).to_bytes(1, byteorder="big")
        b2 += int(bits[6:8], 2).to_bytes(1, byteorder="big")
    return b2


def merge_four_bytes_in_one(b):
    """Merges groups of 4 bytes together (taking the 2 lower bits of each)."""
    b2 = bytes()
    for i in range(len(b) // 4):
        bits = ""
        for j in range(4):
            bitsj = bin(b[4 * i + j]).lstrip("0b").zfill(8)
            bits += bitsj[6:8]
        b2 += int(bits, 2).to_bytes(1, byteorder="big")
    return b2
