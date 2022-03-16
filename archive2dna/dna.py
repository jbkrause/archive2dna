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
import array
import random

#############################
### DNA to bits and bytes ###
#############################

# Correspondence dictionary between DNA bases and bits

bits2dna_dict = {"00": "A", "01": "T", "10": "G", "11": "C"}

dna2bits_dict = {}
for kk, v in bits2dna_dict.items():
    dna2bits_dict[v] = kk


# Conversion from bytes to DNA and reeerse


def bytes2dna(b):
    """Converts a byte (type:bytes) to a DNA sequence (type:string)."""
    out = ""
    if b is not None:
        for x in b:
            tmp = bin(x).lstrip("0b").zfill(8)
            out += bits2dna_dict[tmp[0:2]]
            out += bits2dna_dict[tmp[2:4]]
            out += bits2dna_dict[tmp[4:6]]
            out += bits2dna_dict[tmp[6:8]]
        return out
    else:
        return ""


def dna2bytes(d):
    """Converts a DNA sequence (type:string) to bytes (type:bytes)."""
    if d == None or d == "None":
        return None
    out = b""
    for i in range(int(len(d) / 4)):
        bits = dna2bits_dict[d[4 * i]]
        bits += dna2bits_dict[d[4 * i + 1]]
        bits += dna2bits_dict[d[4 * i + 2]]
        bits += dna2bits_dict[d[4 * i + 3]]
        out += int(bits, 2).to_bytes(len(bits) // 8, byteorder="big")
    return out


def bits2dna(b):
    """Converts a byte containing 2 bits to DNA"""
    if b is None:
        return None
    else:
        tmp = bin(b).lstrip("0b").zfill(8)
        return bits2dna_dict[tmp[6:8]]


def dna2bits(b):
    """Converts a DNA base to two bits (in a byte)"""
    if b is None:
        return None
    else:
        bits = dna2bits_dict[b]
        return int(bits, 2)


def int2bytes(i, n=1):
    """Converts integer (0-255) to 1 byte."""
    return i.to_bytes(n, "big")


def bytes2int(b):
    """Converts integer (0-255) to 1 byte."""
    return int.from_bytes(b, byteorder="big")


def get_bytearray(a):
    """Extracts bytearray from DNA represenntations. Empty cells (None) are filtered out."""
    out = bytearray()
    for x in a:
        if x != None:  # Nones come form data padding required to fit data in 2D array
            if isinstance(x, bytes):
                out += x
            else:
                out.append(x)
    return out


# Merge and split bases in arrays


def merge_bases(a, block_size=7):
    """Converta array of single DNA bases into array of bases grouped by block_size.
    Warning: works if array length is not a multiple of block_size,
    in this case the group of last bases is padded with 0's"""
    blocks = len(a) // block_size
    if len(a) % block_size != 0:
        blocks += 1
    out = array.array("i")
    for i in range(blocks):
        # last block may be shorter
        block_range = min([block_size, len(a) - i * block_size])
        n = 0
        for j in range(block_range):
            n += a[i * block_size + j] << 2 * (block_size - j - 1)
        out.append(n)
    return out


def split_bases(a, block_size=7):
    """Converta array of DNA bases grouped by blocksize into array of single DNA bases"""
    out = array.array("i")
    for i in range(len(a)):
        n = 0
        for j in range(block_size):
            base_j = a[i] >> 2 * (block_size - j - 1)
            base_j = base_j & 3  # mask: keep two last bits
            out.append(base_j)
    return out


# Check DNA sequence vaidity


def stripDna(s):
    """Strips away character at begining and end of segment"""
    return s.strip(".,-\t ;\"'\r")


def isValidDna(s):
    """Checks wether the sting is a valid DNA sequence"""
    for c in set(s):
        if not c in "ATGC":
            return False
    return True


#########################
### Primer management ###
#########################


def complement_primer(dna):
    """Returns primer complement"""
    conv = {"A": "T", "T": "A", "G": "C", "C": "G"}
    out = ""
    for n in dna:
        out = conv[n] + out
    return out


def id2primer(package_id, length=5):
    """Computes a primer on basis of id string (encoded in utf-8).
    Primer lenght defautl of 5 bytes, i.e 5 x 4 = 20 nuclotides
    (allowing for over 1000 billion possibilities)."""
    primer_bytes = hashlib.sha256(bytes(package_id, encoding="utf-8")).digest()[
        -length:
    ]
    return bytes2dna(primer_bytes)


def add_primers(dna, primer1="", primer2=""):
    """Encapsulates a DNA sequence between two primers"""
    return primer1 + dna + primer2


def remove_primers(dna, primer1="", primer2=""):
    """Remove primers form a DNA sequence (at sequence begin and end)"""
    return dna[len(primer1) : -len(primer2)]


#############
### tests ###
#############


def corrupt_dna_segment(dna, error_rate):
    """add errors into one DNA segments"""
    conv = {"A": "T", "T": "A", "G": "C", "C": "G"}
    new_dna = list(dna)
    for i in range(len(new_dna)):
        r = random.random()
        if r < error_rate:
            new_dna[i] = conv[dna[i]]
    t = "".join(new_dna)
    return t
