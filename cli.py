#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of archive2dna.
#
# archive2dna is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Foobar is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>
#
# Author : Jan Krause-Bilvin
# First release: 2022-02-02

import pprint
import configparser
import argparse
import sys

from archive2dna import package
from archive2dna import dna as dna_module

pp = pprint.PrettyPrinter(depth=6, stream=sys.stderr)


def encode(args):
    binary_data = args.infile.read()
    container = createContainer(args)
    container.load_binary(binary_data)
    container.create_logical_redundancy()
    container.convert_to_dna()
    container.compute_segments_sizes()
    text = container.write_dna()
    args.outfile.write(text)
    pp.pprint(container.compute_stats())


def decode(args):
    dna_data = args.infile.read()
    container = createContainer(args)
    container.load_dna(dna_data)
    container.check_and_correct_logical_redundancy()
    binary_data = container.write_binary()
    args.outfile.write(binary_data)
    pp.pprint(container.compute_stats())


def corrupt(args):
    error_rate = args.error_rate / 100
    text = args.infile.read().split("\n")
    corrupted_segments = 0
    number_of_corruptions = 0
    for i in range(len(text)):
        corrupted = dna_module.corrupt_dna_segment(text[i], error_rate)
        corruptions = sum(i != j for i, j in zip(corrupted, text[i]))
        if corruptions > 0:
            corrupted_segments += 1
            number_of_corruptions += corruptions
            text[i] = corrupted
    args.outfile.write("\n".join(text))
    print("Corrupted segments :", corrupted_segments, "/", len(text), file=sys.stderr)
    print("Total number of corruptions :", number_of_corruptions, file=sys.stderr)


def createContainer(args):
    # read config
    config = configparser.ConfigParser()
    config.read("config.ini")
    section = config[args.config]
    primer_length = int(section["primer_length"])
    mi = int(section["mi"])
    mo = int(section["mo"])
    index_length = int(section["index_length"])
    index_positions = int(section["index_positions"])
    N = int(section["N"])
    K = int(section["K"])
    target_redundancy = float(section["target_redundancy"])

    # read technical config
    technical = config["TECHNICAL"]
    auto_zip = not (technical["auto_zip"] == "False")
    representation_type = technical["representation_type"]
    representation_url = technical["representation_url"]
    logging_file = technical["logging_file"]
    logging_level = technical["logging_level"]

    if args.package_id == None:
        primer_length = 0

    return package.Container(
        package_id=args.package_id,
        primer_length=primer_length,
        mi=mi,
        mo=mo,
        index_length=index_length,
        index_positions=index_positions,
        N=N,
        K=K,
        target_redundancy=target_redundancy,
        representation_type=representation_type,
        representation_url=representation_url,
        logging_file=logging_file,
        logging_level=logging_level,
        auto_zip=auto_zip,
    )


def ranged_float(min, max):
    """Returns an argument type function for ArgumentParser checking a float
    with a range between min and max."""

    def ranged_float_checker(arg):
        try:
            f = float(arg)
        except ValueError:
            raise argparse.ArgumentTypeError("must be a floating point number")
        if f < min or f > max:
            raise argparse.ArgumentTypeError(
                f"must be in the range [{str(min)} .. {str(max)}]"
            )
        return f

    return ranged_float_checker


# Main parser.
parser = argparse.ArgumentParser(
    description="Encode/decode information package to DNA."
)
parser.add_argument(
    "--config",
    default="DEFAULT",
    help="Config set to be used, e.g. DEFAULT or BIG (see config.ini)",
)
subparsers = parser.add_subparsers(required=True)

# Encode parser.
encode_parser = subparsers.add_parser("encode", help="encode binary data into dna")
encode_parser.set_defaults(func=encode)
encode_parser.add_argument(
    "--id",
    dest="package_id",
    help="Information package ID, used to generate the primer",
)
encode_parser.add_argument(
    "infile",
    nargs="?",
    type=argparse.FileType("rb"),
    default=sys.stdin.buffer,
    help="input binary file",
)
encode_parser.add_argument(
    "outfile",
    nargs="?",
    type=argparse.FileType("w"),
    default=sys.stdout,
    help="output dna formatted file",
)

# Decode parser.
decoder_parser = subparsers.add_parser(
    "decode", help="decode dna back into the original binary representation"
)
decoder_parser.set_defaults(func=decode)
decoder_parser.add_argument(
    "--id",
    dest="package_id",
    help="Information package ID, used to generate the primer",
)
decoder_parser.add_argument(
    "infile",
    nargs="?",
    type=argparse.FileType("r"),
    default=sys.stdin,
    help="input dna formatted file",
)
decoder_parser.add_argument(
    "outfile",
    nargs="?",
    type=argparse.FileType("wb"),
    default=sys.stdout.buffer,
    help="output binary file",
)

# Corrupt parser.
corrupt_parser = subparsers.add_parser(
    "corrupt", help="corrupt dna for testing purposes"
)
corrupt_parser.set_defaults(func=corrupt)
corrupt_parser.add_argument(
    "-e",
    "--error-rate",
    default=0.5,
    dest="error_rate",
    type=ranged_float(0, 100),
    help="Error rate ER in percentage, used to corrupt the DNA",
)
corrupt_parser.add_argument(
    "infile",
    nargs="?",
    type=argparse.FileType("r"),
    default=sys.stdin,
    help="input dna formatted file",
)
corrupt_parser.add_argument(
    "outfile",
    nargs="?",
    type=argparse.FileType("w"),
    default=sys.stdout,
    help="output corrupted dna formatted file",
)

args = parser.parse_args()
args.func(args)
