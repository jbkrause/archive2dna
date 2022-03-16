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

import os
import pprint
import configparser
import argparse

from archive2dna import package
from archive2dna import dna as dna_module

pp = pprint.PrettyPrinter(depth=6)

# Parse arguments
usage = """cli.py ACTION FILE_IN FILE_OUT [--id package_id] [--er error_rate] [--config config_section]
           action : encode | decode | corrupt\n"""
parser = argparse.ArgumentParser(
    description="Encode/decode information package to DNA."
)
parser.add_argument("action", help="encode | decode | corrupt .")
parser.add_argument("input_file", help="Input file.")
parser.add_argument("output_file", help="Output file.")
parser.add_argument(
    "--id",
    dest="package_id",
    help="Information package ID, used to generate the primer",
)
parser.add_argument(
    "--er",
    dest="error_rate",
    type=float,
    help="Error rate ER in perecentage, used to corrupt the DNA",
)
parser.add_argument(
    "--config", help="Config set to be used, e.g. DEFAULT or BIG (see config.ini)"
)
args = parser.parse_args()

# read config
cfg = configparser.ConfigParser()
cfg.read("config.ini")
if args.config == None:
    cfg_set = "DEFAULT"
else:
    cfg_set = args.config
primer_length = int(cfg[cfg_set]["primer_length"])
mi = int(cfg[cfg_set]["mi"])
mo = int(cfg[cfg_set]["mo"])
index_length = int(cfg[cfg_set]["index_length"])
index_positions = int(cfg[cfg_set]["index_positions"])
N = int(cfg[cfg_set]["N"])
K = int(cfg[cfg_set]["K"])
target_redundancy = float(cfg[cfg_set]["target_redundancy"])

# read technical config
cfg_set = "TECHNICAL"
if cfg[cfg_set]["auto_zip"] == "False":
    auto_zip = False
else:
    auto_zip = True
representation_type = cfg[cfg_set]["representation_type"]
representation_url = cfg[cfg_set]["representation_url"]
logging_file = cfg[cfg_set]["logging_file"]
logging_level = cfg[cfg_set]["logging_level"]

if args.package_id == None:
    primer_length = 0

if args.action == "encode":
    binary = args.input_file
    dna = args.output_file
    binary_data = open(binary, "rb").read()
    c = package.Container(
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
    c.load_binary(binary_data)
    c.create_logical_redundancy()
    c.convert_to_dna()
    c.compute_segments_sizes()
    text = c.write_dna()
    open(dna, "w").write(text)
    pp.pprint(c.compute_stats())

elif args.action == "decode":
    binary = args.output_file
    dna = args.input_file
    c = package.Container(
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
    text = open(dna, "r").read()
    c.load_dna(text)
    c.check_and_correct_logical_redundancy()
    binary_data = c.write_binary()
    open(binary, "wb").write(binary_data)
    pp.pprint(c.compute_stats())

elif args.action == "corrupt":
    dna = args.input_file
    dna_out = args.output_file
    if args.error_rate is not None:
        error_rate = args.error_rate / 100
    else:
        error_rate = 0.5 / 100
    text = open(dna, "r").read().split("\n")
    corrupted_segments = 0
    number_of_corruptions = 0
    for i in range(len(text)):
        corrupted = dna_module.corrupt_dna_segment(text[i], error_rate)
        corruptions = sum(i != j for i, j in zip(corrupted, text[i]))
        if corruptions > 0:
            corrupted_segments += 1
            number_of_corruptions += corruptions
            text[i] = corrupted
    open(dna_out, "w").write("\n".join(text))
    print("Corrupted segments :", corrupted_segments, "/", len(text))
    print("Total number of corruptions :", number_of_corruptions)

else:
    print(usage)
