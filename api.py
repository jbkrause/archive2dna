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

import json
import io

from archive2dna import package

# import main Flask class and request
from flask import Flask, request, jsonify, send_file

# parameters
api_version = "0.0.1"
api_name = "archive2dna"
response_template = {"api": {"name": api_name, "version": api_version}}

# create the Flask app
app = Flask(__name__)


@app.route("/", methods=["GET"])
def process_root():
    if request.method == "GET":
        """general documentation"""
        response = response_template
        us = request.url.split("/")
        u = us[0] + "://" + us[2] + "/"
        response[
            "description"
        ] = "Encodes binary file to DNA and decodes back from DNA to binary."
        response["routes"] = {
            "encode": {
                "route": u + "encode",
                "description": "Encodes binary file to DNA.",
            },
            "decode": {
                "route": u + "decode",
                "description": "Decodes DNA file to binary.",
            },
        }
    return jsonify(response)


@app.route("/encode", methods=["GET", "POST"])
def process_encode():
    if request.method == "GET":
        """encode route documentation"""
        response = response_template
        response["description"] = "Encodes binary file to DNA."
        response["request-verb"] = "POST"
        response["request-type"] = "multipart"
        response["request-parameters"] = {
            "data": "binary file",
            "id": "information package id [optionnal]",
        }
        response[
            "example"
        ] = 'curl -F id="test:1" -F data="@aip.zip" -X POST http://localhost:8080/encode'

    if request.method == "POST":
        """encode binary to DNA"""
        package_id = request.form.get("id")
        data = request.files.get("data").read()
        c = package.Container(package_id=package_id, primer_length=5)
        c.load_binary(data)
        c.create_logical_redundancy()
        c.convert_to_dna()
        c.compute_segments_sizes()
        text = c.write_dna()
        response = response_template
        response["encodedDNA"] = text.split("\n")
        response["statistics"] = c.compute_stats()

    return jsonify(response)


@app.route("/decode", methods=["GET", "POST"])
def process_decode():
    if request.method == "GET":
        """decode route documentation"""
        response = response_template
        response[
            "description"
        ] = "Decodes DNA text file to binary file. DNA text file consists of one DNA segment per line (composed of letters A, T, G and C) and is encoded in utf-8."
        response["request-verb"] = "POST"
        response["request-parameters"] = {
            "data": "dna file",
            "id": "information package id [optionnal]",
        }
        response[
            "example"
        ] = "curl --data-biinary @dna.txt -X POST http://localhost:8080/decode -o aip_decoded.zip"
        return jsonify(response)

    if request.method == "POST":
        """decode DNA to binary"""
        dna = request.get_data().decode("utf-8")
        c = package.Container(package_id="does_not_matter", primer_length=5)
        c.load_dna(str(dna))
        c.check_and_correct_logical_redundancy()
        output_bytes = c.write_binary()
        f = io.BytesIO(output_bytes)
        return send_file(f, attachment_filename="package.zip", as_attachment=False)
        # return send_file(f, attachment_filename='package.zip' ,as_attachment=True)


if __name__ == "__main__":
    # run app in debug mode on port 8080
    app.run(debug=True, port=8080)
