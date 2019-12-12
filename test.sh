#!/usr/bin/env bash

docker run -it --rm -v ~/Documents/MiBuddy2_docker:/data mibuddy2 barcodes.txt -s spalax -debug
