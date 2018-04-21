#!/usr/bin/env bash
montage output/produce_blobplots3/*/phylum/*blobplot.covsum.png \
    -scale 1000x1000 \
    -mode Concatenate \
    -tile 3x3 montage.png
