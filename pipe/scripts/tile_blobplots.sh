#!/usr/bin/env bash
tmpdir=/tmp/tiled_blobplots
mkdir $tmpdir

for img in output/produce_covplots_reduced/*.covsum.png; do
    convert "$img" -resize 300x300 $tmpdir/"$(basename "$img")"
done

montage "$tmpdir/*.covsum.png" \
    -scale 1000x1000 \
    -mode Concatenate \
    -tile 3x3 covplots_tiled.png

rm -rf $tmpdir

montage output/produce_blobplots_reduced/*/phylum/*blobplot.covsum.png \
    -scale 1000x1000 \
    -mode Concatenate \
    -tile 3x3 blobplots_tiled.png
