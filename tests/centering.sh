#!/bin/sh

# Exercise manual centering on an unresampled image.  The fixture's reference
# world coordinate (0, 0) is at input pixel (129, 129), which makes output
# sizes smaller than 258 require a negative CRPIX shift and larger sizes a
# positive one.

set -eu

workdir=${TMPDIR:-/tmp}/swarp-centering-$$
mkdir "$workdir"
trap 'rm -rf "$workdir"' EXIT HUP INT TERM

fits_value()
{
  key=$1
  file=$2
  LC_ALL=C head -c 2880 "$file" | fold -w 80 | awk -v key="$key" '
    substr($0, 1, 8) ~ "^" key " *$" {
      value = substr($0, 11)
      sub(/\/.*/, "", value)
      gsub(/ /, "", value)
      print value
      exit
    }
  '
}

check_center()
{
  size_x=$1
  size_y=$2
  expected_x=$3
  expected_y=$4
  output="$workdir/coadd-${size_x}x${size_y}.fits"

  "$SWARP" "$srcdir/../test/test.fits" \
    -RESAMPLE N \
    -COMBINE Y \
    -SUBTRACT_BACK N \
    -CENTER_TYPE MANUAL \
    -CENTER 0,0 \
    -IMAGE_SIZE "$size_x,$size_y" \
    -HEADER_ONLY Y \
    -IMAGEOUT_NAME "$output" \
    -WEIGHTOUT_NAME "$workdir/weight.fits" \
    -WRITE_XML N \
    -VERBOSE_TYPE QUIET

  actual_x=$(fits_value CRPIX1 "$output")
  actual_y=$(fits_value CRPIX2 "$output")
  awk -v actual="$actual_x" -v expected="$expected_x" \
    'BEGIN { exit actual == expected ? 0 : 1 }' || {
    echo "CRPIX1 for ${size_x}x${size_y}: expected $expected_x, got $actual_x" >&2
    return 1
  }
  awk -v actual="$actual_y" -v expected="$expected_y" \
    'BEGIN { exit actual == expected ? 0 : 1 }' || {
    echo "CRPIX2 for ${size_x}x${size_y}: expected $expected_y, got $actual_y" >&2
    return 1
  }
}

# Cover negative and positive shifts on both axes, with even and odd sizes.
check_center 100 300 50 150
check_center 300 100 150 50
check_center 101 301 51 151
check_center 301 101 151 51
