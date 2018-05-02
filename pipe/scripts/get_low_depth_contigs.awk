#!/usr/bin/gawk -f
## Extract scaffold names with extremely low average coverage from the output
## of blobtools map2cov
{
    if (NR > 6) {
        if ($3 < 3) {
            print $1;
        }
    }
}
