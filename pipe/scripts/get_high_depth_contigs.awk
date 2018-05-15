#!/usr/bin/gawk -f
## Extract scaffold names which do not have low (< 3) coverage from output
## of blobtools map2cov
{
    if (NR > 6) {
        if ($3 > 3) {
            print $1;
        }
    }
}
