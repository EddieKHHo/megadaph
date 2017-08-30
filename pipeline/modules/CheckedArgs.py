#!/usr/bin/env python
import os
from fen_util import any_dont_exist, prefix
from glob import glob

class CheckedArgs:
    def __init__(self, args):
        self.assemblies = self.set_assemblies(args)
        self.fwd_reads, self.rev_reads = self.match_reads(args)

    # Set the assemblies slot
    # Param:
    #   args    A docopt argument dictionary
    # Return:
    #   A list of assembly files if they are valid. Otherwise, raise an 
    #   exception
    def set_assemblies(self, args):
        assemblies = sorted(args['<assembly>'])

        if any_dont_exist(assemblies): 
            raise ValueError("Invalid assembly argument") 
        else:
            return assemblies

    # Set the fwd_reads and rev_reads slots
    # Param:
    #   args    A docopt argument dictionary
    # Return:
    #   A tuple of lists (fwd_read, rev_reads). Exception if any reads cannot 
    #   be found.
    def match_reads(self, args):
        read_dir = args['<read_dir>']
        prefixes = map(prefix, self.assemblies)
        prefixes = map(os.path.basename, prefixes)

        fwd = []
        rev = []
        for pre in prefixes:
            prefixed_read_dir = os.path.join(read_dir, pre)

            if os.path.exists(prefixed_read_dir):

                fwd_reads = glob(prefixed_read_dir + "/" + pre + ".R1.f*")
                rev_reads = glob(prefixed_read_dir + "/" + pre + ".R2.f*")

                if any_dont_exist(fwd_reads + rev_reads):
                    raise ValueError("Reads missing from " + prefixed_read_dir)
                else:
                    fwd = fwd + fwd_reads
                    rev = rev + rev_reads
            else:
                raise ValueError(prefixed_read_dir + " does not exist")

        return (sorted(fwd), sorted(rev))


