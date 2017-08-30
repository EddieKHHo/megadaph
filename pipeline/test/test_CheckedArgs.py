#!/bin/env/python

import os, tempfile
import unittest
from glob import glob
from CheckedArgs import CheckedArgs

class CheckedArgsTestCase(unittest.TestCase):
    def test_slots(self):
        test_assembly_dir = os.path.abspath('testdat/assembly/')
        test_read_dir = os.path.abspath('testdat/reads/')
        
        # We only care about the assembly_files
        assembly_files = glob(test_assembly_dir+"/*.fasta")

        args = {
            '<assembly>' : assembly_files,
            '<read_dir>' : test_read_dir
        }

        test_class = CheckedArgs(args)
        
        expected_assemblies = sorted(['testdat/assembly/FA_SC.scaffolds.fasta', 
                                     'testdat/assembly/FB_SC.scaffolds.fasta',
                                     'testdat/assembly/IA_SC.scaffolds.fasta',
                                     'testdat/assembly/IC_SC.scaffolds.fasta'])

        expected_fwd_reads = sorted(['testdat/reads/FA_SC/FA_SC.R1.fq.gz',
                                    'testdat/reads/FB_SC/FB_SC.R1.fq.gz',
                                    'testdat/reads/IA_SC/IA_SC.R1.fq.gz',
                                    'testdat/reads/IC_SC/IC_SC.R1.fq.gz'])

        expected_rev_reads = sorted(['testdat/reads/FA_SC/FA_SC.R2.fq.gz',
                                    'testdat/reads/FB_SC/FB_SC.R2.fq.gz',
                                    'testdat/reads/IA_SC/IA_SC.R2.fq.gz',
                                    'testdat/reads/IC_SC/IC_SC.R2.fq.gz'])

        expected_assemblies = map(os.path.abspath, expected_assemblies)
        expected_fwd_reads = map(os.path.abspath, expected_fwd_reads)
        expected_rev_reads = map(os.path.abspath, expected_rev_reads)

        self.assertEqual(test_class.assemblies, expected_assemblies)
        self.assertEqual(test_class.fwd_reads, expected_fwd_reads)
        self.assertEqual(test_class.rev_reads, expected_rev_reads)

if __name__ == '__main__':
    unittest.main()

