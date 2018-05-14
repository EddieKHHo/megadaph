#!/usr/bin/env python
"""Find variants which are unique to sample from a multisample gatk variant
table.
"""
from pandas import DataFrame


class VariantTable:
    def __init__(self, df):
        self.df = df

    @property
    def allele_counts(self):
        """Count number of alternative alleles per sample in gatk variant table.

        Returns
        -------
        List[DataFrame]
            A table of allele counts for each sample in `df`

        """
        allele_counts = [DataFrame(self.df[col]
                                   .str.split(',')
                                   .tolist()
                                   ).fillna(0).astype(int)
                         for col in self.adcols]
        return allele_counts

    def __len__(self):
        return len(self.df.values)

    @property
    def coverage(self):
        coverage = [self.df[col] for col in self.dpcols]
        return coverage

    @property
    def adcols(self):
        adcols = [x for x in self.df if '.AD' in x]
        return adcols

    @property
    def dpcols(self):
        dpcols = [x for x in self.df if '.DP' in x]
        return dpcols

    @property
    def nsamples(self):
        return len(self.adcols)

    @property
    def samples(self):
        return [colname.split('.')[0] for colname in self.adcols]

    @property
    def allele_freq(self):
        allele_freq = [counts.divide(coverage.values, axis='index')
                       for counts, coverage
                       in zip(self.allele_counts, self.coverage)]
        return allele_freq
