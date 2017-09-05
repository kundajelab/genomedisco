#!/usr/bin/env python

import sys
import os

import hifive
import h5py
import numpy
import gzip

class Encode_Data(hifive.HiCData):

    def load_data_from_int(self, fendfilename, filelist):
        self.history += "HiCData.load_data_from_int(fendfilename='%s', filelist=%s, ) - " % (fendfilename, str(filelist))
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, \
                ("The fend file %s was not found. No data was loaded.\n") % (fendfilename),
            self.history += "Error: '%s' not found\n" % fendfilename
            return None
        self.fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                       os.path.dirname(self.file)), os.path.basename(fendfilename))
        self.fends = h5py.File(fendfilename, 'r')
        if 'binned' in self.fends['/'].attrs and self.fends['/'].attrs['binned'] is not None:
            self.binned = True
        else:
            self.binned = False
        if 'fends' in self.fends and self.fends['fends'] is not None:
            self.re = True
        else:
            self.re = False
        self.history = self.fends['/'].attrs['history'] + self.history
        self.chr2int = {}
        chroms = self.fends['chromosomes'][...]
        for i, j in enumerate(chroms):
            self.chr2int[j] = i
        self.insert_distribution = numpy.zeros((182, 2), dtype=numpy.int32)
        self.insert_distribution[1:, 1] = numpy.round(numpy.exp(numpy.linspace(3.8, 12.8, 181))).astype(numpy.int32)
        # load data from all files, skipping if chromosome not in the fend file.
        if isinstance(filelist, str):
            filelist = [filelist]
        int_filelist = []
        for filename in filelist:
            int_filelist.append("%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(filename)),
                                           os.path.dirname(self.file)), os.path.basename(filename)))
        self.int_filelist = ",".join(int_filelist)
        fend_pairs = []
        for i in range(len(chroms)):
            fend_pairs.append([])
            for j in range(i + 1):
                fend_pairs[i].append({})
        total_reads = 0
        for fname in filelist:
            data = []
            for i in range(len(chroms)):
                data.append([])
                for j in range(i + 1):
                    data[i].append({})
            if not os.path.exists(fname):
                if not self.silent:
                    print >> sys.stderr, ("The file %s was not found...skipped.\n") % (fname.split('/')[-1]),
                self.history += "'%s' not found, " % fname
                continue
            if not self.silent:
                print >> sys.stderr, ("Loading data from %s...") % (fname.split('/')[-1]),
            a = 0
            input = gzip.open(fname, 'r', 1)
            new_reads = 0
            for line in input:
                temp = line.strip('\n').split('\t')
                temp[0] = temp[0].strip('chr')
                temp[2] = temp[2].strip('chr')
                if temp[0] not in self.chr2int or temp[2] not in self.chr2int:
                    self.stats['chr_not_in_fends'] += 1
                    #print 'dont like chromo'
                    continue
                chrint1 = self.chr2int[temp[0]]
                chrint2 = self.chr2int[temp[2]]
                start1 = int(temp[1])
                start2 = int(temp[3])
                if chrint1 == chrint2:
                    if start1 > start2:
                        start1, start2 = start2, start1
                elif chrint2 < chrint1:
                    chrint1, chrint2, start1, start2 = chrint2, chrint1, start2, start1
                key = (start1, start2)
                if key not in data[chrint2][chrint1]:
                    data[chrint2][chrint1][key] = 1
                else:
                    data[chrint2][chrint1][key] += 1
                a += 1
                if a >= 10000000:
                    for i in range(len(data)):
                        for j in range(len(data[i])):
                            temp = numpy.zeros((len(data[i][j]), 3), dtype=numpy.int32)
                            k = 0
                            for key, count in data[i][j].iteritems():
                                temp[k, :] = (key[0], key[1], count)
                                self.stats['pcr_duplicates'] += count - 1
                                k += 1
                            data[i][j] = temp
                            new_reads += numpy.sum(data[i][j][:, 2])
                    if self.re:
                        self._find_fend_pairs(data, fend_pairs, True)
                    else:
                        self._find_bin_pairs(data, fend_pairs, True)
                    if not self.silent:
                        print >> sys.stderr, ("\r%s\rLoading data from %s...") % (' '*50, fname.split('/')[-1]),
                    for i in range(len(data)):
                        for j in range(len(data[i])):
                            data[i][j] = {}
                    a = 0
            input.close()
            if a > 0:
                for i in range(len(data)):
                    for j in range(len(data[i])):
                        temp = numpy.zeros((len(data[i][j]), 3), dtype=numpy.int32)
                        k = 0
                        for key, count in data[i][j].iteritems():
                            temp[k, :] = (key[0], key[1], count)
                            self.stats['pcr_duplicates'] += count - 1
                            k += 1
                        data[i][j] = temp
                        new_reads += numpy.sum(data[i][j][:, 2])
                # map data to fends, filtering as needed
                if new_reads > 0 and a > 0:
                    if self.re:
                        self._find_fend_pairs(data, fend_pairs, True)
                    else:
                        self._find_bin_pairs(data, fend_pairs, True)
            total_reads += new_reads
            if not self.silent:
                print >> sys.stderr, ("\r%s\r%i validly-mapped reads pairs loaded.\n") % (' ' * 50, new_reads),
        self.stats['total_reads'] = total_reads + self.stats['chr_not_in_fends']
        if self.re:
            self._clean_fend_pairs(fend_pairs)
        total_fend_pairs = 0
        for i in range(len(fend_pairs)):
            for j in range(len(fend_pairs[i])):
                total_fend_pairs += len(fend_pairs[i][j])
        if total_fend_pairs == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            self.history += "Error: no valid data loaded\n"
            return None
        if not self.silent:
            print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i valid fend pairs\n") %\
                             (total_reads, total_fend_pairs),
        # write fend pairs to h5dict
        if self.re and self.binned:
            self._parse_binned_fend_pairs(fend_pairs)
        else:
            self._parse_fend_pairs(fend_pairs)
        self.history += 'Success\n'
        return None

def main():
    int_fname, fend_fname, out_fname = sys.argv[1:4]
    data = Encode_Data(out_fname, 'w')
    data.load_data_from_int(fend_fname, int_fname)
    data.save()


if __name__ == "__main__":
    main()
