
# Copyright 2008 2009 2010 2011 2013 Torbjorn Granlund, Douglas Wikstrom
#
# This file is part of GMP Modular Exponentiation Extension (GMPMEE).
#
# GMPMEE is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GMPMEE is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GMPMEE. If not, see <http://www.gnu.org/licenses/>.


import os
import sys

simexp_command = "./mpz_simexp_time"

def tostr(val):
    return str(val)

def form(l):
    return lambda s: " "*(l - len(s)) + s

class exponent_bitlen_results:
    def __init__(self):
        self.times = []

    def add(self, time):
        self.times.append(time)

    def pretty_print(self):
        min_time = min(self.times)

        for time in self.times:
            if time <= min_time + 0.00001:
                print "%6s"% ("*" + "%.2f" % time),
            else:
                print "%6.2f" % time,

    def len(self):
        return len(self.times)

    def best_block_width(self):
        m = max(self.times)
        mi = 1
        for i in range(0, len(self.times)):
            if self.times[i] < m:
                m = self.times[i]
                mi = i + 1
        return mi

    def max_block_width(self):
        return len(self.times)

class modulus_bitlen_results:
    def __init__(self):
        self.ebls = {}
        
    def add(self, exponent_bitlen, block_width_result):
        self.ebls[exponent_bitlen] = block_width_result

    def pretty_print(self):
        maxlen = 0
        for eblr in self.ebls.values():
            if eblr.len() > maxlen:
                maxlen = eblr.len()

        print "%6s" % "ebl",
        for i in range(1, maxlen + 1):
            print "%6d" % i,
        print ""

        ebl_keys = self.ebls.keys()
        ebl_keys.sort()
        for key in ebl_keys:
            print "%6d" % key,
            self.ebls[key].pretty_print()
            print ""

    def c_best_block_widths_thresholds(self, max_block_width):
        ebl_keys = self.ebls.keys()
        ebl_keys.sort()

        bw = 1
        entries = []
        for key in ebl_keys:
            while self.ebls[key].best_block_width() > bw:
                entry = str(key)
                entries.append(entry)
                bw = bw + 1

        for i in range(0, max_block_width - len(entries)):
            entries.append(str(0))
        return "{" + ", ".join(entries) + "}"

    def max_block_width(self):
        m = 0
        for eblr in self.ebls.values():
            if eblr.max_block_width() > m:
                m = eblr.max_block_width()
        return m

class results:
    def __init__(self):
        self.mbls = {}

    def add(self, modulus_bitlen, ebls):
        self.mbls[modulus_bitlen] = ebls

    def pretty_print(self):
        mbl_keys = self.mbls.keys()
        mbl_keys.sort()
        
        for key in mbl_keys:
            print "mbl: %6d" % key 
            self.mbls[key].pretty_print()
            print "\n"

    def c_modulus_bitlen_array(self, name):
        mbl_keys = self.mbls.keys()
        mbl_keys.sort()
        return  "int %s[%d] = {" % (name, len(mbl_keys)) \
            + ", ".join(map(tostr, mbl_keys)) + "};"

    def max_block_width(self):
        m = 0
        for mblr in self.mbls.values():
            if mblr.max_block_width() > m:
                m = mblr.max_block_width()
        return m

    def c_best_block_widths_thresholds(self, name):
        mbl_keys = self.mbls.keys()
        mbl_keys.sort()
        m = self.max_block_width()
        rows = []
        for key in mbl_keys:
            rows.append(self.mbls[key].c_best_block_widths_thresholds(m))
        return "int %s[%d][%d] = \n{\n" % (name, len(mbl_keys), m) \
            + ",\n".join(rows) + "\n};"



def time_simexp(no_bases, modulus_bitlen, exponents_bitlen, block_width):
    return float(os.popen("%s -ps %s %s %s %s" %
                          (simexp_command,
                           no_bases,
                           modulus_bitlen,
                           exponents_bitlen,
                           block_width)).read().strip())

def time_exponent_bitlen(modulus_bitlen, exponents_bitlen, max_block_width):
    no_bases = test_size / (modulus_bitlen * exponents_bitlen)
    eblr = exponent_bitlen_results()
    for bw in range(1, max_block_width):
        eblr.add(time_simexp(no_bases,
                             modulus_bitlen,
                             exponents_bitlen,
                             bw))
    return eblr

def time_modulus_bitlen(modulus_bitlen,
                        exponents_bitlens,
                        max_block_width):
    mblr = modulus_bitlen_results()
    for ebl in exponents_bitlens:
        mblr.add(ebl, time_exponent_bitlen(modulus_bitlen,
                                           ebl,
                                           max_block_width))
    return mblr

def time(modulus_bitlens,
         exponents_bitlens,
         max_block_width):
    res = results()
    for mbl in modulus_bitlens:

        print "mbl: %6d" % mbl
        sys.stdout.flush()

        entry = time_modulus_bitlen(mbl,
                                    exponents_bitlens,
                                    max_block_width)
        entry.pretty_print()
        print ""
        sys.stdout.flush()

        res.add(mbl, entry)
    return res



test_size = 2000000000

modulus_bitlens = []
for i in range(0, 4):
    modulus_bitlens.append(256*2**i)

exponents_bitlens = []
for i in range(100, 2 * max(modulus_bitlens), 50):
    exponents_bitlens.append(i)

results = time(modulus_bitlens, exponents_bitlens, 12)

#results.pretty_print()

print results.c_modulus_bitlen_array("modulus_bitlens")

print results.c_best_block_widths_thresholds("best_block_widths")
