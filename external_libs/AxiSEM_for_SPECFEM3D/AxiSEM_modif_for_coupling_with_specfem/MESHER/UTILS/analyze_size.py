#!/usr/bin/env python
import csv
import operator


l = []

f = open('fort.735', 'r')

for line in f:
    ll = line.lstrip()
    l.append((ll.split()[0], int(ll.split()[1]), float(ll.split()[2])))

f.close()

sortedlist = sorted(l, key=operator.itemgetter(2), reverse=True)

f = open('largest_arrys.dat', 'w')

for item in sortedlist:
  print >> f, '%-30s %10d %7.2f' % item

f.close()
