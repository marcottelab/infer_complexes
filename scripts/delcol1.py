from __future__ import division
import sys

for line in sys.stdin:
    #python prints double-spaced basically unless you add a comma at the end
    print '\t'.join(line.split('\t')[1:]),
