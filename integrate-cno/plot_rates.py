import numpy as np
import matplotlib.pyplot as plt

f = open('net_rates.dat','r')
ht = f.readline().split('  ')
d = {}
h = []
for hi in ht:
        hi = hi.strip().lstrip()
        if hi != '':
            print hi
            h.append(hi)
            d[hi] = []
for l in f:
    ls = l.split()
    for hi, lsi in zip(h, ls):
        d[hi].append(float(lsi))

for hi in h:
    d[hi] = np.array(d[hi])

for hi in h:
    if hi != 'temp':
        fig = plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(d['temp']/1.0e9, d[hi])
        ax.set_yscale('log')
        ax.set_xscale('log')
        plt.xlabel('T (GK)')
        p1 = hi.find('(')
        p2 = hi.find(')')
        rx = hi[p1+1:p2]
        bc = rx.strip().split(',')
        if bc[0] == '':
                runit = 's^{-1}'
        else:
                runit = 'cm^3~ mol^{-1}~ s^{-1}'
        plt.ylabel('Rate ($\\mathrm{'+runit+'}$)')
        plt.title(hi)
        plt.savefig(hi + '.pdf')

