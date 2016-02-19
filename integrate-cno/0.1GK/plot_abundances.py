import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
# Take the net_history input file as a command line argument
parser.add_argument('nethistoryfile', type=str)
args = parser.parse_args()

f = open(args.nethistoryfile,'r')
h = f.readline().split()
d = {}
for hi in h:
        d[hi] = []
        
for l in f:
    ls = l.split()
    for hi, lsi in zip(h, ls):
        d[hi].append(float(lsi))

for hi in h:
    d[hi] = np.array(d[hi])

for hi in h:
    if hi != 'Time':
        fig = plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(d['Time'], d[hi])
#        ax.set_yscale('log')
        ax.set_xscale('log')
        plt.xlabel('T (s)')
        his = hi.split('_')
        hip = his[0] + '_{' + his[1] + '}'
        plt.ylabel('$\\mathrm{'+hip+'}$')
        plt.savefig(hi + '.pdf')

fig = plt.figure()
ax=fig.add_subplot(111)
legend_handles = []
for hi in h:
        if 'Y_' in hi:
                his = hi.split('_')
                hip = his[0] + '_{' + his[1] + '}'
                ax.plot(d['Time'], d[hi], label='$'+hip+'$')
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('T (s)')
framepos = ax.get_position()
ax.set_position([framepos.x0, framepos.y0, framepos.width*0.8, framepos.height])
plt.legend(bbox_to_anchor=(1.4,1))
plt.ylabel('Y')
plt.savefig('net_history.pdf')
