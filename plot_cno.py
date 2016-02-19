import numpy as np
import matplotlib.pyplot as plt

# Temperature in units of K
t = np.linspace(0.5e9, 0.6e9, 1000, endpoint=True)
t9 = t/1.0e9

rate_n13_beta = np.ones(1000)*np.exp(-6.76)
ctemp1 = np.array([1.813560e+01, 0.000000e+00, -1.516760e+01, 9.551660e-02, 3.065900e+00, -5.073390e-01, -6.666670e-01])
ctemp2 = np.array([1.099710e+01, -6.126020e+00, 1.571220e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, -1.500000e+00])
def c2r(ct, tgk):
    s = [ct[i]*(tgk**((2.0*i-5.0)/3.0)) for i in xrange(1,6)]
    return np.exp(ct[0] + np.sum(s) + ct[6]*np.log(tgk))
r1 = [c2r(ctemp1, ti) for ti in t9]
r2 = [c2r(ctemp2, ti) for ti in t9]
rate_n13_pg = (np.array(r1) + np.array(r2))*1.4*(1e-4)*(0.7)
for i in xrange(len(rate_n13_pg)):
    if rate_n13_pg[i] == 0.0:
        print i
        print rate_n13_pg[i]
tau_pg = 1/rate_n13_pg
tau_beta = np.log(2.0)/rate_n13_beta

# Plot
fig = plt.figure()
ax = fig.add_subplot(111)
#plt.xticks([])
#plt.yticks([])
ax.set_xlim([0.5, 0.6])
bh, = ax.plot(t9, tau_beta, 'b', label='Beta Decay')
ph, = ax.plot(t9, tau_pg, 'r', label='Proton Capture')
#ax.set_yscale('log')
#fig.set_size_inches(5,3)

plt.xlabel('T (GK)')
plt.ylabel('$\\tau ~ \mathrm{(s)}$')
#ax.text(1.025, -0.05, '$E/E_F$', transform=ax.transAxes)
#ax.text(-0.01, 1.05, '$F$', transform=ax.transAxes)
plt.legend(handles=[bh, ph], loc='upper right')
#plt.title(which_den+' density URCA', y = 1.01)

plt.savefig('cno_temperature.pdf')    
