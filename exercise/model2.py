#!/usr/bin/python
# -*- coding: latin-1 -*-
import numpy
import numpy.ma
import matplotlib.pyplot
import math
import pylab
from mpl_toolkits.mplot3d import Axes3D

# Partial plotting
minFreq = 0.0
maxFreq = 3.0
minkpar = 1.0e-4
minkperp = 1.0e-4
maxkpar = 1.0
maxkperp = 10.0

# Sampling
kpar_start = 1.0e-4
kpar_end = 1.0
kpar_samples = 40.0
kperp_start = 1.0e-4
kperp_end = 10.0
kperp_samples = 40.0
# Next one: do not start from 0.0, WHAMP does not like that.
w_start = 0.001
w_end = 3.0
w_range = numpy.arange(w_start, w_end, 0.4)

# Parameters
ionTemperature = 0.001 # keV
electronTemperature = 2.0e-3 # keV
electronGyrofreq = 200.0 # kHz
density = 2.0e8 # m^-3

# Graphics
DPI=300

# Constants
c = 299792458.0
kb = 1.3806505e-23
mp = 1.67262171e-27
me = 9.1093826e-31
q = 1.60217653e-19
mu0 = 4*math.pi*1.0e-7
epsilon0 = 8.85418781762e-12
gamma = 5.0/3.0
eV = 1.60217733e-19 # J

v_the = math.sqrt(2.0*kb * electronTemperature / me)
r_Larmore = v_the / (electronGyrofreq * 1000.0 * 2.0 * math.pi)

print ("Preparing for WHAMP...")
## Generate input file
fileWHAMPinput = open("WHAMPinput", "w")
fileWHAMPinput.write("%e %e 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n" % (density, density)) # Number density of species
fileWHAMPinput.write("%e %e 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n" % (electronTemperature, ionTemperature)) # Temperatures in keV
fileWHAMPinput.write("1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n")
fileWHAMPinput.write("1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n")
fileWHAMPinput.write("0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n")
fileWHAMPinput.write("0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n") # Species (1 H, 0 electron)
fileWHAMPinput.write("0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n")
fileWHAMPinput.write("%e\n" % (electronGyrofreq)) # electron gyrofrequency in kHz
fileWHAMPinput.write("0\n")
fileWHAMPinput.close()

fileWHAMP_CLIinput = open("WHAMP_CLIinput", "w")
fileWHAMP_CLIinput.write('f=0.01,z=0.01,p=0.01\n')#M(2)=0.0,
fileWHAMP_CLIinput.write('zpf\n')

## p perpendicular
## z parallel
#for kpar in numpy.arange(kpar_start, kpar_end, (kpar_end-kpar_start)/kpar_samples):
for kpar in numpy.logspace(numpy.log10(kpar_start), numpy.log10(kpar_end), num=kpar_samples):
   #for kperp in numpy.arange(kperp_start, kperp_end, (kperp_end-kperp_start)/kperp_samples):
   for kperp in numpy.logspace(numpy.log10(kperp_start), numpy.log10(kperp_end), num=kperp_samples):
      for f in w_range:
         fileWHAMP_CLIinput.write('f=')
         fileWHAMP_CLIinput.write(str(f))
         fileWHAMP_CLIinput.write('z=')
         fileWHAMP_CLIinput.write(str(kpar))
         fileWHAMP_CLIinput.write('p=')
         fileWHAMP_CLIinput.write(str(kperp))
         fileWHAMP_CLIinput.write('\n')

fileWHAMP_CLIinput.write('s\n\n')
fileWHAMP_CLIinput.close()

from subprocess import call
call(["sync"])

fileWHAMPOut = open('WHAMP_CLI_output.txt', 'w')
fileWHAMP_CLIinput = open("WHAMP_CLIinput", "r")

print ("Running WHAMP...")

from subprocess import Popen
proc = Popen(["./whamp" , "-file", "WHAMPinput"], stdin=fileWHAMP_CLIinput, stdout=fileWHAMPOut)
proc.wait()

print ("WHAMP done.")

print ("Loading WHAMP data...")
# old style is just load
WHAMP=numpy.fromregex('WHAMP_CLI_output.txt', r' +([0-9]+\.[0-9]+) +([0-9]+.[0-9]+) +(-?[0-9]+\.[0-9]+E[\+\-][0-9]+) +(-?[0-9]+\.[0-9]+E[\+\-][0-9]+)', numpy.dtype('f8'))
print ("Loading done.")

print ("Processing WHAMP data...")
mask1D = numpy.ma.masked_outside(WHAMP[:,2], w_start, w_end).mask
mask2Dfreq = numpy.repeat(mask1D, 4).reshape(WHAMP.shape)
WHAMPfreqrange = numpy.ma.masked_array(WHAMP, mask2Dfreq)

# Now we're also done with WHAMP

fig = matplotlib.pyplot.figure(num=None, facecolor='w', edgecolor='k')
fig.set_size_inches(5.0, 4.0)
fig.set_dpi(DPI)

axes = fig.add_subplot(111, projection='3d')

# Plot all in black 
pl3 = axes.scatter(numpy.log10(WHAMPfreqrange[::-1,1]), numpy.log10(WHAMPfreqrange[::-1,0]), WHAMPfreqrange[::-1,2], edgecolor=None, s=2, c=WHAMPfreqrange[::-1,2], cmap='jet', vmin=minFreq, vmax=maxFreq); # marker='d', cmap=matplotlib.pyplot.cm.get_cmap('Blues')
cbar = fig.colorbar(pl3, shrink=0.8)
cbar.ax.set_xlabel('$\omega/\omega_\mathrm{ce}$')
cbar.ax.tick_params(size=4)
#cbar.set_ticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])



matplotlib.rcParams.update({'font.size': 8})
matplotlib.rcParams.update({'lines.linewidth': 1})
matplotlib.rcParams.update({'axes.linewidth': 1})
matplotlib.rcParams.update({'xtick.major.size': 8})
matplotlib.rcParams.update({'ytick.major.size': 8})
matplotlib.rcParams.update({'lines.markeredgewidth': 2})

axes.set_xlabel('\n$k_\perp\cdot r_\mathrm{Le}$', fontsize=16, linespacing=3)
axes.set_ylabel('\n$k_\parallel\cdot r_\mathrm{Le}$', fontsize=16, linespacing=3)
axes.set_zlabel('\n$\omega/\omega_\mathrm{ce}$', fontsize=16, linespacing=2)
axes.set_xmargin(0.0)
axes.set_ymargin(0.0)
axes.set_zmargin(0.0)
axes.set_xlim(numpy.log10(minkperp), numpy.log10(maxkperp))
axes.set_ylim(numpy.log10(minkpar), numpy.log10(maxkpar))

axes.tick_params(axis='x', direction='out')
axes.tick_params(axis='y', direction='out')
axes.tick_params(axis='z', direction='out')

axes.w_xaxis.gridlines.set_linestyle('-')
axes.w_yaxis.gridlines.set_linestyle('-')
axes.w_zaxis.gridlines.set_linestyle('-')

axes.w_xaxis.gridlines.set_lw(0.5)
axes.w_yaxis.gridlines.set_lw(0.5)
axes.w_zaxis.gridlines.set_lw(0.5)
axes.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
axes.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
axes.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))


axes.azim = -120
axes.elev = 10

print ("Et voilà !")
matplotlib.pyplot.show()

#matplotlib.pyplot.savefig('plot2.png', dpi=DPI)
