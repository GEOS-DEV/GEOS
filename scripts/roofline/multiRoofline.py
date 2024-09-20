import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FormatStrFormatter


SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 24

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rcParams["figure.autolayout"] = True

def generateEvelope( peakFlop, bandwidth ):
  elbow = peakFlop/bandwidth

  ai = [ 1, elbow, 100 ]
  flops = [ peakFlop - (ai[1] - ai[0]) * ( bandwidth) , peakFlop, peakFlop ]

  return ai, flops



achitectureSpecSheetData = dict( [('V100', (7.8, 0.9)),
                                  ('A100', (8.777, 2)),
                                  ('MI250X', (24, 1.64)),
                                  ('H100-PCIe', (26, 2)),
                                  ('H100-SMX', (34, 3.35))
                                  ] )

achitectureSOLData = dict( [('V100', (6.57, 0.877)),
                            ('A100', (8.79, 2.03)),
                            ('MI250X', (24, 1.64))
                            # ('H100-PCIe', (26, 2)),
                            # ('H100-SMX', (34, 3.35))
                            ] )


fig = plt.figure()
fig.set_size_inches(8, 4)
#fig.tight_layout()
ax = fig.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Arithmetic Intensity [FLOP/Byte]')
ax.set_ylabel('GFLOP/s')
ax.set_xlim([ 1.0  , 100.0 ])
ax.set_ylim([ 1.0  , 100.0 ])
ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
ax.yaxis.set_major_formatter(FormatStrFormatter("%3.3g"))
#ax.xaxis.set_minor_formatter(FormatStrFormatter("%d"))
#ax.yaxis.set_minor_formatter(FormatStrFormatter("%3.3g"))
ax.grid(which='major', axis='both', linestyle='-', linewidth=2)
ax.grid(which='minor', axis='both', linestyle=':', linewidth=1)

roofline1, = ax.plot( *generateEvelope( *achitectureSOLData['V100']),   'k',label='V100 SOL' )
roofline2, = ax.plot( *generateEvelope( *achitectureSOLData['A100']),   'b',label='A100 SOL' )
roofline3, = ax.plot( *generateEvelope( *achitectureSOLData['MI250X']), 'r',label='MI250X SOL' )
# ax.plot( *generateEvelope( *achitectureSOLData['H100-PCIe']), 'c' )
# ax.plot( *generateEvelope( *achitectureSOLData['H100-SMX']), 'r' )
# plt.text(6, 50, 'SOL Kernel Rooflines', fontsize = 16)
# plt.arrow(10, 50, 4, -26, head_width=10.0, head_length=10.0,)


performanceData = dict( [('V100', (10.59, 4.78)),
                         ('A100', ( 24.05, 7.72)),
                         ('MI250X', (24.02, 9.126))
                        #  ('H100-PCIe', (26, 2)),
                        #  ('H100-SMX', (34, 3.35))
                         ] )

data1, = ax.plot( *performanceData['V100'],   linestyle='None', color='k', marker='o' ,label='V100 MF A(x)' )
data2, = ax.plot( *performanceData['A100'],   linestyle='None', color='b', marker='o' ,label='A100 MF A(x)' )
data3, = ax.plot( *performanceData['MI250X'], linestyle='None', color='r', marker='o' ,label='MI250X MF A(x)' )



legend1 = ax.legend( handles=[roofline1, roofline2, roofline3], loc='upper left', borderaxespad=0.1 )
ax.add_artist(legend1)
legend2 = ax.legend( handles=[data1, data2, data3], loc='lower right', borderaxespad=0.1 )
#plt.legend(bbox_to_anchor=(1.0, 0.0), loc='lower right', borderaxespad=0.1)

#fig.savefig('roofline.eps')
plt.show()
