import os
import argparse
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

# some helpers
descr = 'Post processing for Class09 Cases \n'
parser = argparse.ArgumentParser(description=descr)
parser.add_argument('-f', nargs=1, help='specify the logfilei path to load')
parser.add_argument('-t', nargs=1, help='specify drain or hyst for ref loading (drain|hyst)')
args = parser.parse_args()

def extract_from_log(prefix, name):
    cmd = f'grep Dissolved {name} | awk "' + 'BEGIN{s = 0;}{print \\$8\\" \\"\\$9\\", \\"\\$12\\" \\"\\$13;}' + f' "> {prefix}_dissolved'
    # cmd1 = f'grep Trapped {name} | awk "' + '{print \\$7\\" \\"\\$8;}' + f' "> {prefix}_trapped'
    cmd2 = f'grep \'Phase mass\' {name} | awk "' + '{print \\$11\\" \\"\\$12;}' + f' "> {prefix}_mass'
    if  not os.system(cmd2) and not os.system(cmd):
        print('success')
    else:
        raise BaseException("grep|awk not working")

    t0 = np.loadtxt(f'{prefix}_dissolved', delimiter=',')[:,2]
    # t1 = np.loadtxt(f'{prefix}_trapped', delimiter=',')[:,0]
    t2 = np.loadtxt(f'{prefix}_mass', delimiter=',')[:,0]

    return t0,t2

def extract_time_stat(prefix, name):
    cmd = f' grep -v ^Rank {name} | grep -B1 Statistics | grep Time | awk "' + '{print \\$2;}' + '" | sed ' + ' "s/s,//g" ' + f' > {prefix}_timestep'
    # print(cmd)
    if not os.system(cmd):
        print('success')
    else:
        raise BaseException("grep|awk not working")

    time = np.loadtxt(f'{prefix}_timestep')

    return time


#fetching ref values
if (args.t[0] == 'hyst'):
    ref = pd.read_csv(os.getcwd() + '/article-results/data_hyst.csv')
    ofs = [str(i) for i in [1,3,5,7,9,11]]
elif (args.t[0] == 'drain'):
    ref = pd.read_csv(os.getcwd() + '/article-results/data_drain.csv')
    ofs = [str(i) for i in [9,11,5,7,13,15]]
else:
    raise NotImplemented

if __name__ == "__main__":

	# some conversion
	YEAR_IN_SEC = 86400 * 365

	#time and values from refs
	tags = ['GEM', 'GPRS', 'SLB']
	ref_co2g_ti = {}
	ref_co2g = {}
	ref_co2l_ti = {}
	ref_co2l = {}

	ix = 0
	for i, tg in enumerate(tags):
		ref_co2g_ti[tg] = ref['gas_'+tg].to_numpy()
		ref_co2g[tg] = ref['Unnamed: '+ofs[ix]].to_numpy()
		ix = ix + 1
		ref_co2l_ti[tg] = ref['diss_'+tg].to_numpy()
		ref_co2l[tg] = ref['Unnamed: '+ofs[ix]].to_numpy()
	##

	# raise clearer exception if input files doesn t exist
	if not os.path.exists(args.f[0]):
		raise ValueError(f'Log file does not exists at {args.f[0]}')

	# get values
	sco2_l, sco2_g = extract_from_log(args.t[0], args.f[0])

	#matplotlib structs
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	dt_list = extract_time_stat(args.t[0], args.f[0])
	plt.plot(dt_list/YEAR_IN_SEC, sco2_l, '--r')
	plt.plot(dt_list/YEAR_IN_SEC, sco2_g, '-r')
	plt.plot(dt_list/YEAR_IN_SEC, (sco2_l+sco2_g), ':r')

	lgd = ['m_co2_l', 'm_co2_g', 'm_co2_tot']
	colors = [ 'b', 'g', 'm']
	for i, tg in enumerate(tags):
		plt.plot(ref_co2g_ti[tg], ref_co2g[tg], colors[i])
		plt.plot(ref_co2l_ti[tg], ref_co2l[tg], '--'+colors[i])
		plt.plot(ref_co2l_ti[tg], (ref_co2g[tg] + ref_co2l[tg]), ':'+colors[i])
		lgd.append('ref( '+tg+' )')
		lgd.append('_nolegend_')
		lgd.append('_nolegend_')

	ax.legend(lgd, loc='upper center', bbox_to_anchor=(0.5,1.15),
				ncol=2, fancybox=True, shadow=True)
	
	plt.xlabel('time[years]')
	plt.ylabel('Co_2 mass [kg]')
	plt.xlim([0, 50])
	plt.ylim([0, 14e9])
	plt.savefig(f'figure_compare_{args.t[0]}.png')
