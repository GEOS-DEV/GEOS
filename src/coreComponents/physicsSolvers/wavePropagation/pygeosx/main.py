'''
Created on 9/02/2021

@author: macpro
'''

from acquisition import EQUISPACEDAcquisition
from client import *
from AcousticShot import *
from utils import parse_args

args = parse_args()

acq = EQUISPACEDAcquisition(boundary=[[0,2000],[0,2000],[0,2000]],
                            dt=0.05,
                            start_source_pos    = [751, 1001],
                            end_source_pos      = [1501, 1001],
                            start_receivers_pos = [[21, 1001]],
                            end_receivers_pos   = [[1981, 1001]],
                            number_of_sources   = 3,
                            number_of_receivers = 10,
                            source_depth = 101,
                            receivers_depth = 51)


cluster = SLURMCluster(job_name="seismicAcquisition", nodes=2, cores=8, 
                       env_extra=["physics/geosx_deps","physics/pygeosx"], python="$Python3_EXECUTABLE",
                       node_name="bora")


client = Client(cluster)
client.scale(2)

acqs = acq.split(3)

futures = client.map(acoustic_shots, acqs, args.xml, cores=6, x_partition=6)

result = client.gather(futures)

print(result)
