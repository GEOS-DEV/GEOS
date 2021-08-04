'''
Created on 9/02/2021

@author: macpro
'''

from acquisition import EQUISPACEDAcquisition
from client import Client, SLURMCluster, Future
from AcousticShot import acoustic_shots
from utils import parse_args

args = parse_args()

acq = EQUISPACEDAcquisition(boundary=[[0,13520],[0,13520],[0,4200]],
                            dt=0.002,
                            velocity_model="/beegfs/jbesset/codes/SEP_REDUCE_Model/338x338x105_velModel.geos",
                            start_source_pos    = [4001, 7001],
                            end_source_pos      = [10001, 7001],
                            start_receivers_pos = [[101, 7001]],
                            end_receivers_pos   = [[13001, 7001]],
                            number_of_sources   = 2,
                            number_of_receivers = 100,
                            source_depth = 4099,
                            receivers_depth = 4149)


acq.add_xml(args.xml)
acq.limitedAperture(3000, 2000)

acqs = acq.split(2)

cluster = SLURMCluster(job_name="seismicAcquisition", nodes=1, cores=16, 
                       env_extra=["physics/geosx_deps","physics/pygeosx"], python="$Python3_EXECUTABLE",
                       node_name="bora")


client = Client(cluster)
client.scale(2)

#future = client.submit(acoustic_shots, acq, args.xml, cores=6, x_partition=6)
futures = client.map(acoustic_shots, acqs, cores=16, x_partition=8, y_partition=2)

#result = future.result()
results = client.gather(futures)

#print(result)
print(results)
