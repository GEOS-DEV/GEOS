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
                            start_source_pos    = [7001, 7001],
                            end_source_pos      = [12001, 7001],
                            start_receivers_pos = [[21, 7001]],
                            end_receivers_pos   = [[13501, 7001]],
                            number_of_sources   = 4,
                            number_of_receivers = 676,
                            source_depth = 4099,
                            receivers_depth = 4149)


acq.add_xml(args.xml)
acq.limitedAperture(3000)

acqs = acq.split(4)
cluster = SLURMCluster(job_name="seismicAcquisition", nodes=1, cores=16, 
                       env_extra=["physics/geosx_deps","physics/pygeosx"], python="$Python3_EXECUTABLE",
                       node_name="bora", walltime="03:00:00")


client = Client(cluster)
client.scale(4)
"""
future = client.submit(acoustic_shots, acq, cores=16, x_partition=8, y_partition=2)
result = future.result()
print(result)
"""
futures = client.map(acoustic_shots, acqs, cores=16, x_partition=8, y_partition=2)
results = client.gather(futures)
print(results)

