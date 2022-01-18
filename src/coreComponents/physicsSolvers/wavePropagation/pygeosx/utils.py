import json
import argparse
from munch import munchify
import os
from tempfile import mkstemp, mkdtemp
import uuid


def parse_args():
    """Get parameters from python command and put them into args variable
    Return
    ------
    args : list
        Variable containing command parameters
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--xml", help = "input xml file", required = True)
    args = parser.parse_args()

    return args


def args_to_json(args):
    tempfile = str(uuid.uuid4()) + ".json"
    dict = {}
    i=0
    if isinstance(args,(list,tuple)):
        for arg in args:
            dict["arg"+str(i)] = obj_to_dict(arg)
            i+=1
    else:
        dict["arg0"] = obj_to_dict(args)

    with open(tempfile, 'w+') as f:
        json.dump(dict,f)
        f.close()

    return tempfile


def obj_to_dict(obj):
    if not hasattr(obj,"__dict__"):
        return obj
    result = {}
    for key, val in obj.__dict__.items():
        if key.startswith("_"):
            continue
        element = []
        if isinstance(val, list):
            for item in val:
                element.append(obj_to_dict(item))
        else:
            element = obj_to_dict(val)
        result[key] = element
    return result


def ricker(maxT, dt, f0):
    """ Source function
    Parameters
    ----------
    maxT : float
        The max time for simulation
    dt : float
        The time step for simulation
    f0 : float
        Intensity
    Return
    ------
    fi :
        np array containing source value at all timestep
    """

    T0 = 1.0/f0;
    fi = [0.0] * int(maxT/dt)

    for t in range(int(maxT/dt)):
        t0 = dt*t
        if t0 <= -0.9*T0 or t0 >= 2.9*T0:
            fi[t] = 0.0;
        else:
            tmp      = f0*t0-1.0
            f0tm1_2  = 2*(tmp*m.pi)*(tmp*m.pi)
            gaussian = m.exp( -f0tm1_2 )
            fi[t]    = -(t0-1)*gaussian

    return fi


def obj_to_json(obj):
    tempfile = str(uuid.uuid4()) + ".json"
    dic = obj_to_dict(obj)
    with open(tempfile, 'a') as jfile:
        json.dump(dic, jfile)

    return tempfile


def json_to_dict(jfile):
    dic = json.load(jfile)
    return dic
