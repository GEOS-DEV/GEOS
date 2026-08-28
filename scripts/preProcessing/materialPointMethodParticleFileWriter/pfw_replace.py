import numpy as np                   
import os                                  
import importlib
import argparse
import re
import math
import sys

# The purpose of this script is to provide a simple parametric method of generating input scripts

class ParameterSet:
    def __init__(self, pairings, name_format=None):
        self.pairings = pairings
        self.numPairs = len(list(self.pairings.values())[0])
        self.name_format = name_format
        if self.name_format is None:
            self.name_format = [ "" for i in range(self.numPairs)]
        
        # All pairings must have same number of values
        for key, value in pairings.items():
            assert len(value) == self.numPairs, "All value arrays must have the same length!"
    
    def size(self):
        return self.numPairs


class Combination:
    def __init__(self, parameterSets):
        self.parameterSets = parameterSets
        self.numParamterSets = len(self.parameterSets)

        self.numCombinations = 1
        self.setSizes = [None for d in range(self.numParamterSets)]
        self.m = [None for dd in range(self.numParamterSets)]
        for d in range(self.numParamterSets):
            self.setSizes[d] = self.parameterSets[d].size()
            self.numCombinations *= self.setSizes[d]
            self.m[d] = self.numCombinations

        # print(self.numParamterSets, self.numCombinations, self.setSizes)

    def getParameterIndices(self, c):
        pi = [None for d in range(self.numParamterSets)]
        for d in range(self.numParamterSets):
            if d == 0:
                pi[d] = c % self.m[d]
            else:
                pi[d] = int(math.floor((c % self.m[d]) / self.m[d - 1]))
        return pi
    
    def generateCombinations(self):
        combinations = {}
        for s in self.parameterSets:
            for key, value in s.pairings.items():
                combinations[key] = [None for n in range(self.numCombinations)]

        # print(combinations)
        
        names = ["" for c in range(self.numCombinations)]
        for c in range(self.numCombinations):
            indices = self.getParameterIndices(c)
            for i, s in enumerate(self.parameterSets):
                names[c] += s.name_format[indices[i]]
                if i < self.numParamterSets-1:
                    names[c] += "_"
                for key, value in s.pairings.items():
                    combinations[key][c] = value[indices[i]]

        return combinations, names


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Manipulate PFW input files')
    parser.add_argument('file', help="file with job variation specs")
    parser.add_argument('-o', '--output_dir', default=None, help="directory to output variants")
    args = parser.parse_args()
    print("Command line args",args)

    assert os.path.isfile(args.file), "Could not find file, please check the path!"
    file_dir = os.path.dirname(args.file)

    sys.path.append(file_dir)
    inputFile = __import__(re.sub(".py", "", os.path.basename(args.file)))
    

    template_name = os.path.join(file_dir, inputFile.template_name)
    combObj = inputFile.combObj

    combinations, names = combObj.generateCombinations()
    instance_base_name = re.sub("_template.py", "", os.path.basename(template_name))

    if args.output_dir == None:
        output_dir = file_dir
    else:
        assert os.path.isdir(args.output_dir), "Could not find output dir!"
        output_dir = args.output_dir

    # Read template and replace parameters
    lines = ""
    with open(template_name,'r') as f_in:
        lines = f_in.readlines()
        f_in.close()

    for i in range(combObj.numCombinations):
        instance_name = os.path.join(output_dir, instance_base_name + "_" + names[i] + ".py")
        instance_txt = ""

        print("Generating instance:", instance_name)

        # Replace parameters from template script
        for line in lines:
            L = line
            for key, values in combinations.items():
                pattern  =  r"([\s\t]*" + key.replace('"', '\"').replace('[','\[').replace(']','\]') + "[\s\t]*=[\s\t]*)(.*)"
                result = re.search(pattern,L)
                if result is not None:
                    value = values[i]
                    L = result.group(1)
                    L += str(value)
                    # if type(value) is float or type(value) is int:
                    #     L += str(value)
                    # else:
                    #     L += '"' + value + '"'

                    print(L)
                    L += '\n'
                    break
                # print(result, pattern, L, result.group(2) if result is not None else "not found")

            instance_txt += L

        # Write instance text to file
        with open(instance_name, 'w') as f_out:
            f_out.write(instance_txt)
            f_out.close()
        
        print()

    print("Job names for runClean submission")
    print("\n".join([re.sub('pfw_input_', "",instance_base_name + "_" + n) for n in names]))

