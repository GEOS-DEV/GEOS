#!/usr/bin/env python

import os, sys, glob, getopt
from roofline import roofline


def process(fin, caliperKernelName, LABELS, FLOPS, AIL1, AIL2, AIHBM):

    lTime = []
    lFlop = []
    lL1 = []
    lL2 = []
    lDram = []

    nsComputeFile = fin + ".nscompute"
    res = open(nsComputeFile, "r")
    prevLine = ""

    for line in res:

        ### Time
        if "smsp__cycles_elapsed.sum.per_second" in line:
            linesp = line.split(",")
            tmpRate = float(linesp[len(linesp) - 1].strip("\n").strip('"'))
            linesp = prevLine.split(",")
            tmpTotal = float(linesp[len(linesp) - 1].strip("\n").strip('"'))
            # print( "tmpRate = ", tmpRate)
            # print( "tmpTotal = ", tmpTotal)
            lTime.append(tmpTotal / tmpRate)

        ### DP FLOP
        # add + mul
        if "smsp__sass_thread_inst_executed_op_dadd_pred_on.sum" in line:
            linesp = line.split(",")
            lFlop.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')))
        if "smsp__sass_thread_inst_executed_op_dmul_pred_on.sum" in line:
            linesp = line.split(",")
            lFlop.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')))
        # fma
        if "smsp__sass_thread_inst_executed_op_dfma_pred_on.sum" in line:
            linesp = line.split(",")
            lFlop.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')) * 2.0)

        ### L1 transactions
        # global
        if "l1tex__t_sectors_pipe_lsu_mem_global_op_ld.sum" in line:
            linesp = line.split(",")
            lL1.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')))
        if "l1tex__t_sectors_pipe_lsu_mem_global_op_st.sum" in line:
            linesp = line.split(",")
            lL1.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')))
        # local
        if "l1tex__t_sectors_pipe_lsu_mem_local_op_ld.sum" in line:
            linesp = line.split(",")
            lL1.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')))
        if "l1tex__t_sectors_pipe_lsu_mem_local_op_st.sum" in line:
            linesp = line.split(",")
            lL1.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')))
        # shared
        if "l1tex__data_pipe_lsu_wavefronts_mem_shared_op_ld.sum" in line:
            linesp = line.split(",")
            lL1.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')))
        if "l1tex__data_pipe_lsu_wavefronts_mem_shared_op_st.sum" in line:
            linesp = line.split(",")
            lL1.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')))
        # atomic
        if "l1tex__t_set_accesses_pipe_lsu_mem_global_op_atom.sum" in line:
            linesp = line.split(",")
            lL1.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')))
        if "l1tex__t_set_accesses_pipe_lsu_mem_global_op_red.sum" in line:
            linesp = line.split(",")
            lL1.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')))
        if "l1tex__t_set_accesses_pipe_tex_mem_surface_op_atom.sum" in line:
            linesp = line.split(",")
            lL1.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')))
        if "l1tex__t_set_accesses_pipe_tex_mem_surface_op_red.sum" in line:
            linesp = line.split(",")
            lL1.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')))

        ### L2 transactions
        # read and write
        if ("lts__t_sectors_op_read.sum" in line) or (
            "lts__t_sectors_op_write.sum" in line
        ):
            linesp = line.split(",")
            lL2.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')))
        # atomic
        if ("lts__t_sectors_op_red.sum" in line) or (
            "lts__t_sectors_op_atom.sum" in line
        ):
            linesp = line.split(",")
            lL2.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')) * 2.0)

        ### DRAM transactions
        if ("dram__sectors_write.sum" in line) or ("dram__sectors_read.sum" in line):
            linesp = line.split(",")
            lDram.append(float(linesp[len(linesp) - 1].strip("\n").strip('"')))

        prevLine = line  # end of for

    res.close()

    caliperFile = fin + ".caliper"
    res = open(caliperFile, "r")
    prevLine = ""

    for line in res:
        if caliperKernelName in line:
            linesp = line.split()
            caliperTime = float(linesp[1])

    res.close()

    transactionSize = 32.0
    GIGA = 1.0e9

    time = min(caliperTime, sum(lTime))
    flop = sum(lFlop)
    bytesL1 = sum(lL1) * transactionSize
    bytesL2 = sum(lL2) * transactionSize
    bytesDram = sum(lDram) * transactionSize

    tmpLABEL = fin
    tmpFLOPS = flop / time / GIGA
    tmpAIL1 = flop / bytesL1
    tmpAIL2 = flop / bytesL2
    tmpAIHBM = flop / bytesDram

    print(tmpLABEL)
    print("Time\t{}".format(time))
    print("FLOP\t{}".format(flop))
    print("FLOPS\t{}".format(tmpFLOPS))
    print("AI_L1\t{}".format(tmpAIL1))
    print("AI_L2\t{}".format(tmpAIL2))
    print("AI_HBM\t{}\n".format(tmpAIHBM))

    LABELS.append(tmpLABEL)
    FLOPS.append(tmpFLOPS)
    AIL1.append(tmpAIL1)
    AIL2.append(tmpAIL2)
    AIHBM.append(tmpAIHBM)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("process.py <name_of_kernel_in_caliper> <fileroot0> <fileroot...> ")
        sys.exit(2)

    caliperKernelName = sys.argv[1]

    fileStart = 2
    numFiles = len(sys.argv) - fileStart
    print("numFiles = {}".format(numFiles))

    # Get all profiling files
    files = []
    for f in range(fileStart, numFiles + fileStart):
        files.append(sys.argv[f])
    files.sort()

    print(files)
    if not files:
        raise RuntimeError("No profiling data found")

    LABELS = []
    FLOPS = []
    AIL1 = []
    AIL2 = []
    AIHBM = []

    # Process profiling data
    for f in files:
        print("processing {}".format(f))
        process(f, caliperKernelName, LABELS, FLOPS, AIL1, AIL2, AIHBM)

    print(LABELS)
    print(FLOPS)

    # Generate roofline plot
    roofline(LABELS, FLOPS, AIL1, AIL2, AIHBM)
