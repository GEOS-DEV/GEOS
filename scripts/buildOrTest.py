#!/bin/sh
"exec" "python" "-u" "-B" "$0" "$@"
#python scripts/config-build.py -hc host-configs/darwin-clang.cmake
#python scripts/config-build.py -hc host-configs/darwin-clang37.cmake
#python scripts/config-build.py -hc host-configs/darwin-gcc.cmake
#python scripts/config-build.py -bp build-xcode -hc host-configs/darwin-clang.cmake -x

darwinHosts = ["clang", "clang37", "gcc"]
chaosHosts = ["gcc@4.9.3"]

platforms = {"darwin": darwinHosts, "chaos_5_x86_64_ib": chaosHosts}
trueNames = {
    "darwin": "darwin",
    "osx": "darwin",
    "lc": "chaos_5_x86_64_ib",
    "chaos_5_x86_64_ib": "chaos_5_x86_64_ib",
    "chaos": "chaos_5_x86_64_ib",
    "linux2": "chaos_5_x86_64_ib"
}
import os
import subprocess
import tempfile
import sys


def executeSubProcess(command, workingDirectory=os.getcwd(), verbose=2, stdin=sys.stdin):

    if verbose == -1:
        verbose = globalVerbosity
    if verbose > 1:
        print("Executing: " + command + "\n\t Working Directory: " + workingDirectory)
    #***************************************************************************************************************
    #Note: Even though python's documentation says that "shell=True" opens up a computer for malicious shell commands,
    # it is needed to allow users to fully utilize shell commands, such as cd.
    #***************************************************************************************************************

    # TODO: Followable commands aren't working in Windows right now - initially there were some pickling difficulties,
    # but now we are seeing behaviors that look like multiprocessing subprocesses are being launched in incorrect directories.
    # To be troubleshooted later.
    if False:
        with tempfile.NamedTemporaryFile() as tmpFile:
            launcher = FollowableCommand(command, workingDirectory, tmpFile, stdin)
            launcher.run(startStreaming=3.0)
            tmpFile.seek(0)
            output = tmpFile.read()
            if verbose > 1 and launcher.stopFollowing == 0:
                print(output.strip())
            process = launcher.finishedProcesses.get()
    else:
        with tempfile.TemporaryFile() as tmpFile:
            process = subprocess.Popen(command,
                                       cwd=workingDirectory,
                                       shell=(os.name != "nt"),
                                       stdout=tmpFile.fileno(),
                                       stderr=subprocess.STDOUT)
            process.wait()
            tmpFile.seek(0)
            output = tmpFile.read()
            if verbose > 1:
                print(output.strip())

    process.output = output
    if process.returncode != 0 and verbose > 1:
        print("Command '" + command + "': exited with error code " + str(process.returncode))
    return process


def execute(cmd, dryRun):
    if dryRun:
        print cmd
    else:
        process = executeSubProcess(cmd)
        if (process.returncode):
            exit(process.returncode)


def main(platform, build, local, dryRun=False):
    hosts = platforms[trueNames[platform.lower()]]
    platform = trueNames[platform]
    for host in hosts:
        if local:
            execute("cd src/thirdparty && chairajabuild", dryRun)
        if build:
            execute("scripts/config-build.py -hc host-configs/%s-%s.cmake" % (platform, host), dryRun)
        cmd = "make " if build else "make  test"
        execute("cd build-%s-%s* && %s" % (platform, host, cmd), dryRun)


usage = "scripts/buildOrTest [build|test] [<platform> [--local]] "
if __name__ == "__main__":
    local = False
    if ("scripts" not in sys.argv[0]):
        print "USAGE:" + usage
        exit(1)
    try:
        build = sys.argv[1].lower() == "build"
        test = sys.argv[1].lower() == "test"

        if (not (build or test)):
            print "USAGE:" + usage
            exit(1)
    except IndexError:
        print "USAGE:" + usage
        exit(1)
    try:
        platform = sys.argv[2].lower()
        try:
            local = sys.argv[3]
        except IndexError:
            pass
    except IndexError:
        platform = sys.platform
    build = build or not test
    main(platform, build, local)
