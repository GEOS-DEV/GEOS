# FILE ....... parseRstNarrative.py

# PURPOSE .... This script parses a .cpp file and creates an .rst file in a
#              narrative format. Code comments in the .cpp become standard text,
#              while actual code becomes a verbatim block in the rst. The target
#              use case is to transform a unit or integrated test into a quick user
#              tutorial by adding longer comments in the implementation.

# CREATED .... Sep 2018

import sys, string

# ... DETECT RST BLOCK START


def detectStart(rawLine):
    startDetected = ("BEGIN_RST_NARRATIVE" in rawLine)
    if (startDetected):
        splitLine = rawLine.split()    # split using whitespace as separator
        outputFileName = splitLine[-1]    # -1 = last string element
    else:
        outputFileName = ''
    return (startDetected, outputFileName)


# ... DETECT RST BLOCK FINISH


def detectFinish(rawLine):
    finishDetected = ("END_RST_NARRATIVE" in rawLine)
    return finishDetected


# ... MAIN PARSER


def main():

    # confirm number onf input arguments
    assert len(sys.argv) == 2, "Script was expecting one input argument"

    # open desired input file
    inputFileName = str(sys.argv[1])
    print "Attempting to parse file:", inputFileName
    inputFile = open(inputFileName, "r")

    # define state variables
    activeRstBlock = False
    activeVerbatimBlock = False

    # loop over file line by line and apply parsing logic
    for rawLine in inputFile:

        # if not activeRstBlock, look for a start flag.  once one is found, open
        # a new RST file with the desired output filename.  we also save an
        # uncommented version of the code in rawCode so we can print a clean
        # version at the end to allow users to easily see the overall code flow.
        if (not activeRstBlock):
            (startDetected, outputFileName) = detectStart(rawLine)
            if (startDetected):
                activeRstBlock = True
                print "Creating new RST file:", outputFileName
                outputFile = open(outputFileName, "w")
                codeCount = 0
                rawCode = ""

        # we are in an activeRstBlock state.  either we will encounter a finish flag
        # or a line we want to write to output
        else:
            finishDetected = detectFinish(rawLine)

            # if we find an end flag, we first print an uncommented version of the
            # code, close the file, and switch back to an inactiveRstBlock state
            if (finishDetected):
                activeRstBlock = False
                outputFile.write("\n**Clean Code:**\n\n")
                outputFile.write(".. code:: cpp\n\n")
                outputFile.write(rawCode)
                outputFile.close()

            # if we get here, this means we have a line we want to write into the
            # output RST file.  just need to decide if it is code or comment.
            else:
                splitLine = rawLine.strip().split()
                emptyLine = (len(splitLine) == 0)
                commentLine = (not emptyLine and "//" in splitLine[0])

                # handle empty lines
                if (emptyLine):
                    outputFile.write("\n")
                    if (activeVerbatimBlock):
                        rawCode += "\n"

                # handle comment lines, stripping off //.
                # a comment line also means that we can't be in an active
                # verbatim code block
                elif (commentLine):
                    if (activeVerbatimBlock):
                        outputFile.write("\n")
                        activeVerbatimBlock = False
                    cleanLine = rawLine.strip().replace("//", "").strip()
                    outputFile.write(cleanLine)
                    outputFile.write("\n")

                # handle code lines.  if this is the first code line after a
                # comment block, we need to start the verbatim block. also, rst
                # code blocks do not preserve leading white space, and therefore
                # the code indentation between blocks cannot be preserved.  as a
                # quick workaround, we begin the code line with a line number.  we
                # should explore other strategies down the road.
                else:
                    if (not activeVerbatimBlock):
                        outputFile.write("\n.. code:: cpp \n\n")
                        activeVerbatimBlock = True
                    outputFile.write("\t" + str(codeCount).zfill(3) + ")  " +
                                     rawLine)
                    rawCode += "\t" + rawLine
                    codeCount += 1

    # close file before quitting
    inputFile.close()
    print "Completed parsing file:", inputFileName


# ... EXECUTE MAIN

main()
