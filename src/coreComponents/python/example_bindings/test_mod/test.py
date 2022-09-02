def printStringFromPython(testString):
    print 'Recieved string from C++: %s' % (testString)
    return testString + '_modified_by_python'


def modifyNumpyArray(inputArray):
    import sys

    print 'Recieved array from C++'
    print inputArray

    inputArray += 0.123
