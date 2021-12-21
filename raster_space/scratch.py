import os

def CreateScratchFolder(workspace):

    # check if the current path maps to geodatabase
    n = len(workspace)
    if n > 4:
        end = workspace[n-4: n] # extract last 4 letters
        if end == ".gdb":  # geodatabase
            workspace = os.path.dirname(workspace)

    i = 1
    scratch = str()
    while i != 0:
        if os.path.isdir(workspace + '/scratch' + str(i)) == False:
            os.mkdir(workspace + '/scratch' + str(i))
            scratch = workspace + '/scratch' + str(i)
            i = 0

        else:
            i += 1

    return scratch