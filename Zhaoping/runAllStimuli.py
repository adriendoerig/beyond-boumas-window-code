import numpy, os
import matplotlib.pyplot as plt
from createInputAndFilters import createInput
from createNetwork import *
from multiprocessing import Pool, cpu_count

# Parameters of the simulation
inputNames  = sorted(os.listdir('./stimuli'))
inputNames = [x for x in inputNames if x != 'Thumbs.db']
# print inputNames
swapLR = 1 # mirror image to make L stim from a R stim or vice-versa
justCheckWhatIsDone = 0 # if true, the program will display the run that were already run
testRun     = 0         # if true, the program will only run 1 stimuli, without using parallel processing
plotFun     = 1
plotFilters = 0
slotSize    = 7    # 10 # 11
filterSize  = 11   #  5 #  8
nOri        = 12   # 12
sigmaX      = 1.0  # 1.0
sigmaY      = 40.0 # 40.0
oLambda     = 40.0 # 40.0
phi         = 0
maxRange    = 10 # max connection distance (manhattan-like)
gain        = 2
dt          = 0.1
nTimeSteps  = int(15/dt)           # int(15/dt)
integrationWindowSize = int(15/dt) # int(15/dt)
resizeFactor          = 1.0

# Main loop
def runOneStim(inputName):

    # Whether a run that just checks what is done or a true run
    if justCheckWhatIsDone:
        if inputName[:-4]+'.npy' in os.listdir('results'):
            print 'hey'
            return
        return
    else:
        print "current stimulus = "+inputName[:-4]

    # Create an input for the Zhaoping network from the original image
    # inputVector, nRows, nCols = createInput(inputName, slotSize, filterSize, nOri, sigmaX, sigmaY, oLambda, phi, gain, plotFun, plotFilters, resizeFactor, swapLR)
    inputVector, nRows, nCols = createInput(inputName, slotSize, filterSize, nOri, sigmaX, sigmaY, oLambda, phi, gain, plotFun, plotFilters, resizeFactor, swapLR)
    # return

    # Create the network and its connections
    alphaX, alphaY, Tx, Ly, g1, g2, Ic, I0, J0 = createParams()
    J, W = createConnections(nOri, maxRange)
    x, y = initializeNetwork(nRows, nCols, nOri, I0, Ic)

    # Run the network, then produces figures
    generateOutput(x, y, inputVector, I0, Ic, alphaX, alphaY, Tx, Ly, g1, g2, nOri, nRows, nCols, J0, J, W, dt, nTimeSteps, integrationWindowSize, filterSize, sigmaX, sigmaY, oLambda, phi, slotSize, inputName)

if __name__ == '__main__':

    if testRun:
        runOneStim(inputNames[1])
    else:
        pool = Pool(cpu_count()-1)
        pool.map(runOneStim, inputNames)
    