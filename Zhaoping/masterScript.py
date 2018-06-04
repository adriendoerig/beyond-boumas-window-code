import numpy
import matplotlib.pyplot as plt
from createInput import createInput
from createNetwork import *

# Parameters of the simulation
inputName   = 'zhaoping 1HBig'
plotFun     = False
plotFilters = False
slotSize    = 22   # 11
filterSize  = 8    # 8
nOri        = 12   # 12
sigmaX      = 1.0  # 1.0
sigmaY      = 40.0 # 40.0
oLambda     = 40.0 # 40.0
phi         = 0
maxRange    = 10 # max connection distance (manhattan-like)
gain        = 2
dt          = 0.1
nTimeSteps  = int(15/dt)
integrationWindowSize = int(5/dt)

# Create an input for the Zhaoping network from the original image
input, nRows, nCols = createInput(inputName, slotSize, filterSize, nOri, sigmaX, sigmaY, oLambda, phi, gain, plotFun, plotFilters)

# Create the network and its connections
alphaX, alphaY, Tx, Ly, g1, g2, Ic, I0, J0 = createParams()
J, W = createConnections(nOri, maxRange)
x, y = initializeNetwork(nRows, nCols, nOri, I0, Ic)

# Run the network, then produces figures
thisStimulusOutput = generateOutput(x, y, input, I0, Ic, alphaX, alphaY, Tx, Ly, g1, g2, nOri, nRows, nCols, J0, J, W, dt, nTimeSteps, integrationWindowSize, filterSize, sigmaX, sigmaY, oLambda, phi, slotSize, inputName)
