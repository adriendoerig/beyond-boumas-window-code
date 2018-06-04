# Imports
import numpy

# Function 1: Take filter parameters and build 2 oriented filters with different polarities for connection pattern from the LGN to V1
# Usage : filters1, filters2 = createFilters(numOrientations=8, size=4, sigma2=0.75, Olambda=4)
def createFilters(numOrientations, size, sigmaX, sigmaY, Olambda, phi):

    # Initialize the filters
    filters1 = numpy.zeros((numOrientations, size, size))
    filters2 = numpy.zeros((numOrientations, size, size))

    # Fill them with gabors
    midSize = (size-1.)/2.
    maxValue = -1
    for k in range(0, numOrientations):
        theta = numpy.pi*(k+1)/numOrientations + phi
        for i in range(0, size):
            for j in range(0, size):
                x = (i-midSize)*numpy.cos(theta) + (j-midSize)*numpy.sin(theta)
                y = -(i-midSize)*numpy.sin(theta) + (j-midSize)*numpy.cos(theta)
                filters1[k][i][j] = numpy.exp(-((x*x)/sigmaX + (y*y)/sigmaY)) * numpy.sin(2*numpy.pi*x/Olambda)
                filters2[k][i][j] = -filters1[k][i][j]

    # #  Rescale filters so max value is 1.0
    # for k in range(0, numOrientations):
    #     maxValue = numpy.amax(numpy.abs(filters1[k]))
    #     filters1[k] /= maxValue
    #     filters2[k] /= maxValue
    #     filters1[k][numpy.abs(filters1[k]) < 0.3] = 0.0
    #     filters2[k][numpy.abs(filters2[k]) < 0.3] = 0.0

    # Tuning 8 orientation filters
    for k in range(numOrientations):
        sumXX = numpy.sum(filters1[k]*filters1[k])
        if numOrientations == 8 and k in [1,3,5,7]:
            sumXX *= 1.2
        filters1[k] /= 1.0*sumXX
        filters2[k] /= 1.0*sumXX

    import matplotlib.pyplot as plt
    for k in xrange(numOrientations):
        if k in [0,1,3]:
            plt.figure()
            plt.imshow(filters1[k], interpolation='nearest', cmap='Greys_r')
    plt.show()

    return filters1, filters2


# Function 2: Take filter parameters and build connection pooling and connection filters arrays
# Usage (for V1 e.g.) : V1PoolingFilters, V1PoolingConnections1, V1PoolingConnections2 = createPoolConnAndFilters(numOrientations=8, VPoolSize=3, sigma2=4.0,  Olambda=5, phi=0.0)
# Usage (for V2 e.g.) : V2PoolingFilters, V2PoolingConnections1, V2PoolingConnections2 = createPoolConnAndFilters(numOrientations=8, VPoolSize=7, sigma2=26.0, Olambda=9, phi=0.0)
def createPoolingConnectionsAndFilters(numOrientations, VPoolSize, sigma2, Olambda, phi):

    # Build the angles in radians, based on the taxi-cab distance
    angles = []
    for k in range(numOrientations/2):
        taxiCabDistance = 4.0*(k+1)/numOrientations
        try:
            alpha = numpy.arctan(taxiCabDistance/(2 - taxiCabDistance))
        except ZeroDivisionError:
            alpha = numpy.pi/2 # because tan(pi/2) = inf
        angles.append(alpha + numpy.pi/4)
        angles.append(alpha - numpy.pi/4)

    # This is kind of a mess, but we could code it better
    for k in range(len(angles)):
        if angles[k] <= 0.0:
            angles[k] += numpy.pi
        if numOrientations == 2: # special case ... but I could not do it otherwise
            angles[k] += numpy.pi/4

    # Sort the angles, because they are generated in a twisted way (hard to explain, but we can see that on skype)
    angles = numpy.sort(angles)

    # Set up orientation kernels for each filter
    midSize = (VPoolSize-1.0)/2.0
    VPoolingFilters = numpy.zeros((numOrientations, VPoolSize, VPoolSize))
    for k in range(0, numOrientations):
        theta = angles[k] + phi
        for i in range(0, VPoolSize):
            for j in range(0, VPoolSize):

                # Transformed coordinates: rotation by an angle of theta
                x =  (i-midSize)*numpy.cos(theta) + (j-midSize)*numpy.sin(theta)
                y = -(i-midSize)*numpy.sin(theta) + (j-midSize)*numpy.cos(theta)

                # If the rotated x value is zero, that means the pixel(i,j) is exactly at the right angle
                if numpy.abs(x) < 0.001:
                    VPoolingFilters[k][i][j] = 1.0

    # Set layer23 pooling connections (connect to points at either extreme of pooling line ; 1 = to the right ; 2 = to the left)
    VPoolingConnections1 = VPoolingFilters.copy()
    VPoolingConnections2 = VPoolingFilters.copy()

    # Do the pooling connections
    for k in range(numOrientations):

        # Want only the end points of each filter line (remove all interior points)
        for i in range(1, VPoolSize - 1):
            for j in range(1, VPoolSize - 1):
                VPoolingConnections1[k][i][j] = 0.0
                VPoolingConnections2[k][i][j] = 0.0

        # Segregates between right and left directions
        for i in range(VPoolSize):
            for j in range(VPoolSize):
                if j == (VPoolSize-1)/2:
                    VPoolingConnections1[k][0][j] = 0.0
                    VPoolingConnections2[k][VPoolSize-1][j] = 0.0
                elif j < (VPoolSize-1)/2:
                    VPoolingConnections1[k][i][j] = 0.0
                else:
                    VPoolingConnections2[k][i][j] = 0.0

    # Returns all the needed filters to the program
    return VPoolingFilters, VPoolingConnections1, VPoolingConnections2


# Set up filters for Brightness/Darkness/SurfaceSeg/BoundarySeg spreading and boundary blocking
def createBrightnessSurfaceBoundarySpreadingFilters(numOrientations, numBrightnessFlows, brightnessSpreadingSpeeds, surfaceSegmentationSpreadingSpeeds, boundarySegmentationSpreadingSpeeds):

    # Set up filters for Brightness/Darkness filling-in stage (spreads in various directions) and boundary blocking
    brightnessFlowFilter = []
    H = numOrientations - 1  # Vertical orientation index
    V = numOrientations / 2 - 1  # Horizontal orientation index
    notStopFlowOrientation = [V, H]
    if numBrightnessFlows == 8:
        UpRight = 1
        DownRight = 5
        notStopFlowOrientation = [V, H, UpRight, DownRight]

    brightnessBoundaryBlockFilter = []
    for k in range(len(brightnessSpreadingSpeeds)):
        if numBrightnessFlows == 4:      # Right, Down, Left, Up
            brightnessFlowFilter.append([[brightnessSpreadingSpeeds[k], 0], [0, brightnessSpreadingSpeeds[k]],
                                         [-brightnessSpreadingSpeeds[k], 0],
                                         [0, -brightnessSpreadingSpeeds[k]]])
        if numBrightnessFlows == 8:      # Up, Right, Up, Left, UpRight, UpLeft, DownRight, DownLeft
            brightnessFlowFilter.append([[brightnessSpreadingSpeeds[k], 0], [0, brightnessSpreadingSpeeds[k]],
                                         [-brightnessSpreadingSpeeds[k], 0], [0, -brightnessSpreadingSpeeds[k]],
                                         [brightnessSpreadingSpeeds[k], -brightnessSpreadingSpeeds[k]],
                                         [-brightnessSpreadingSpeeds[k], -brightnessSpreadingSpeeds[k]],
                                         [brightnessSpreadingSpeeds[k], brightnessSpreadingSpeeds[k]],
                                         [-brightnessSpreadingSpeeds[k], brightnessSpreadingSpeeds[k]]])
        brightnessBoundaryBlockFilter.append([])
        directionBBF1 = []
        directionBBF2 = []
        directionBBF3 = []
        directionBBF4 = []

        # One direction
        for d in range(1, (
            brightnessSpreadingSpeeds[k]+1)):  # First index indicates the only orientation that does NOT block flow
            directionBBF1.append([notStopFlowOrientation[0], -d,     0])     # Up
            directionBBF1.append([notStopFlowOrientation[0], -d,    -1])     # Up
            directionBBF2.append([notStopFlowOrientation[1], -1,    -d])     # Right
            directionBBF2.append([notStopFlowOrientation[1],  0,    -d])     # Right
            directionBBF3.append([notStopFlowOrientation[0], (d-1),  0])     # Up
            directionBBF3.append([notStopFlowOrientation[0], (d-1), -1])     # Up
            directionBBF4.append([notStopFlowOrientation[1], -1,    (d-1)])  # Left
            directionBBF4.append([notStopFlowOrientation[1],  0,    (d-1)])  # Left
        brightnessBoundaryBlockFilter[k].append(directionBBF1)
        brightnessBoundaryBlockFilter[k].append(directionBBF2)
        brightnessBoundaryBlockFilter[k].append(directionBBF3)
        brightnessBoundaryBlockFilter[k].append(directionBBF4)

        if numBrightnessFlows == 8:  # Includes main diagonals
            directionBBF5 = []
            directionBBF6 = []
            directionBBF7 = []
            directionBBF8 = []
            for d in range(1,(brightnessSpreadingSpeeds[k]+1)):                  # First index indicates the only orientation that does NOT block flow
                directionBBF5.append([notStopFlowOrientation[2], -d,    -d])     # UpRight
                directionBBF6.append([notStopFlowOrientation[3], (d-1), -d])     # UpLeft
                directionBBF7.append([notStopFlowOrientation[3], -d,    (d-1)])  # DownRight
                directionBBF8.append([notStopFlowOrientation[2], (d-1), (d-1)])  # DownLeft
            brightnessBoundaryBlockFilter[k].append(directionBBF5)
            brightnessBoundaryBlockFilter[k].append(directionBBF6)
            brightnessBoundaryBlockFilter[k].append(directionBBF7)
            brightnessBoundaryBlockFilter[k].append(directionBBF8)

    # Set up filters for Surface segmentation spreading stage (spreads in various directions) and boundary blocking
    surfaceSegFlowFilter          = []  # Down, Left, Up, Right
    surfaceSegBoundaryBlockFilter = []
    for k in range(len(surfaceSegmentationSpreadingSpeeds)):
        surfaceSegFlowFilter.append(
            [[surfaceSegmentationSpreadingSpeeds[k], 0], [0, surfaceSegmentationSpreadingSpeeds[k]]])  # Right, Down
        surfaceSegBoundaryBlockFilter.append([])
        directionBBF1 = []
        directionBBF2 = []
        # One direction
        for d in range(1,(surfaceSegmentationSpreadingSpeeds[k]+1)):   # First index indicates the only orientation that does NOT block flow
            directionBBF1.append([notStopFlowOrientation[1], -d, 0])   # Right
            directionBBF1.append([notStopFlowOrientation[1], -d, -1])  # Right
            directionBBF2.append([notStopFlowOrientation[0], -1, -d])  # Down
            directionBBF2.append([notStopFlowOrientation[0], 0, -d])   # Down
        surfaceSegBoundaryBlockFilter[k].append(directionBBF1)
        surfaceSegBoundaryBlockFilter[k].append(directionBBF2)

    # Set up filters for Boundary segmentation spreading stage (spreads in various directions) and boundary blocking
    boundarySegFlowFilter = []  # Down, Left, Up, Right
    for k in range(len(boundarySegmentationSpreadingSpeeds)):
        boundarySegFlowFilter.append(
            [[boundarySegmentationSpreadingSpeeds[k], 0], [0, boundarySegmentationSpreadingSpeeds[k]]])  # Right, Down

    return brightnessFlowFilter, brightnessBoundaryBlockFilter, surfaceSegFlowFilter, surfaceSegBoundaryBlockFilter, boundarySegFlowFilter