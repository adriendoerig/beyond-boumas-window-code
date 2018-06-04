import numpy
import matplotlib.pyplot as plt
from readAndCropJPG import readAndCropJPG
from scipy import signal


def createInput(inputName, slotSize, filterSize, nOri, sigmaX, sigmaY, oLambda, phi, gain, plotFun=False, plotFilters=False, resizeFactor=1.0, swapLR=0):

    # Create the filters
    filters = createFilters(nOri, filterSize, sigmaX, sigmaY, oLambda, phi)
    filters = filters*filters
    filterRadius = int(filterSize/2)

    # Plot the filters
    if plotFilters:
        for k in xrange(nOri):
            plt.figure(k)
            plt.imshow(filters[k])
        plt.show()

    # Read the original image
    import scipy.ndimage.interpolation
    inputImage, imageRows, imageCols = readAndCropJPG("./stimuli/"+inputName, resizeFactor)
    # inputImage = scipy.ndimage.interpolation.zoom(inputImage, (0.5,0.5), order=3)
    if swapLR:
        inputImage = numpy.fliplr(inputImage)
    inputToReturn = numpy.zeros((imageRows, imageCols, nOri))

    # Loop to create the input image for the Zhaoping network

    # MUST ADD ZERO PADDING TO THE INPUT IMAGE HERE

    fun = numpy.zeros((imageRows*filterSize, imageCols*filterSize))
    for i in xrange(filterRadius, imageRows-filterRadius):
        for j in xrange(filterRadius, imageCols-filterRadius):

            # Find max response for each slot
            maxMatch = 0.0
            maxTheta = 0
            slotImage = inputImage[i-filterRadius:i-filterRadius+filterSize, j-filterRadius:j-filterRadius+filterSize]
            if slotImage[slotImage>0].any():
                for k in xrange(nOri):
                    match = numpy.sum(numpy.multiply(filters[k], slotImage))
                    if match > maxMatch:
                        maxMatch = match
                        maxTheta = k

            # Give constant weight to max orientation (you can remove this if you want smooth input)
            useThreshold = 1
            if useThreshold:
                threshold = 3000.0 # 2500.0
                if maxMatch < threshold:
                    maxMatch = 0.0
                # else:
                #     print maxMatch

            # Create filtered input (max for iso-orientation, less and less for other orientations)
            for k in xrange(nOri):
                weight = numpy.exp(-numpy.abs(maxTheta-k)/(numpy.pi/8))
                inputToReturn[i,j,k] = weight*maxMatch

            # Visualisation of the input edges (like in Zhaoping 2003)
            fun[i*filterSize:(i+1)*filterSize, j*filterSize:(j+1)*filterSize] = maxMatch*filters[maxTheta]

    # Plot the input edges
    if plotFun:
        plt.figure(figsize=(int(imageCols*filterSize/100.0),int(imageRows*filterSize/100.0)))
        plt.imshow(fun)
        # plt.show()
        plt.savefig('results/'+inputName[:-4]+'Input.jpg')

    # Return the input array
    print 'Input created...'
    return gain*numpy.transpose(numpy.reshape(inputToReturn, (imageRows*imageCols, nOri))), imageRows, imageCols


def createFilters(numOrientations, size, sigmaX, sigmaY, Olambda, phi):

    # Initialize the filters
    filters = numpy.zeros((numOrientations, size, size))

    # Fill them with gabors
    midSize = (size-1.)/2.
    for k in range(0, numOrientations):
        theta = numpy.pi*(k+1)/numOrientations + phi
        for i in range(0, size):
            for j in range(0, size):
                x = (i-midSize)*numpy.cos(theta) + (j-midSize)*numpy.sin(theta)
                y = -(i-midSize)*numpy.sin(theta) + (j-midSize)*numpy.cos(theta)
                filters[k][i][j] = numpy.exp(-((x*x)/sigmaX + (y*y)/sigmaY)) * numpy.sin(2*numpy.pi*x/Olambda)

    # Tuning 8 orientation filters
    for k in range(numOrientations):
        sumXX = numpy.sum(filters[k]*filters[k])
        if k in [5, 11]: # vertical / horizontal orientations
            sumXX *= 1.1
        if k in [2, 8]:  # diagonal orientations
            sumXX *= 0.67
        else:            # intermediate orientations
            sumXX *= 0.7
        filters[k] /= 1.0*sumXX

    # Return the filters
    return filters
