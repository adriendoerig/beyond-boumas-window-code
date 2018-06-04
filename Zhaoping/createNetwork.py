import sys
import numpy
from scipy import signal, spatial
import itertools
import matplotlib.pyplot as plt
from createInputAndFilters import createFilters


def createParams():

    alphaX = 1
    alphaY = 1
    Tx = 1
    Ly = 1.2
    g1 = 0.21
    g2 = 2.5
    Ic = 1.0
    I0 = 0.85
    J0 = 0.8

    return alphaX, alphaY, Tx, Ly, g1, g2, Ic, I0, J0


def gX(x, Tx):

    output = numpy.zeros(x.shape)
    for i in range(len(x)):
        if x[i] < Tx:
            pass
        elif Tx <= x[i] <= (Tx+1):
            output[i] = (x[i]-Tx)
        else:
            output[i] = 1.0

    return output


def gY(y, Ly, g1, g2):

    output = numpy.zeros(y.shape)
    for i in range(len(y)):
        if y[i] < 0:
            pass
        elif 0 <= y[i] <= Ly:
            output[i] = g1*y[i]
        else:
            output[i] = g1*Ly + g2*(y[i]-Ly)

    return output


def psi(k1, y, Ly, g1, g2, nOri):

    output = numpy.zeros(y[k1,:].shape)
    for k2 in xrange(nOri):
        deltaTheta =((numpy.pi*numpy.abs(k1-k2)/nOri)%(2*numpy.pi))*180/numpy.pi
        if deltaTheta == 0:
            output += 1.0*gY(y[k2,:], Ly, g1, g2)
        elif deltaTheta == 1*180/nOri:
            output += 0.8*gY(y[k2,:], Ly, g1, g2)
        elif deltaTheta == 2*180/nOri:
            output += 0.7*gY(y[k2,:], Ly, g1, g2)
        else:
            pass

    return output


def descartesDistance(A,B,nRows,nCols):

    (i1, j1) = numpy.unravel_index(A, (nRows, nCols))
    (i2, j2) = numpy.unravel_index(B, (nRows, nCols))

    return numpy.sqrt((i1-i2)**2 + (j1-j2)**2)


def descartesAngle(i1, j1, i2, j2):

    if j1 == j2:
        return numpy.sign(i1-i2)*numpy.pi/2
    else:
        return -numpy.arctan(float(i1-i2)/float(j1-j2))


def normalizationIndexes(i, D, nRows, nCols):

    # Build a filter and take border effects into account
    (i1, j1) = numpy.unravel_index(i, (nRows, nCols))
    indices = [] # numpy.zeros((2*D+1, 2*D+1))
    for ii in range(2*D+1):
        for jj in range(2*D+1):
            if numpy.sqrt((ii-D)**2 + (jj-D)**2) <= D:
                if 0 <= i1+(ii-D) < nRows and 0 <= j1+(jj-D) < nCols:
                    indices.append(numpy.ravel_multi_index((i1+(ii-D), j1+(jj-D)), (nRows, nCols)))

    return indices


def INoise(nRows, nCols, nOri):

    return 1*numpy.random.normal(loc=0.0, scale=(0.1)**2/(2*0.1), size=(nOri, nRows*nCols))


def INormalization(x, Tx, nOri, nRows, nCols):

    outputToGive = numpy.zeros(x.shape)
    for i in xrange(nRows*nCols):
        js = normalizationIndexes(i, 2, nRows, nCols)
        for k1 in xrange(nOri):
            outputCurrent = 0.0
            for k2 in xrange(nOri):
                if k1 != k2:
                    outputCurrent += numpy.sum(gX(x[k2,js], Tx))/len(js)
            outputToGive[k1,i] = -2.0*outputCurrent**2
    return outputToGive


def createConnections(nOri, maxRange):

    excitatory = numpy.zeros((nOri, nOri, 2*maxRange+1, 2*maxRange+1))
    inhibitory = numpy.zeros((nOri, nOri, 2*maxRange+1, 2*maxRange+1))
    for i in xrange(2*maxRange+1):
        for j in xrange(2*maxRange+1):

            d = numpy.sqrt((i-maxRange)**2 + (j-maxRange)**2)
            if d <= 10.0:

                gamma = descartesAngle(i, j, maxRange, maxRange) # between -pi/2 and pi/2 (I checked)
                for k1 in xrange(nOri):
                    for k2 in xrange(nOri):

                        # theta1
                        theta = (numpy.pi*(k1+1)/nOri)%(2*numpy.pi) # between 0 and pi
                        theta1 = gamma - theta
                        if theta1 <= -numpy.pi/2:
                            theta1 += numpy.pi
                        if theta1 > numpy.pi/2:
                            theta1 -= numpy.pi

                        # theta2
                        thetaPrime = (numpy.pi*(k2+1)/nOri)%(2*numpy.pi) # between 0 and pi
                        theta2 = gamma - thetaPrime
                        if theta2 <= -numpy.pi/2:
                            theta2 += numpy.pi
                        if theta2 > numpy.pi/2:
                            theta2 -= numpy.pi

                        # deltaTheta
                        delta = theta-thetaPrime
                        if delta <= -numpy.pi/2:
                            delta += numpy.pi
                        if delta > numpy.pi/2:
                            delta -= numpy.pi

                        # according to Zhaoping advice, make sure that theta1 is smaller than theta2
                        if theta1 > theta2:
                            temp = theta1
                            theta1 = theta2
                            theta2 = temp
                            delta = -delta

                        # beta
                        beta = 2*numpy.abs(theta1) + 2*numpy.sin(numpy.abs(theta1+theta2))

                        # excitatory connections
                        if d!=0 and ((beta<numpy.pi/2.69) or (beta<numpy.pi/1.1 and numpy.abs(theta1)<numpy.pi/5.9 and numpy.abs(theta2)<numpy.pi/5.9)):
                            excitatory[k1,k2,i,j] = 0.126*numpy.exp(-(beta/d)**2-2*(beta/d)**7-d**2/90)

                        # inhibitory connections
                        if d==0 or beta<numpy.pi/1.1 or numpy.abs(delta)>=numpy.pi/3 or numpy.abs(theta1)<numpy.pi/11.999:
                            pass
                        else:
                            inhibitory[k1,k2,i,j] = 0.14*(1-numpy.exp(-0.4*(beta/d)**1.5))*numpy.exp(-numpy.sign(delta)*(numpy.abs(delta)/(numpy.pi/4))**1.5)

    # for k in xrange(nOri):
    #     plt.figure(k)
    #     plt.imshow(excitatory[k, (k+4)%nOri, :, :] - inhibitory[k, (k+4)%nOri, :, :])
    # plt.show()

    return excitatory, inhibitory


def initializeNetwork(nRows, nCols, nOri, I0, Ic):

    x = I0 + INoise(nRows, nCols, nOri)
    y = Ic + INoise(nRows, nCols, nOri)

    return x, y


def integrationStep(x, y, input, IBackExc, IBackInh, alphaX, alphaY, Tx, Ly, g1, g2, nOri, nRows, nCols, J0, J, W, dt, currentTimeStep, integrationWindowSize):

    for k1 in xrange(nOri):

        postSynExcitActivity = numpy.zeros((nRows*nCols,))
        postSynInhibActivity = numpy.zeros((nRows*nCols,))
        for k2 in xrange(nOri):

            crossOriPreSynActivity = numpy.reshape(gX(x[k2,:], Tx), (nRows, nCols))
            postSynExcitActivity += numpy.reshape(signal.convolve2d(crossOriPreSynActivity, J[k1,k2,:,:], mode='same', boundary='wrap'), postSynExcitActivity.shape)
            postSynInhibActivity += numpy.reshape(signal.convolve2d(crossOriPreSynActivity, W[k1,k2,:,:], mode='same', boundary='wrap'), postSynInhibActivity.shape)

        # Generate a live plot of the mean activity over integration windows
        if currentTimeStep%integrationWindowSize == 0 and currentTimeStep > 0:

            plt.figure((k1)%nOri)
            plt.imshow(numpy.reshape(postSynExcitActivity, (nRows, nCols)) - numpy.reshape(postSynInhibActivity, (nRows, nCols)))

        x[k1,:] += dt*(-alphaX*x[k1,:] - psi(k1, y, Ly, g1, g2, nOri) + J0*gX(x[k1,:], Tx) + postSynExcitActivity + input[k1,:] + IBackExc[k1,:])
        y[k1,:] += dt*(-alphaY*y[k1,:] + gX(x[k1,:], Tx) + postSynInhibActivity + IBackInh[k1,:])

    return x, y


def generateOutput(x, y, input, I0, Ic, alphaX, alphaY, Tx, Ly, g1, g2, nOri, nRows, nCols, J0, J, W, dt, nTimeSteps, integrationWindowSize, filterSize, sigmaX, sigmaY, oLambda, phi, slotSize, name):

    # Create the filters and slot parameters
    filters = createFilters(nOri, filterSize, sigmaX, sigmaY, oLambda, phi)
    filters = filters*filters
    nFunRows = nRows*filterSize
    nFunCols = nCols*filterSize

    # plt.figure()
    # plt.ion()
    output = numpy.zeros(x.shape)
    for timeStep in xrange(nTimeSteps+1):

        print str(timeStep)+' time steps on '+str(nTimeSteps)

        # Create background currents with noise
        IBackExc = I0 + INormalization(x, Tx, nOri, nRows, nCols) + INoise(nRows, nCols, nOri)
        IBackInh = Ic + INoise(nRows, nCols, nOri)

        # Simulate the network
        x, y = integrationStep(x, y, input, IBackExc, IBackInh, alphaX, alphaY, Tx, Ly, g1, g2, nOri, nRows, nCols, J0, J, W, dt, timeStep, integrationWindowSize)

        # Generate a live plot of the mean activity over integration windows
        for k in xrange(nOri):
            output[k,:] += gX(x[k,:], Tx)/integrationWindowSize
        if timeStep%integrationWindowSize == 0 and timeStep > 0:

            fun = numpy.zeros((nFunRows, nFunCols))
            for i in xrange(nRows*nCols):

                maxMatch = numpy.max(output[:,i])
                maxTheta = numpy.argmax(output[:,i])
                (iR, iC) = numpy.unravel_index(i, (nRows, nCols))
                fun[iR*filterSize:(iR+1)*filterSize, iC*filterSize:(iC+1)*filterSize] = maxMatch*filters[maxTheta]

            # Save output as a .npy file BEFORE it is reinitialized to zero for the next run
            numpy.save('results/'+name[:-4],output)
            output = numpy.zeros(x.shape)

            # Save fun array as a .jpg file, for visual
            plt.figure(1000, figsize=(int(nCols*filterSize/100.0),int(nRows*filterSize/100.0)))
            plt.imshow(fun)
            plt.savefig('results/'+name[:-4]+'Figure.jpg')

            # plt.pause(0.001)
            # plt.show()