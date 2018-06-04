import numpy as np
from scipy import signal
import os
from readAndCropJPG import readAndCropJPG
import matplotlib.pyplot as plt
from makeHumanData import makeHumanData

# choose what you want to plot ; can be ['circles', ...]
for condition in ['circles', 'gestalts', 'hexagons', 'irreg1', 'irreg2', 'malania', 'octagons', 'pattern2', 'patternIrregular', 'patternStars', 'squares', 'stars']:

    # choose correct vernier for normalization of the model results
    vernierType = 'small'
    templateSize = 40 # not useful anymore
    if condition is 'malania':
        vernierType = 'malania'
        templateSize = 120
    if condition is 'squares'  or 'pattern' in condition:
        vernierType = 'squares'
        templateSize = 40
    if condition is 'gestalts':
        vernierType = 'gestalt'
        templateSize = 60
    vernierType = vernierType.lower()
    vernierType = vernierType[0].upper() + vernierType[1:]
    vernierType = 'vernier'+vernierType

    # create stimuli results list (enter the stimuli you wish to plot)
    outputs = os.listdir('results/left')

    # a dict with stimuli names and corresponding crowding results
    crowdingResults = {}

    # a dict with stimuli names and corresponding human data
    humanData = makeHumanData(condition, vernierType)

    for stim in outputs:
        if stim[-3:] == 'npy' and (condition in stim or vernierType in stim):

            left  = np.reshape(np.load('results/left/'  + stim), (12, 257, 257))
            right = np.reshape(np.load('results/right/' + stim), (12, 257, 257))

            templateSize *= 1
            (nRows, nCols) = (left.shape[1], left.shape[2])
            # (cornerRow, cornerCol) = (int((nRows-templateSize)/2), int((nCols-templateSize)/2))
            # left  = left [:, cornerRow:cornerRow+templateSize, cornerCol:cornerCol+templateSize]
            # right = right[:, cornerRow:cornerRow+templateSize, cornerCol:cornerCol+templateSize]
            # filterRow = np.atleast_2d(signal.gaussian(nRows, std=templateSize))
            # filterCol = np.atleast_2d(signal.gaussian(nCols, std=templateSize/2))
            # filterImg = np.dot(filterRow.T, filterCol)
            # left  *= filterImg
            # right *= filterImg


            thisResult = np.sum(np.abs((left - right)))/np.sum(left + right)
            # thisResult = 1.0/np.sum(np.abs(left-right))
            if condition in stim:
                crowdingResults[stim[:-4]] = thisResult
            if vernierType in stim:
                crowdingVernier = thisResult
            print stim[:-4] + ': ' + str(thisResult)

            # plt.figure()
            # plt.imshow((left-right))
            # plt.show()

    maxModel = max(max(crowdingResults.values()), crowdingVernier)
    for d in crowdingResults.keys():
        crowdingResults[d] = 1.1*maxModel-crowdingResults[d]
    humanVernier = humanData[vernierType]
    modelVernier = 1.1*maxModel-crowdingVernier
    for d in crowdingResults.keys():
        crowdingResults[d] *= humanVernier/modelVernier

    ####### PLOT RESULTS #######

    N = len(crowdingResults)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.35       # the width of the bars

    fig, ax = plt.subplots()
    # print sorted(humanData.values())

    rects1 = ax.bar(ind, [value for (key, value) in sorted(humanData.items())[:len(crowdingResults)]], width, color='b')
    rects2 = ax.bar(ind + width, [value for (key, value) in sorted(crowdingResults.items())], width, color=(1.0,1.0,0.0))

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Thresholds')
    ax.set_title('Stimuli = '+condition)
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(sorted(humanData.keys()))
    ax.plot([-0.3,len(crowdingResults)], [humanVernier, humanVernier], 'r--')
    ax.legend((rects1[0], rects2[0]), ('Human Data', 'Model'))

    # def autolabel(rects):
    #     """
    #     Attach a text label above each bar displaying its height
    #     """
    #     for rect in rects:
    #         height = rect.get_height()
    #         ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
    #                 '%d' % int(height),
    #                 ha='center', va='bottom')

    # autolabel(rects1)
    # autolabel(rects2)

    # plt.show()
    plt.savefig('results/plots/'+condition+'.png')