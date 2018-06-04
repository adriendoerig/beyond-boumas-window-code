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

    vStim = [a for a in outputs if (a[-3:] == 'npy' and vernierType in a and 'vernier' in a)][0]
    vLeft  = np.reshape(np.load('results/left/'  + vStim), (12, 257, 257))
    vRight = np.reshape(np.load('results/right/' + vStim), (12, 257, 257))

    for stim in outputs:

        if stim[-3:] == 'npy' and (condition in stim or vernierType in stim):

            left  = np.reshape(np.load('results/left/'  + stim), (12, 257, 257))
            right = np.reshape(np.load('results/right/' + stim), (12, 257, 257))

            thisResult = np.sum(left*vLeft + right*vRight)/(np.sum(left)+np.sum(right))
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
    rects1 = ax.bar(ind, [value for (key, value) in sorted(humanData.items())[:len(crowdingResults)]], width, color='b')
    rects2 = ax.bar(ind + width, [value for (key, value) in sorted(crowdingResults.items())], width, color=(1.0,1.0,0.0))

    # add some text for labels, title and axes ticks, and save figure
    ax.set_ylabel('Thresholds')
    ax.set_title('Stimuli = '+condition)
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(sorted(humanData.keys()))
    ax.plot([-0.3,len(crowdingResults)], [humanVernier, humanVernier], 'r--')
    ax.legend((rects1[0], rects2[0]), ('Human Data', 'Model'))
    plt.savefig('results/plots/'+condition+'.png')