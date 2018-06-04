# Imports
from PIL import Image
from scipy.misc import imresize
import numpy as np
import matplotlib.pyplot as plt

# Functions
def readAndCropJPG(imName, resizeFactor=1.0, onlyZerosAndOnes=False, addCropx=0, addCropy=0):

    # Read the image and transforms it into a [0 to 255] array
    im = Image.open(imName)
    ## im = im.convert('L')
    im = np.array(im)
    if len(np.shape(im)) > 2:
        im = np.mean(im,2) # mean each pixels to 8bit if 24bit
    im = np.divide(im,255.0) # array of floats between 0.0 and 1.0
    if onlyZerosAndOnes:
        im = np.round(np.multiply(im,1)) # uncomment this line if you want only black/white pixels 0.8
    white = np.max(im)

    # # Remove the upper white background (to crop the stimulus image)
    # ## while all(im[0,:]==1.0):
    # ##     im = np.delete(im,0,axis=0)
    # usefulImg = np.where(im[:,int(len(im[0,:])/2)]!=white)[0]
    # indexesToRemove = range(usefulImg[0])
    # im = np.delete(im,indexesToRemove,axis=0)
    #
    # # Remove the left and right white background
    # usefulImg = np.where(im[0,:]!=white)[0]
    # indexesToRemove = np.append(np.array(range(usefulImg[0])),np.array(range(usefulImg[-1]+1,len(im[0,:]))))
    # im = np.delete(im,indexesToRemove,axis=1)
    #
    # # Remove the bottom white background
    # usefulImg = np.where(im[:,0]!=white)[0]
    # indexesToRemove = range(usefulImg[-1]+1,len(im[:,0]))
    # im = np.delete(im,indexesToRemove,axis=0)

    # Crops the image even more, if wanted (less computations)
    if addCropx != 0:
        indexesToRemove = np.append(np.array(range(addCropx)),np.array(range(len(im[0,:])-addCropx,len(im[0,:]))))
        im = np.delete(im,indexesToRemove,axis=1) # horizontal crop (remove columns on both sides)
    if addCropy != 0:
        indexesToRemove = np.append(np.array(range(addCropy)),np.array(range(len(im[:,0])-addCropy,len(im[:,0]))))
        im = np.delete(im,indexesToRemove,axis=0) # vertical crop (remove lines on both sides)

    # Resize the image if needed
    im = imresize(im, resizeFactor, interp='nearest')

    return im, np.shape(im)[0], np.shape(im)[1]