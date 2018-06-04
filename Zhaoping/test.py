import sys
import numpy
from scipy import signal, spatial
import itertools
import matplotlib.pyplot as plt
from createFilters import createFilters
from PIL import Image
from readAndCropJPG import readAndCropJPG
from createInput import createInput
import os

nameList = os.listdir('../stimuli')

for i in nameList:
    print i[:-4]