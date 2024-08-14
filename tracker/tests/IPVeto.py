try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterator, *args, **kwargs):
        return iterator
import os, sys, glob, warnings, glob    
import numpy as np
import scipy as sp
from scipy import constants
from pylab import *
import joblib
import importlib
from importlib import reload


sys.path.append("/project/rrg-mdiamond/owhgabri/MATHUSLA_JupyterNoteBooks/tracker")
sys.path.insert(1, "/project/rrg-mdiamond/owhgabri/MATHUSLA_JupyterNoteBooks/tracker")
os.chdir('/project/rrg-mdiamond/owhgabri/MATHUSLA_JupyterNoteBooks/')

eff = sys.argv[1]; layers = sys.argv[2]; noise = sys.argv[3]; DIST = sys.argv[4];

from tests import Vetoes as VT

sys.path.append("/project/rrg-mdiamond/owhgabri/pyTracker")
sys.path.insert(1, "/project/rrg-mdiamond/owhgabri/pyTracker")
os.chdir('/project/rrg-mdiamond/owhgabri/pyTracker')
print(os.getcwd())
print(joblib.__version__)

import scipy
import copy as cp

# ROOT
import ROOT as root

# Matplotlib
import matplotlib.pyplot as plt
from matplotlib import collections, colors, transforms

import pprint

from tracker import kalmanfilter as KF
from tracker import utilities as Util
from tracker import trackfinder as TF
from tracker import datatypes
from tracker.datatypes import *


reload(TF)
reload(Util)



#data_top_dir = f"/project/rrg-mdiamond/data/MATHUSLA/simulation/run-2024-07-mathusla40-full/DigiOutput"
data_top_dir = f"/project/rrg-mdiamond/data/MATHUSLA/simulation/LLPrun-2024-07-08/DigiOutput"
pathList=[]

for rootFile, dirs, files in os.walk(data_top_dir):
    for filename in files:
        if "stat0io_MuSim-noise" + noise + "-layers" + layers + "-eff" + eff in filename:
            pathList.append(os.path.join(rootFile, filename))


print("Running IP Veto")
print("noise: " + noise)
print("layers: " + layers)
print("efficiency: " + eff)
print("cutoff: " + DIST)


nEvents = 0
nVetoed = 0
for f in pathList:
    events=joblib.load(f)
    file_hits = events["hits"]
    file_tracks = events["tracks"]
    file_vertices = events["vertices"]
    for i in range(len(file_hits)):
        hits= file_hits[i]
        tracks = file_tracks[i]
        vertices = file_vertices[i]
        if len(vertices) == 0:
            continue
        nEvents += 1
        if VT.IPVeto.DistanceVeto(vertices, tracks, hits, float(DIST)):
            nVetoed += 1

print("Number of events:", nEvents)
print("Survival rate:", 1 - nVetoed/nEvents)

