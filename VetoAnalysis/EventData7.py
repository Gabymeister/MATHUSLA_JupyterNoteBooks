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


#sys.path.append("/project/rrg-mdiamond/owhgabri/VetoAnalysis")
#sys.path.insert(1, "/project/rrg-mdiamond/owhgabri/MATHUSLA_JupyterNoteBooks/tracker")
#os.chdir('/project/rrg-mdiamond/owhgabri/MATHUSLA_JupyterNoteBooks/')

eff = sys.argv[1]; layers = sys.argv[2]; noise = sys.argv[3];

import Vetoes as VT

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
data_top_dir = f"/project/rrg-mdiamond/data/MATHUSLA/simulation/run-2024-07-mathusla40-full/DigiOutput"
#data_top_dir = f"/project/rrg-mdiamond/data/MATHUSLA/simulation/LLPrun-2024-08-13/DigiOutput"
pathList=[]

for rootFile, dirs, files in os.walk(data_top_dir):
    for filename in files:
        if "stat0io_MuSim-noise" + noise + "-layers" + layers + "-eff" + eff in filename:
            pathList.append(os.path.join(rootFile, filename))
print(len(pathList))

pathList = pathList[:400]
# Triggered: >=3 different layers hit
nEvents23Tracks = 0 # Number of events with 2-3 tracks
nEvents46Tracks = 0 # Number of events with 4-6 tracks
nEvents7Tracks = 0 # Number of events with 7 or more tracks
nEvents23Tracks1Vertices = 0 # Number of events with >=1 vertex and 2-3 tracks
nEvents45Tracks1Vertices = 0 # Number of events with >=1 vertex and 4-6 tracks
nEvents7Tracks1Vertices = 0 # Number of events with >=1 vertex and >= 7 tracks
nEvents23Triggered = 0 # Number of events with 2-3 tracks that was triggered
nEvents46Triggered = 0 # Number of events with 4-6 tracks that was triggered
nEvents7Triggered = 0 # Number of events with >=7 tracks that was triggered
n23FWVetoed = 0 # Number of events with 2-3 tracks vetoed by FW veto
n46FWVetoed = 0 # Number of events with 4-6 tracks vetoed by FW veto
n7FWVetoed = 0 # Number of events with >=7 tracks vetoed by FW veto
n23FWVetoedTriggered = 0 # Number of events with 2-3 tracks vetoed by FW veto
n46FWVetoedTriggered = 0 # Number of events with 4-6 tracks vetoed by FW veto
n7FWVetoedTriggered = 0 # Number of events with >=7 tracks vetoed by FW veto
n23IPVetoed = 0 # Number of events with 2-3 tracks vetoed by IP veto
n46IPVetoed = 0 # Number of events with 4-6 tracks vetoed by IP veto
n7IPVetoed = 0 # Number of events with >=7 tracks vetoed by IP veto
n23IPVetoedTriggered = 0 # Number of events with 2-3 tracks vetoed by IP veto
n46IPVetoedTriggered = 0 # Number of events with 4-6 tracks vetoed by IP veto
n7IPVetoedTriggered = 0 # Number of events with >=7 tracks vetoed by IP veto
n23UVVetoed = 0 # Number of events with 2-3 tracks vetoed by UV veto
n46UVVetoed = 0 # Number of events with 4-6 tracks vetoed by UV veto
n7UVVetoed = 0 # Number of events with >=7 tracks vetoed by UV veto
n23UVVetoedTriggered = 0 # Number of events with 2-3 tracks vetoed by UV veto
n46UVVetoedTriggered = 0 # Number of events with 4-6 tracks vetoed by UV veto
n7UVVetoedTriggered = 0 # Number of events with >=7 tracks vetoed by UV veto
n23MTVetoed = 0 # Number of events with 2-3 tracks vetoed by MT veto
n46MTVetoed = 0 # Number of events with 4-6 tracks vetoed by MT veto
n7MTVetoed = 0 # Number of events with >=7 tracks vetoed by MT veto
n23MTVetoedTriggered = 0 # Number of events with 2-3 tracks vetoed by MT veto
n46MTVetoedTriggered = 0 # Number of events with 4-6 tracks vetoed by MT veto
n7MTVetoedTriggered = 0 # Number of events with >=7 tracks vetoed by MT veto
n23COVetoed = 0 # Number of events with 2-3 tracks vetoed by CO veto
n46COVetoed = 0 # Number of events with 4-6 tracks vetoed by CO veto
n7COVetoed = 0 # Number of events with >=7 tracks vetoed by CO veto
n23COVetoedTriggered = 0 # Number of events with 2-3 tracks vetoed by CO veto
n46COVetoedTriggered = 0 # Number of events with 4-6 tracks vetoed by CO veto
n7COVetoedTriggered = 0 # Number of events with >=7 tracks vetoed by CO veto
# Organize from least effective to most effective:
# 1. UpVertex
# 2. Material
# 3. Cutoff
# 4. FW Veto 
# 5. IP Veto 
# 6. Track Number


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
        nEvents +=1
        if VT.ConeVeto.Veto(vertices, tracks, hits):
            nVetoed +=1

print("Number of events:", nEvents)
print("Survival rate:", 1 - nVetoed/nEvents)

