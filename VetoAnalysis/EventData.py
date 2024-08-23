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

#pathList=pathList[:500]
print(len(pathList))

#pathList = pathList[:500]
nEventsnTracks = 0 # Number of events with  tracks
nEventsnTracks1Vertices = 0 # Number of events with >=1 vertex and  tracks
nnUVSurvive = 0 # Number of events with  tracks surviving after  UV veto
nnMTSurvive = 0 # Number of events with  tracks surviving after  MT veto
nnCOSurvive = 0 # Number of events with  tracks surviving after  CO veto
nnFWSurvive = 0 # Number of events with  tracks surviving after  FW veto
nnIPSurvive = 0 # Number of events with  tracks surviving after  IP veto
nnTNSurvive = 0 # Number of events with  tracks surviving after  TN veto
# Organize from least effective to most effective:
# 1. UpVertex
# 2. Material
# 3. Cutoff
# 4. FW Veto 
# 5. IP Veto 
# 6. Track Number

nEvents = 0 # Total Number of events
MatVet = VT.MaterialVeto(int(layers))
materials = MatVet.MakeMaterials()
for f in pathList:
    events=joblib.load(f)
    file_hits = events["hits"]
    file_tracks = events["tracks"]
    file_vertices = events["vertices"]
    nEvents += len(file_hits)
    hitlayers = []
    for i in range(len(file_hits)):
        hits= file_hits[i]
        tracks = file_tracks[i]
        vertices = file_vertices[i]
        nEventsnTracks += 1

        for hit in hits:
            if hit.layer not in hitlayers:
                hitlayers.append(hitlayers)
                if len(hitlayers) >= 3:
                    break

        if len(vertices) > 0:
            nEventsnTracks1Vertices += 1
        else:
            continue
        if not VT.UpVertexVeto(vertices, tracks):
            nnUVSurvive += 1
            if not VT.CutoffVeto(vertices, int(layers)):
                nnCOSurvive += 1
                if not VT.FW_Veto.Chi2Veto(vertices, tracks, hits, 50):
                    nnFWSurvive += 1
                    if not VT.IPVeto.DistanceVeto(vertices, tracks, hits, 100):
                        nnIPSurvive += 1
                        if not VT.TrackNumberVeto(vertices, 3):
                            nnTNSurvive += 1

                            
            

print("Number of events:", nEvents)
print("Number of events with  tracks:", nEventsnTracks)
print("Number of events with  tracks and at least one vertex:", nEventsnTracks1Vertices)
print("Number of events with  tracks and at least one vertex and triggered:", nEventsnTracks1Vertices)
print("Number of events with  tracks and passing UV veto:", nnUVSurvive)
print("Number of events with  tracks and passing MT veto:", nnMTSurvive)
print("Number of events with  tracks and passing CO veto:", nnCOSurvive)
print("Number of events with  tracks and passing FW veto:", nnFWSurvive)
print("Number of events with  tracks and passing IP veto:", nnIPSurvive)
print("Number of events with  tracks and passing TN veto:", nnTNSurvive)

