# --------------------------------------------------
# User IO for tracker
# --------------------------------------------------
## The following functions needs to be defined:
## 1. load(filename, *args, **kwargs):
##      """
##      Load the data into a dictionary. 
##
##      INPUT:
##      ---
##      filename: str
##          the filename of the input file
##      *args, **kwargs: 
##          additional optional arguments
##
##      Return:
##      ---
##          data: dict
##              The key of the dictionary is the group of hits.
##              For example, there could be three groups, "inactive", "1" and "2"
##              Group "inactive" will not be processed, and can be used to store floor hits
##              Group "1" and "2" can be used to save hits in top layers and wall tracking layers
##          metadata: dict
##              This dictionary contains metadata. The necessary metadata is the direction of each group
##              Because the tracker takes "y" axis as the default layer direction, it is needed to rotate the coordinates
##              if the layer is actually facing another direction.
##              So, you should provide the coordinate rotation for each group
##              For example: metadata["groups"]["1"] = {"flip_index":[1,2]}
##      """


import math, os, sys, types
import importlib
from collections import namedtuple


import ROOT
# import joblib
import numpy as np

# sys.path.append("..")
import tracker.datatypes as datatypes

# Test:
# import importlib
# io_user = importlib.machinery.SourceFileLoader("*", "io_mu-simulation.py").load_module()
# io_user.load("/project/rrg-mdiamond/tomren/mudata/LLP_geometry_40m///SimOutput/MATHUSLA_LLPfiles_HXX/deweighted_LLP_bb_35/20240410/173948/reconstruction_v2.ROOT")  

# class datatypes:
#     Hit = namedtuple("Hit", ["x", "y", "z", "t", "x_err", "y_err", "z_err", "t_err","layer", "ind"])
#     Track = namedtuple("Track", ["x0", "y0", "z0", "t0", "Ax", "Ay", "Az", "At", "cov", "chi2", "ind", "hits", "hits_filtered"])
#     Vertex = namedtuple("Vertex", ["x0", "y0", "z0", "t0", "cov", "chi2", "tracks"])
#     VertexSeed = namedtuple("VertexSeed",["x0", "y0", "z0", "t0",  "cov", "chi2","dist", "Ntracks", "trackind1","trackind2", "score"])



det_width        = 3.5 # 4.5cm per bar
det_height       = 1 # [cm]
time_resolution  = 1 # [ns], single channel
refraction_index = 1.58

unc_trans = det_width/math.sqrt(12)                  
unc_long = time_resolution*29.979/math.sqrt(2)/refraction_index
unc_time = time_resolution/math.sqrt(2) # ns
unc_thick = det_height/math.sqrt(12) # uncertainty in thickness, cm



def load(filename, printn=2000, start_event=0, end_event=-1, *args, **kwargs):
    # Open the file
    tfile = ROOT.TFile.Open(filename)
    tree_name = tfile.GetListOfKeys()[0].GetName()
    Tree = tfile.Get(tree_name)
    branches = [Tree.GetListOfBranches()[i].GetName() for i in range(len(Tree.GetListOfBranches()))]
    entries = Tree.GetEntries()
    print("\nOpening ROOT file", filename)
    print("  Tree:", tree_name)
    # print("  Branches:", branches)
    print("  Entries:", entries)

    if end_event==-1: end_event = entries
    entries_run = [start_event, min(end_event, entries)]
    print("  Entries to read:", entries_run)
    
    # Tell the version of the data
    if 'Digi_bar_direction' in branches:
        version = 1
    else:
        version = 0
        
    # version = 0
        
        
    

    # Make a dictionary
    data = {
        "inactive":[],
        "1":[],
        "2":[],
    }

    # Read the file
    for entry in range(*entries_run):
        if (entry+1)%printn==0:  print("  event is ", entry+1)

        for key in data:
            data[key].append([])

        hits = root_hits_extractor(Tree, entry, version)
        
        # Group the hits into layer groups
        for hit in hits:
            
            # if  hit.layer<=1 or hit.layer>5:
            #     data["inactive"][-1].append(hit)            
            # elif hit.layer<=5:
            #     data["1"][-1].append(hit)
            
            if  hit.layer<=3: # Two wall and two floor layers
                data["inactive"][-1].append(hit)            
            elif hit.layer<=7: # Top 4 celling tracking layers
                
                data["1"][-1].append(hit)  
            elif hit.layer<=11: #  4 wall tracking layers
                data["2"][-1].append(hit)                  

    print("Finished loading file")

    # Make metadata
    metadata={"groups":{}}
    metadata["groups"]["1"] = {"flip_index":None}
    metadata["groups"]["2"] = {"flip_index":[1,2]}

    # size = getsize(data)
    # print(f"Memory usage of loaded data {size/1e6:.2f} [MB]")

    return data, metadata


def dump(data, filename):
    
    with open(output_filename,"wb") as f:
        pickle.dump(results, f)    
    # joblib.dump(data,filename+".joblib")

def exists(filename):
    return os.path.exists(filename+".joblib")
    
    


# ----------------------------------------------------------
# Helper functions
# ----------------------------------------------------------

def open(filename):
    tfile = ROOT.TFile.Open(filename)
    return tfile

def get_tree(tfile):
    Tree = tfile.Get(tfile.GetListOfKeys()[0].GetName())
    def keys(Tree):
        branch_list = [Tree.GetListOfBranches()[i].GetName() for i in range(len(Tree.GetListOfBranches()))]
        return branch_list
    Tree.keys = types.MethodType(keys,Tree)
    # print(Tree.keys())
    # Tree.keys = branch_list
    
    return Tree    
    
def c2list(x):
    return [x.at(i) for i in range(x.size())]



def make_hits(x, y, z, t, ylayers):
    hits=[]
    for i, layer in enumerate(ylayers):
        if layer%2==1:
            hits.append(datatypes.Hit(x[i], y[i], z[i], t[i], unc_trans, unc_thick, unc_long, unc_time, layer, i, 0, 0))
        else:
            hits.append(datatypes.Hit(x[i], y[i], z[i], t[i], unc_long, unc_thick, unc_trans, unc_time, layer, i, 0, 0))         
            
    return hits  

def make_hits_newsim(x, y, z, t, layer_id, layer_direction, bar_direction, det_id, digi_type):
    hits=[]
    # layer_direction_mod ={0:1, 1:2, 2:0}
    
    for i, layer in enumerate(layer_id):
        unc = [0,0,0]
        other_direction = 0+1+2-layer_direction[i]-bar_direction[i]
        unc[layer_direction[i]] = unc_thick
        unc[bar_direction[i]] = unc_long
        unc[other_direction] = unc_trans
        hits.append(datatypes.Hit(x[i], y[i], z[i], t[i], unc[0], unc[1], unc[2], unc_time, layer, i,det_id[i], digi_type))
            
    return hits 

def root_hits_extractor(Tree, entry, version):
    Tree.GetEntry(entry)
    Digi_x = c2list(Tree.Digi_x)
    Digi_y = c2list(Tree.Digi_y)
    Digi_z = c2list(Tree.Digi_z)
    Digi_t = c2list(Tree.Digi_time)
    Digi_layer = c2list(Tree.Digi_layer_id)
<<<<<<< HEAD
    Digi_type = c2list(Tree.Digi_type)
    return make_hits(Digi_x, Digi_y, Digi_z, Digi_t, Digi_layer, Digi_type)
=======
    
    if version==0:
        return make_hits(Digi_x, Digi_y, Digi_z, Digi_t, Digi_layer)
    elif version==1:
        Digi_layer_direction = c2list(Tree.Digi_layer_direction)
        Digi_bar_direction = c2list(Tree.Digi_bar_direction)
        Digi_det_id = c2list(Tree.Digi_det_id)
        Digi_type = c2list(Tree.Digi_type)
        return make_hits_newsim(Digi_x, Digi_y, Digi_z, Digi_t, Digi_layer, Digi_layer_direction, Digi_bar_direction, Digi_det_id, Digi_type)
        
        
>>>>>>> 9374c40979b31b733b78c3671a20fd2ef4cf37db


from types import ModuleType, FunctionType
from gc import get_referents
def getsize(obj):
    """sum size of object & members."""
    # Custom objects know their class.
    # Function objects seem to know way too much, including modules.
    # Exclude modules as well.

    BLACKLIST = type, ModuleType, FunctionType
    if isinstance(obj, BLACKLIST):
        raise TypeError('getsize() does not take argument of type: '+ str(type(obj)))
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size    
