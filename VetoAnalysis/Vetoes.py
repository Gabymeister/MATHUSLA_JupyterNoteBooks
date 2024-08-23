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
import copy as cp

from scipy.spatial import ConvexHull
import itertools

sys.path.append("/project/rrg-mdiamond/owhgabri/pyTracker/")
sys.path.insert(1, "/project/rrg-mdiamond/owhgabri/pyTracker/")
os.chdir('/project/rrg-mdiamond/owhgabri/pyTracker/')
print(os.getcwd())
print(joblib.__version__)

# ROOT
import ROOT as root

import pprint

from tracker import kalmanfilter as KF
from tracker import utilities as Util
from tracker import trackfinder as TF
from tracker import vertexfinder as VF
from tracker import datatypes
from tracker.datatypes import *

reload(TF)
reload(VF)
reload(Util)


# ---Constants------------------------------------------------------------------ 
c = 29.979
wallMid1 = 6896.6
wallMid2 = 6998.2
floorMid1 = 8550.8
floorMid2 = 8632.4
# -----6 layer--------------------
y_layer_starts = [9730, 9811.6, 9893.2, 9974.8, 10056.4, 10138] # 6 layer 
z_wall_starts = [10900, 10981.6, 11063.2, 11144.8, 11226.4, 11308] # 6 layer

x_mod_starts = [-1950, -950, 50, 950, 1050]
n_floors = 2
z_mod_starts = [7000, 8000, 9000, 10000]
beam_dim = 10
mod_length = 900
front_wall_starts = [6895.8,6997.4]
c = 29.979

wallMid1 = 6896.6
wallMid2 = 6998.2
floorMid1 = 8550.8
floorMid2 = 8632.4

IP = (0,0,0)

y_bottoms = [8550, 8631.6] # Floor Vetos
z_fronts = [6895.8, 6997.4] # Wall Vetos
x_lims = (-1950, 1950)
thickness = 1.6

# ---Utility Functions----------------------------------------------------------
def ConvertVelocity(track):
    """
    Converts the velocity from either dk/dy or dk/dz
    to dk/dt
    """
    if track.Ay == 1: # Horizontal layer track
        A_1x = track.Ax/track.At; A_1y = 1/track.At; A_1z = track.Az/track.At
    else: # Vertical layer track
        A_1x = track.Ax/track.At; A_1y = track.Ay/track.At; A_1z = 1/track.At
    return np.array((A_1x, A_1y, A_1z))


def SortByTime(points):
    """
    Sorts list hits by time
    A hit is of the form (x,y,z,t)
    """
    if len(points) <= 1:
        return points
    pivot = points[len(points) // 2]
    left = []
    middle = []
    right = []
    for point in points:
        if point[-1] < pivot[-1]:
            left.append(point)
        elif point[-1] > pivot[-1]:
            right.append(point)
        elif point[-1] == pivot[-1]:
            middle.append(point)
    return SortByTime(left) + middle + SortByTime(right)


def GoingUp(points):
    """
    Takes a sorted list of points and check if it is going up.
    If it is, return True
    """
    for i in range(1, len(points)):
        if points[i-1][1] > points[i][1]:
            return False
    return True


def ProjectionTime(track, points, layer):
    """
    Takes a track and its SORTED list of points
    closest is true if you want to check in the closest wall/floor,
    false if you wanna check in the farther wall/floor
    Returns the time at which the track is 
    projected to enter the wall or floor
    Given point is (x,y,z,t)
    """
    if layer == 0:# Walls
        if track.Ay == 1: # Horizontal layer track 
            return track.At/track.Az*(wallMid1 - points[0][2]) + points[0][3]
        else: # Vertical layer track
            return track.At*(wallMid1 - points[0][2]) + points[0][3]
    elif layer == 1:
        if track.Ay == 1: # Horizontal layer track
            return track.At/track.Az*(wallMid2 - points[0][2]) + points[0][3]
        else: # Vertical layer track
            return track.At*(wallMid2 - points[0][2]) + points[0][3]
    elif layer == 2: # Floors
        if track.Ay == 1: # Horizontal layer track
            return track.At*(floorMid1 - points[0][1]) + points[0][3]
        else: # Vertical layer track
            return track.At/track.Ay*(floorMid1 - points[0][1]) + points[0][3]
    elif layer == 3: 
        if track.Ay == 1: # Horizontal layer track
            return track.At*(floorMid2 - points[0][1]) + points[0][3]
        else: # Vertical layer track
            return track.At/track.Ay*(floorMid2 - points[0][1]) + points[0][3]

def ProjectionPoint(track, points, t):
    """
    Returns where the track is projected to be at time t.
    points is the list of hits from the track sorted by time
    Returns an empty array if out of the detector
    """
    if track.Ay == 1: #Horizontal layer track
        xf = points[0][0] + track.Ax/track.At*(t - points[0][3])
        yf = points[0][1] + 1/track.At * (t - points[0][3])
        zf = points[0][2] + track.Az/track.At*(t - points[0][3])
    else: # Vertical layer track 
        xf = points[0][0] + track.Ax/track.At*(t - points[0][3])
        yf = points[0][1] + track.Ay/track.At * (t - points[0][3])
        zf = points[0][2] + 1/track.At*(t - points[0][3])
    return np.array((xf,yf,zf))

def GetDistance(point1, point2):
    """
    Get the distance between two points (x,y,z)
    """
    return np.sqrt((point1[0]-point2[0])**2 + (point1[1]-point2[1])**2 + (point1[2]-point2[2])**2)

#----------------------------------------------------------------------------------------
class FW_Veto:
    
    def ProjectionCovariance(track, points, t):
            """
            Returns the projected covariance of x,y,z at time t.
            covariance is x,z,t,Ax,Az,At if horizontal track, or 
            x,y,t,Ax,Ay,At if vertical track.
            points is sorted by time
            """
            dt = t - points[0][3]
            if track.Ay == 1: # Horizontal track
                Jac = np.array([[1, 0, -track.Ax/track.At, dt/track.At, 0, -track.Ax*dt/track.At**2],
                                [0, 0, -1/track.At, 0, 0, -dt/track.At**2],
                                [0, 1, -track.Az/track.At, 0, dt/track.At, -track.Az*dt/track.At**2]])
            else: 
                Jac = np.array([[1, 0, -track.Ax/track.At, dt/track.At, 0, -track.Ax*dt/track.At**2],
                                [0, 1, -track.Ay/track.At, 0, dt/track.At, -track.Ay*dt/track.At**2],
                                [0, 0, -1/track.At, 0, 0, -dt/track.At**2]])    
            return Jac@track.cov@Jac.T
    
    def GetChiSquared(hit, track):
        """
        Get the chi squared between a track projection and the hit
        """
        TrackPoints = SortByTime(track.hits_filtered)
        hit_pos = np.array((hit.x, hit.y, hit.z))
        hit_err = np.array((hit.x_err, hit.y_err, hit.z_err))
        proj_time = hit.t
        proj_point = ProjectionPoint(track,TrackPoints, proj_time)
        proj_cov = FW_Veto.ProjectionCovariance(track, TrackPoints, proj_time)
        proj_cov += np.diag(hit_err)
        residual = proj_point - hit_pos 
        try:
            inv_cov = np.linalg.inv(proj_cov)
        except:
            return None
        return residual.T @ inv_cov @ residual

    def MinChi2(vertices, tracks, hits):
        """
        Gets the minimum chi2 of the event
        """
        min_chi2 = None
        wf_hits = []
        for hit in hits:
            if hit.layer < 4:
                wf_hits.append(hit)
        if len(wf_hits) == 0: # Event passes
            return min_chi2
        for vertex in vertices:
            v_tracks = []
            for t_index in vertex.tracks:
                v_tracks.append(tracks[t_index])
            for track in v_tracks:
                for hit in wf_hits:
                    cur_chi2 = FW_Veto.GetChiSquared(hit,track)
                    if cur_chi2 is not None and (min_chi2 is None or cur_chi2 < min_chi2):
                        min_chi2 = cur_chi2
        return min_chi2
    
    def MinDistance(vertices, tracks, hits):
        """
        Gets the minimum distance of the event between
        all vertex-track projections and f/w hits
        """
        min_dist = None
        wf_hits = []
        for hit in hits:
            if hit.layer < 4:
                wf_hits.append(hit)
        if len(wf_hits) == 0: # Event passes 
            return min_dist
        for vertex in vertices:
            v_tracks = []
            for t_index in vertex.tracks:
                v_tracks.append(tracks[t_index])
            for track in v_tracks:
                for hit in wf_hits:
                    hit_pos = np.array((hit.x, hit.y, hit.z))
                    TrackPoints = SortByTime(track.hits_filtered)
                    proj_time = hit.t
                    proj_point = ProjectionPoint(track,TrackPoints, proj_time)
                    dist = GetDistance(hit_pos, proj_point)
                    if (min_dist is None or dist < min_dist):
                        min_dist = dist
        return min_dist

    def Chi2Veto(vertices, tracks, hits, cutoff):
        """
        Decide if the event is vetoed based on the chi2 cutoff
        Return True if vetoed. Return false if not vetoed
        """
        min_chi2 = FW_Veto.MinChi2(vertices, tracks, hits)
        if min_chi2 is None:
            return False
        if min_chi2 > 0 and min_chi2 < cutoff:
            return True
        return False
 
    def DistanceVeto(vertices, tracks, hits, cutoff):
        """
        Decide if the event is vetoed based on the distance cutoff
        Return True if vetoed. Return false if not vetoed
        """
        min_dist = FW_Veto.MinDistance(vertices, tracks, hits)
        if min_dist is None:
            return False
        if min_dist < cutoff:
            return True
        return False
    
#----------------------------------------------------------------------------------------
class IPVeto:

    def IP_track_cov(hit):
        """
        Construct the covariance of an "IP track" connecting from IP to a given hit
        Two points:
        IP: [0, 0, 0, t-dt], Hit: [x, y, z, t] 
        Returns the covariance of the tracklet
        Covariance will be of the form (x,y,z,Ax,Ay,Az)
        Refer to journal/report for explanation
        """
        x,y,z,t = (hit.x, hit.y, hit.z, hit.t)
        err = (hit.x_err, hit.y_err, hit.z_err, hit.t_err)
        dt= np.linalg.norm((x,y,z))/c
        dt3c2 = dt**3 * c**2
        
        # tracklet Jacobian #
        # dependent on the difinition of the tracklet #
        jac = np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1],
                        [1/dt - x**2/dt3c2, -x*y/dt3c2, -x*z/dt3c2, 0],
                        [-x*y/dt3c2, 1/dt - y**2/dt3c2, -z*y/dt3c2, 0],
                        [-x*z/dt3c2, -z*y/dt3c2, 1/dt - z**2/dt3c2, 0]])
        
        track_cov = jac @ np.diag(np.array(err)**2) @ jac.T
        return track_cov   

    def MakeTracklet(hit):
        """
        Create a track from a hit in the floor/wall
        Assuming two points:
        IP: [0, 0, 0, t-dt], Hit: [x, y, z, t] 
        Returns a Track including covariance.
        No chi2, ind, hits, or hits_filtered included
        """
        x,y,z,t = (hit.x, hit.y, hit.z, hit.t)
        err = (hit.x_err, hit.y_err, hit.z_err, hit.t_err)
        cov = IPVeto.IP_track_cov(hit)
        dt = np.linalg.norm((x,y,z))/c
        Ax = x/dt; Ay = y/dt; Az = z/dt; At = 1
        # Ax = x/t; Ay = y/t; Az = z/t
        tracklet = datatypes.Track(x,y,z,t,Ax,Ay,Az,At,cov,0,0,[],[])
        return tracklet
        
    def TrackMinimumTime(track, tracklet):
        """
        Get the time at which the distance between
        track1 and track2 is minimum
        track2 is the tracklet. It's velocities are parametrized by time
        """
        if track.Ay == 1: # Horizontal track
            Ax1 = track.Ax/track.At; Ay1 = 1/track.At; Az1 = track.Az/track.At
        elif track.Az == 1: # Vertical track
            Ax1 = track.Ax/track.At; Ay1 = track.Ay/track.At; Az1 = 1/track.At
        A1 = np.array((Ax1, Ay1, Az1))
        point1 = np.array((track.x0, track.y0, track.z0))
        
        Ax2 = tracklet.Ax; Ay2 = tracklet.Ay; Az2 = tracklet.Az
        A2 = np.array((Ax2, Ay2, Az2))
        point2 = np.array((tracklet.x0, tracklet.y0, tracklet.z0))
        
        A_diff = A1 - A2
        P_diff = point1 - point2
        return (-np.dot(A_diff, P_diff) + np.dot(A_diff, A1*track.t0 - A2*tracklet.t0))/np.dot(A_diff, A_diff)
                  
    def TrackDistance(track, tracklet, minTime):
        """
        Get the distance between the track and tracklet at minTime
        """
        if track.Ay == 1: # Horizontal track
            Ax1 = track.Ax/track.At; Ay1 = 1/track.At; Az1 = track.Az/track.At
        elif track.Az == 1: # Vertical track
            Ax1 = track.Ax/track.At; Ay1 = track.Ay/track.At; Az1 = 1/track.At
        track1Point = np.array((track.x0 + Ax1*(minTime - track.t0), 
                               track.y0 + Ay1*(minTime - track.t0),
                               track.z0 + Az1*(minTime - track.t0)))
        trackletPoint = np.array((tracklet.x0 + tracklet.Ax*(minTime - tracklet.t0),
                                 tracklet.y0 + tracklet.Ay*(minTime - tracklet.t0),
                                 tracklet.z0 + tracklet.Az*(minTime - tracklet.t0)))
        return np.linalg.norm(track1Point - trackletPoint)
        
    def TrackletProjectionPoint(tracklet, t):
        """
        Get the projection (x,y,z) of the tracklet at time t
        """
        return np.array((tracklet.x0 + tracklet.Ax*(t-tracklet.t0),
                        tracklet.y0 + tracklet.Ay*(t-tracklet.t0),
                        tracklet.z0 + tracklet.Az*(t-tracklet.t0)))
        
    def TrackletProjectionCovariance(tracklet, t):
        """
        Get the projected covariance of the tracklet at time t
        Covariance for tracklets is of the form (x,y,z,Ax,Ay,Az)
        where A's are already with respect to time
        """
        x0,y0,z0,t0 = (tracklet.x0, tracklet.y0, tracklet.z0, tracklet.t0)
        jac = np.array([[1, 0, 0, t-t0, 0, 0],
                        [0, 1, 0, 0, t-t0, 0],
                        [0, 0, 1, 0, 0, t-t0]])
        return jac@track.cov@jac.T
    
    def MinDistance(vertices, tracks, hits):
        """
        Calculates the minimum distance of closest approach
        between all vertex-tracks and IP tracks in the event
        """
        min_dist = None
        wf_hits = []
        for hit in hits:
            if hit.layer < 4:
                wf_hits.append(hit)
        if len(wf_hits) == 0:
            return False
        for vertex in vertices:
            v_tracks = []
            for t_index in vertex.tracks:
                v_tracks.append(tracks[t_index])
            for track in v_tracks:
                for hit in wf_hits:
                    tracklet = IPVeto.MakeTracklet(hit)
                    t_shortest = IPVeto.TrackMinimumTime(track, tracklet)
                    dist = IPVeto.TrackDistance(track, tracklet, t_shortest)
                    if min_dist is None or dist < min_dist:
                        min_dist = dist
        return min_dist

    def DistanceVeto(vertices, tracks, hits, cutoff):
        """
        Decide if the event is vetoed based on the distance cutoff
        Return True if vetoed. Return false if not vetoed
        """ 
        min_dist = IPVeto.MinDistance(vertices, tracks, hits)
        if min_dist is None:
            return False
        if min_dist < cutoff:
            return True
        return False

#----------------------------------------------------------------------------------------
class MaterialVeto:
    """
    This only vetoes if all vertices are within the cutoff range
    """
    def __init__(self, nLayers):
        """
        Set the number of layers for the material veto
        """
        self.y_layers = y_layer_starts[len(y_layer_starts) - nLayers:]
        self.z_walls = z_wall_starts[:nLayers]
        self.wall_top = y_layer_starts[nLayers - n_floors]

    def MakeLayers(self):
        """
        Returns a list of the layer modules of 
        the form
        ((x_min, x_max),(y_min,y_max),(z_min, z_max))
        """
        layers = []
        for i in range(len(self.y_layers)):
            for j in range(len(x_mod_starts)):
                for k in range(len(z_mod_starts)):
                    y_min = self.y_layers[i]; y_max = self.y_layers[i] + thickness
                    x_min = x_mod_starts[j]; x_max = x_mod_starts[j] + mod_length
                    z_min = z_mod_starts[k]; z_max = z_mod_starts[k] + mod_length
                    layers.append(np.array(((x_min,x_max),(y_min,y_max),(z_min,z_max))))
        return layers
    
    def MakeBacks(self):
        """
        Returns a list of the back modules of 
        the form
        ((x_min, x_max),(y_min,y_max),(z_min, z_max))
        """
        backs = []
        y_min = self.wall_top - mod_length; y_max = self.wall_top
        for j in range(len(x_mod_starts)):
            for k in range(len(self.z_walls)):
                x_min = x_mod_starts[j]; x_max = x_mod_starts[j] + mod_length
                z_min = self.z_walls[k]; z_max = self.z_walls[k] + thickness
                backs.append(np.array(((x_min, x_max),(y_min,y_max),(z_min,z_max))))
        return backs
    
    def MakeWalls(self):
        """
        Returns a list of the hermetic walls of the form
        ((x_min, x_max),(y_min,y_max),(z_min, z_max))
        """
        walls = []
        y_min = y_bottoms[0]; y_max = self.wall_top
        x_min = x_mod_starts[0]; x_max = x_mod_starts[-1] + mod_length
        for i in range(len(front_wall_starts)):
            z_min = front_wall_starts[i]; z_max = front_wall_starts[i] + thickness
            walls.append(np.array(((x_min,x_max),(y_min,y_max),(z_min,z_max))))
        return walls
    
    def MakeFloors(self):
        """
        Returns a list of the hermetic walls of the form
        ((x_min, x_max),(y_min,y_max),(z_min, z_max))
        """
        floors = []
        z_min = z_mod_starts[0]; z_max = z_mod_starts[-1] + mod_length
        x_min = x_mod_starts[0]; x_max = x_mod_starts[-1] + mod_length
        for i in range(len(y_bottoms)):
            y_min = y_bottoms[i]; y_max = y_bottoms[i] + thickness
            floors.append(np.array(((x_min,x_max),(y_min,y_max),(z_min,z_max))))
        return floors
        
    def MakeBeams(self):
        """
        Returns a list with beams of the form 
        ((x_min, x_max),(y_min,y_max),(z_min, z_max))
        """
        beams = []
        all_layers = y_bottoms + self.y_layers
        for i in range(1, len(all_layers) - 1):
            for j in range(len(x_mod_starts)):
                for k in range(len(z_mod_starts)):
                    y_min = all_layers[i] + thickness; y_max = all_layers[i+1]
                    
                    x_min = x_mod_starts[j]; x_max = x_mod_starts[j] + beam_dim
                    z_min = z_mod_starts[k]; z_max = z_mod_starts[k] + beam_dim
                    beams.append(np.array(((x_min, x_max),(y_min,y_max),(z_min,z_max))))
                    
                    x_min = x_mod_starts[j]; x_max = x_mod_starts[j] + beam_dim
                    z_min = z_mod_starts[k] + mod_length - beam_dim; z_max = z_mod_starts[k] + mod_length
                    beams.append(np.array(((x_min, x_max),(y_min,y_max),(z_min,z_max))))
    
                    x_min = x_mod_starts[j] + mod_length - beam_dim; x_max = x_mod_starts[j] + mod_length
                    z_min = z_mod_starts[k]; z_max = z_mod_starts[k] + beam_dim
                    beams.append(np.array(((x_min, x_max),(y_min,y_max),(z_min,z_max))))
    
                    x_min = x_mod_starts[j] + mod_length - beam_dim; x_max = x_mod_starts[j] + mod_length
                    z_min = z_mod_starts[k] + mod_length - beam_dim; z_max = z_mod_starts[k] + mod_length
                    beams.append(np.array(((x_min, x_max),(y_min,y_max),(z_min,z_max))))
        return beams

    def MakeGround(self):
        """
        Returns a list with the ground of the form
        ((x_min, x_max),(y_min,y_max),(z_min, z_max))
        """
        ground = []
        x_min = -9999999; x_max = 9999999
        y_min = -9999999; y_max = 8547
        z_min = -9999999; z_max = 9999999
        ground.append(np.array(((x_min,x_max), (y_min,y_max),(z_min,z_max))))
        return ground
    
    def MakeMaterials(self):
        """
        Returns all the materials concatenated together
        """
        beams = self.MakeBeams()
        layers = self.MakeLayers()
        backs = self.MakeBacks()
        walls = self.MakeWalls()
        floors = self.MakeFloors()
        ground = self.MakeGround()
        return np.concatenate((beams,layers,backs,walls,floors, ground))
    
    def GetResidual(self, vertex, rectangle):
        """
        Get the distance between the vertex and a rectangle
        The rectangle is of the form 
        ((x_min, x_max),(y_min,y_max),(z_min, z_max))
        """
        x_min, x_max = rectangle[0]
        y_min, y_max = rectangle[1]
        z_min, z_max = rectangle[2]
        x_v,y_v,z_v = np.array((vertex.x0, vertex.y0, vertex.z0))
        x_clam = max(x_min, min(x_v, x_max))
        y_clam = max(y_min, min(y_v, y_max))
        z_clam = max(z_min, min(z_v, z_max))
        return np.array(((x_v - x_clam), (y_v - y_clam), (z_v - z_clam)))
    
    def GetChi2(self, vertex, rectangle):
        """
        Get the chi2 between the vertex and the rectangle
        """
        res = self.GetResidual(vertex, rectangle)
        trimmed_cov = np.delete(vertex.cov, 3, 0)
        trimmed_cov = np.delete(trimmed_cov, 3, 1)
        return res.T@np.linalg.inv(trimmed_cov)@res
    
    def MinDistance(self, vertex, materials):
        """
        Returns the minimum distance to a material
        """
        min_dist = None
        for material in materials:
            res = self.GetResidual(vertex, material)
            dist = np.sqrt(res[0]**2 + res[1]**2 + res[2]**2)
            if min_dist is None or dist < min_dist:
                min_dist = dist
        return min_dist
    
    def MinChi2(self, vertex):
        """
        Returns the minimum chi2 of the vertex
        """
        materials = self.MakeMaterials()
        min_chi2 = None
        for material in materials:
            chi2 = self.GetChi2(vertex, material)
            if min_chi2 is None or chi2 < min_chi2:
                min_chi2 = chi2
        return min_chi2
    
    def MinMinDistance(self, vertices, materials):
        """
        Get the minimum of the minimum distances
        in the event (strict version)
        """
        min_dist = None 
        for vertex in vertices:
            dist = self.MinDistance(vertex, materials)
            if min_dist is None or dist < min_dist:
                min_dist = dist
        return min_dist

    def MaxMinDistance(self, vertices, materials):
        """
        Get the maximum of the minimum distances
        in the event (lenient version)
        """
        max_dist = None 
        for vertex in vertices:
            dist = self.MinDistance(vertex, materials)
            if max_dist is None or dist > max_dist:
                max_dist = dist
        return max_dist

    def DistanceVeto(self, vertices, materials, cutoff, strict=False):
        """
        Vetoes an event if all vertices are within the cutoff distance from a material.
        Return True if vetoed. Return false if not vetoed
        if strict: veto if a single vertex is within the cutoff distance from a material
        """
        max_dist = None # Max of the minimum distances
        min_dist = None #
        for vertex in vertices:
            dist = self.MinDistance(vertex, materials)
            if max_dist is None or dist > max_dist:
                max_dist = dist
            if min_dist is None or dist < min_dist:
                min_dist = dist
        if not strict and max_dist is None or max_dist > cutoff:
            return False
        elif strict and min_dist is None or min_dist > cutoff:
            return False
        return True

#----------------------------------------------------------------------------------------
def CutoffVeto(vertices, nLayers, strict=False):
    """
    Vetoes an event based on vertex being in decay volume
    Return True if vetoed. Return false if not vetoed
    if strict: veto if a single vertex is outside the cutoff
    """
    for vertex in vertices:
        if not (vertex.y0 < y_layer_starts[len(y_layer_starts) - nLayers] and \
        vertex.y0 > y_bottoms[0] + 20 and vertex.x0 > x_mod_starts[0] + 20 and \
        VErtex.x0 < x_mod_starts[-1] + mod_length - 20 and \
        vertex.z0 > front_wall_starts[0] + 20 and \
        vertex.z0 < z_wall_starts[-1] - 20): # One outside 
            if strict:
                return True
        else: # One inside
            if not strict:
                return False
    if strict: # None outside
        return False
    else: # All outside
        return True

#----------------------------------------------------------------------------------------
def TrackNumberVeto(vertices, n, strict=False):
    """
    Vetoes an event based on track number cutoff n
    Return True if vetoed. Return false if not vetoed
    if strict: veto if a single vertex has less than n tracks
    """
    for vertex in vertices:
        if len(vertex.tracks) < n:             
            if strict: # one vertex less than n tracks
                return True
        else:
            if not strict: # one vertex >= than n tracks
                return False
    if strict: # All vertices >= n tracks
        return False
    if not strict: # No vertices >= n tracks
        return True

#----------------------------------------------------------------------------------------
def UpVertexVeto(vertices, tracks):
        """
        Vetoes an event based on whether there is a vertex with
        downward going tracks. Return True if vetoed. Return false if not vetoed
        """
        for i in range(len(vertices)):
            vertex = vertices[i]
            isupward=True
            for i in vertex.tracks:
                # All tracks need to go upwards
                # All tracks need to be higher than the vertex
                if (tracks[i].Ay/tracks[i].At<0) or (tracks[i].y0 < vertex.y0):
                    return True


#----------------------------------------------------------------------------------------
class ConeVeto:

    def MaxTracks(vertex, tracks):
        """
        Get the track indices that have the largest angle of separation
        """
        inds = None
        min_angle = None
        for i in range(len(vertex.tracks)):
            ind_1 = vertex.tracks[i]
            t1 = tracks[ind_1]
            A1 = ConvertVelocity(t1)
            for j in range(i + 1, len(vertex.tracks)):
                ind_2 = vertex.tracks[j]
                t2 = tracks[ind_2]
                A2 = ConvertVelocity(t2)
                angle = np.arccos(np.dot(A1, A2)/(np.linalg.norm(A1)*np.linalg.norm(A2)))
                if min_angle is None or angle < min_angle:
                    min_angle = angle
                    inds = (ind_1, ind_2)
        return inds 
    
    def ClosestApproach(point, track):
        """
        Gets the point at which the track is closest to point
        """
        A = ConvertVelocity(track)
        P_0 = np.array((track.x0, track.y0, track.z0))
        t_c = np.dot((point - P_0), A)/np.linalg.norm(A)**2 + track.t0
        return P_0 + A*(t_c - track.t0)
    
    def EllipseAngle(p1, p2):
        """
        Get the angle of the ellipse. This is always with 
        respect to the x axis
        """
        if p1[1] - p2[1] < 0.001: # Floor
            dz = p2[2] - p1[2]; dx = p2[0] - p1[0]
            slope = dz/dx
            angle = np.arctan(dz/dx)
        else: # Wall
            dy = p2[1] - p1[1]; dx = p2[0] - p1[0]
            slope = dy/dx
            angle = np.arctan(dy/dx)        
        return angle
            
    def GetEllipse(vertex, tracks, layer):
        """
        Gets the center of the backwards cone of the vertex
        """
        
        ind_1, ind_2 = ConeVeto.MaxTracks(vertex, tracks)
        t1 = tracks[ind_1]; t2 = tracks[ind_2]
        t1Points = SortByTime(t1.hits_filtered); t2Points = SortByTime(t2.hits_filtered)
        time1 = ProjectionTime(t1, t1Points, layer); time2 = ProjectionTime(t2, t2Points, layer)
        p1 = ProjectionPoint(t1, t1Points, time1); p2 = ProjectionPoint(t2, t2Points, time2)
        midpoint = (p1 + p2)/2
        c1 = midpoint[0]
        if layer < 2: # Wall
            c2 = midpoint[1]
            a = GetDistance(p1, midpoint) # Axis 1
            b = GetDistance(ConeVeto.ClosestApproach(midpoint, t1), midpoint) # Axis 2
            angle = ConeVeto.EllipseAngle(p1, p2)
            def Ellipse(x, y):
                """
                Returns whether (x,y) lies within the ellipse
                """
                x_p = (x-c1)*np.cos(angle) + (y-c2)*np.sin(angle)
                y_p = (y-c2)*np.cos(angle) - (x-c1)*np.sin(angle)
                return ((x_p)**2/a**2 + (y_p)**2/b**2) <= 1
        else: # Floor
            c2 = midpoint[2]
            a = GetDistance(p1, midpoint) # Axis 1
            b = GetDistance(ConeVeto.ClosestApproach(midpoint, t1), midpoint) # Axis 2
            angle = ConeVeto.EllipseAngle(p1, p2)
            def Ellipse(x, z):
                """
                Returns whether (x,z) lies within the ellipse
                """
                x_p = (x-c1)*np.cos(angle) + (z-c2)*np.sin(angle)
                z_p = (z-c2)*np.cos(angle) - (x-c1)*np.sin(angle)
                return ((x_p)**2/a**2 + (z_p)**2/b**2) <= 1
        return Ellipse
            
    def Ellipses(vertex, tracks):
        """
        Gets the Ellipses for each layer
        """
        ells = []
        for i in range(len(y_bottoms) + len(z_fronts)):
            ells.append(ConeVeto.GetEllipse(vertex,tracks,i))
        return ells

    def Veto(vertices, tracks, hits):
        """
        Vetoes an event based on track number cutoff n
        Return True if vetoed. Return false if not vetoed
        """
        wf_hits = []
        for hit in hits:
            if hit.layer < 4:
                wf_hits.append(hit)
        for vertex in vertices:
            layer_ells = ConeVeto.Ellipses(vertex, tracks)
            for hit in wf_hits:
                layer = hit.layer
                if hit.layer < 2 and layer_ells[layer](hit.x, hit.y): # Wall
                    return True
                elif layer_ells[layer](hit.x,hit.z): # Floor
                    return True
        return False



class Detector:
    # Use cm as unit
    cm = 1      

    OutBoxLimits =[[-1950.0, 1950.0],  
                   [8547.0, 8547.0+1754.200],  
                   [7000.0, 7000.0+3900]]
    
    BoxMargin = 20 # Leave 20 cm on all sides
    BoxDecayHeight = 10.968*100
    BoxLimits   = [[-1950.0 + BoxMargin, 1950.0 - BoxMargin],  
                   [8547.0 +84.6 + BoxMargin, 8547.0 + BoxDecayHeight - BoxMargin],  
                   [7000.0 + BoxMargin, 7000.0+3900 - BoxMargin]    ]
    
    WallLimits  = [ BoxLimits[0], [BoxLimits[1][0],BoxLimits[1][0] + 2000.0 ] ] # [x,y] [cm,cm]

    scintillator_height_all = 1.6 # 2cm +0.3*2Al case
    # 2023-04-25 Tom: changing to 6 layers. New numbers pull from simulation


    LayerYLims= [[8547.000,8547.000+scintillator_height_all],
                 [8628.600,8628.600+scintillator_height_all],
                 [9890.200,9890.200+scintillator_height_all],
                 [9971.800,9971.800+scintillator_height_all],
                 [10053.400,10053.400+scintillator_height_all],
                 [10135.000,10135.000+scintillator_height_all]]    
        
    
    
    y_floor = LayerYLims[2][1]
    z_wall = BoxLimits[2][0] + 3 # [cm] add 3 cm to account for wall width

    ModuleXLims = [ [-4950. + 1000.*n, -4050. + 1000*n] for n in range(10) ]
    ModuleZLims = [ [7000.  + 1000.*n,  7900. + 1000*n] for n in range(10) ]


    def __init__(self):
        # print("Detector Constructed")
        self.time_resolution=1e-9 #1ns=1e-9s
        
        self.n_top_layers = 5
        self.x_edge_length = 39.0 # decrease from 99->39 -- Tom
        self.y_edge_length = 39.0 # decrease from 99->39 -- Tom
        self.n_modules = 4

        self.x_displacement = 70.0
        self.y_displacement = -49.5
        self.z_displacement = 6001.5

        self.layer_x_edge_length = 9.0
        self.layer_y_edge_length = 9.0

        self.scint_x_edge_length = 4.5
        self.scint_y_edge_length = 0.045
        self.scintillator_length = self.scint_x_edge_length
        self.scintillator_width = self.scint_y_edge_length
        self.scintillator_height = 0.02
        self.scintillator_casing_thickness = 0.003

        self.steel_height = 0.03

        self.air_gap = 30

        self.layer_spacing = 0.8
        self.layer_count   = 8

        self.module_x_edge_length = 9.0
        self.module_y_edge_length = 9.0
        self.module_case_thickness = 0.02

        self.full_layer_height =self.layer_w_case = self.scintillator_height + 2*self.scintillator_casing_thickness

        self.layer_z_displacement = [                          0.5*self.layer_w_case,
                                        1*self.layer_spacing + 1.5*self.layer_w_case,
                                    5 +                        0.5*self.layer_w_case,
                                    5 + 1*self.layer_spacing + 1.5*self.layer_w_case,
                                    5 + 2*self.layer_spacing + 2.5*self.layer_w_case,
                                    5 + 3*self.layer_spacing + 3.5*self.layer_w_case,
                                    5 + 4*self.layer_spacing + 4.5*self.layer_w_case,
                                    5 + 5*self.layer_spacing + 5.5*self.layer_w_case]

        self.module_x_displacement = [i*(self.module_x_edge_length + 1.0) -0.5 * self.x_edge_length + 0.5*self.module_x_edge_length for i in range(self.n_modules)]
        self.module_y_displacement = [i*(self.module_y_edge_length + 1.0) -0.5 * self.y_edge_length + 0.5*self.module_y_edge_length for i in range(self.n_modules)]
        pass
    def xLims(self):
        return self.BoxLimits[0]

    def yLims(self):
        return self.BoxLimits[1]

    def zLims(self):
        return self.BoxLimits[2]

    def LayerY(self, n):
        return self.LayerYLims[n]

    def LayerYMid(self, n):
        return (self.LayerYLims[n][0] + self.LayerYLims[n][1])/2.

    def numLayers(self):
        return len(self.LayerYLims)

    def DrawColor(self):
        return "tab:gray"

    def inLayer(self, yVal):
        for layerN, layerLims in enumerate(self.LayerYLims):
            if yVal > layerLims[0] and yVal < layerLims[1]:
                return layerN
        return -1

    def nextLayer(self, yVal):
        for n in range(len(self.LayerYLims)-1):
            if yVal > self.LayerYLims[n][1] and yVal < self.LayerYLims[n+1][0]:
                return n+1
        return 999


    def inLayer_w_Error(self, yVal, yErr):
        for layerN, layerLims in enumerate(self.LayerYLims):

            lower = yVal - yErr
            upper = yVal + yErr

            if lower < layerLims[0] and upper > layerLims[1]:
                return layerN

            if lower < layerLims[0] and upper > layerLims[0]:
                return layerN

            if lower < layerLims[1] and upper > layerLims[1]:
                return layerN


        return -1

    def inBox(self, x, y, z):
        if x > self.xLims()[0] and x < self.xLims()[1]:
            if y > self.yLims()[0] and y < self.yLims()[1]:
                if z > self.zLims()[0] and z < self.zLims()[1]:
                    return True
        return False
    
det=Detector()  




    
class Cut:
    
    CONST_CMS_LOC = np.array([0,0,0])
    CONST_VETO_LAYER = np.array([0,1,2,3])
    CONST_NOISE_RATE_SUPRESS = 0.1
    RNG = np.random.default_rng()
    
    @staticmethod
    def angle(v1,v2):
        """Find angle between two vectors"""
        mag1 = np.linalg.norm(v1)
        mag2 = np.linalg.norm(v2)
        opening_angle = np.arccos(np.dot(v1, v2)/( mag1*mag2 ))
        return opening_angle
    


    @staticmethod
    def find_cone(unit_directions):
        """
        Find the minimum cone that encompasses multiple unit vectors
        
        INPUT:
        unit_directions:
            list of unit direction. MUST HAVE AT LEAST 4 vectors.

        RETURN:
        ---
        centorid: direction
        max_angle: opening angle

        """

        # Step 2: Compute the convex hull
        hull = ConvexHull(unit_directions)


        # Step 3: Compute the centroid of the convex hull (cone axis)
        hull_vertices = unit_directions[hull.vertices]
        centroid = np.mean(hull_vertices, axis=0)
        centroid /= np.linalg.norm(centroid)  # Normalize to unit vector

        # Step 4: Calculate the maximum angle (cone aperture) between centroid and hull vertices
        angles = np.arccos(np.clip(np.dot(hull_vertices, centroid), -1.0, 1.0))
        max_angle = np.max(angles)

        return centroid, max_angle    

    @staticmethod
    def vertex_updown(vertices, tracks):
        """
        Find index of vertex that has all tracks going up
        """
        upward_inds = []
        for i in range(len(vertices)):
            vertex = vertices[i]
            isupward=True
            for i in vertex.tracks:
                # All tracks need to go upwards
                if tracks[i].At<0:
                    isupward=False
                    break
            if isupward:
                upward_inds.append(i)

        return upward_inds

    @staticmethod
    def vertex_opening_angle(vertices, tracks, return_all=False):
        """
        Get the largest opening angle of all vetices, and the axis of the corresponding cone

        RETURN:
        ---
        open_angle: float
        direction: [kx, ky, kz]
        """
        open_angles = []
        directions = []
        for k1 in range(len(vertices)):
            vertex = vertices[k1]
            trackindices = vertex.tracks
            open_angles_i = []
            directions_i = []
            
            ## This is wrong when more than 3 tracks. Max angle between any two is not the encompassing cone
            if len(trackindices)<=3:
                # Calculate opening angle for all combination of tracks within vertex
                combolist = list(itertools.combinations(trackindices, 2))
                for combo in combolist:
                    combo = [int(combo[0]),int(combo[1])]

                    track1, track2 = tracks[combo[0]], tracks[combo[1]]
                    v1 = [track1.Ax/track1.At, 1/track1.At, track1.Az/track1.At]
                    v2 = [track2.Ax/track2.At, 1/track2.At, track2.Az/track2.At]

                    mag1 = np.linalg.norm(v1)
                    mag2 = np.linalg.norm(v2)
                    cos_opening_angle = np.dot(v1, v2)/( mag1*mag2 )

                    direction = v1/mag1 + v2/mag2
                    direction = direction/np.linalg.norm(direction)
                    open_angles_i.append(np.arccos(cos_opening_angle))
                    directions_i.append(direction)

                ind_i=np.argmax(open_angles_i)
                open_angles.append(open_angles_i[ind_i])
                directions.append(directions_i[ind_i])
            
            ## Correct way for >=4 tracks
            else:
                track_directions = [[tracks[ind].Ax/tracks[ind].At, tracks[ind].Ay/tracks[ind].At, tracks[ind].Az/tracks[ind].At] for ind in trackindices]
                track_directions = np.array([a/np.linalg.norm(a) for a in track_directions])
                
                try:
                    centroid, max_angle = find_cone(track_directions)
                except KeyboardInterrupt:
                    break
                except:
                    centroid=[0,0,0]
                    max_angle=0
                    

                open_angles.append(max_angle*2)
                directions.append(centroid)            

        if return_all:
            return open_angles, directions
        else:
            # Get largest opening angle
            ind=np.argmax(open_angles)
            return open_angles[ind], directions[ind]
    
    @staticmethod
    def vertex_conatins_cms(vertices, tracks, margin=0):
        """check if any vertex contains the direction of CMS IP
        Return true if ANY vertex contains the CMS direction
        """    
        open_angles, directions = Cut.vertex_opening_angle(vertices, tracks, return_all=True)
        inds = []
        for i in range(len(directions)):
            direction_cms = np.array([vertices[i].x0, vertices[i].y0, vertices[i].z0]) - Cut.CONST_CMS_LOC
            open_angle = open_angles[i]
            angle_diff = Cut.angle(direction_cms, directions[i])
            if angle_diff<open_angle+margin:
                inds.append(i)
            
        return inds
            

    @staticmethod
    def vertex_n_tracks(vertices, tracks):
        ntrack=[len(vertex.tracks) for vertex in vertices]
        return max(ntrack)

    @staticmethod
    def vertex_min_chi2(vertices, tracks):
        chi2s=[vertex.chi2 for vertex in vertices]
        return min(chi2s)

    @staticmethod
    def vertex_inbox(vertices, tracks):
        inbox_inds=[]
        for i in range(len(vertices)):
            vertex=vertices[i]
            if det.inBox(vertex.x0,vertex.y0,vertex.z0):
                inbox_inds.append(i)

        return inbox_inds
    
    @staticmethod
    def vertex_veto_hits_in_cone(vertices, tracks, hits, opening_angles=None, directions=None, vertex_inds=None, hit_speed_cut = [25, 35]):
        """
        Check the number of hits in veto layers
        """
        if opening_angles is None or directions is None:
            opening_angles, directions = Cut.vertex_opening_angle(vertices, tracks, return_all=True)
        if vertex_inds is None:
            vertex_inds = np.arange(len(vertices))
            
        n_veto_hits = []
        # For each vertex, count number of hits that is within the cone
        for i in vertex_inds:
            vertex_position = np.array([vertices[i].x0, vertices[i].y0, vertices[i].z0])
            vertex_time = vertices[i].t0
            vertex_direction = directions[i]
            vertex_halfangle = opening_angles[i]/2
            nhits=0
            for hit in hits:
                if hit.layer not in Cut.CONST_VETO_LAYER:
                    continue
                hit_direction = (vertex_position - np.array([hit.x, hit.y, hit.z])) / (vertex_time - hit.t)
                hit_speed = np.linalg.norm(hit_direction)
                
                # # Randomly drop the noise hits
                # hit_type = hit.type
                
                if (Cut.angle(hit_direction, vertex_direction) < vertex_halfangle) \
                    and hit_speed>hit_speed_cut[0] and  hit_speed<hit_speed_cut[1]:
                    nhits+=1
            n_veto_hits.append(nhits)
        return n_veto_hits 
    
    @staticmethod
    def vertex_veto_hits_dist(vertices, tracks, hits):
        default_dist = float('inf')
        vertex_inds = np.arange(len(vertices))

        veto_dists = []
        for i in vertex_inds: 
            # Loop all tracks
            v = vertices[i]
            
            dists_this_vertex = []
            for itrack in v.tracks:
                track = tracks[itrack]
                dists_this_track = []
                for hit in hits:
                    if hit.layer not in Cut.CONST_VETO_LAYER:
                        continue
                    dists_this_track.append(Cut.track_to_point_dist(Cut.track_format(track), [hit.x, hit.y, hit.z, hit.t]))
                
                min_dist_track = min(dists_this_track) if len(dists_this_track)>0 else default_dist
                dists_this_vertex.append(min_dist_track)
            min_dist_vert = min(dists_this_vertex) if len(dists_this_vertex)>0 else default_dist
            veto_dists.append(min_dist_vert)
        return veto_dists
    
    
    
            
                    
    @staticmethod
    def track_format(track_km):
        """Change the parameter format of the track from (x0,y0,z0,t0,Ax, Ay, Az, At) to (x0,y0,z0,vx, vy, vz,t0)"""
        return [track_km.x0, track_km.y0, track_km.z0, track_km.Ax/track_km.At, track_km.Ay/track_km.At, track_km.Az/track_km.At, track_km.t0]
    
    @staticmethod
    def closest_approach(tr1, tr2):
        """
        Input:
        ---
        tr1,tr2:
            lists of track paramters [x0, y0, x0, vx, vy, vz, t0]

        return:
        ---
        midpoint([x,y,z,t]), distance
        """

        r1, v1, t1 = tr1[:3], tr1[3:6], tr1[:6]
        r2, v2, t2 = tr2[:3], tr2[3:6], tr2[:6]
        
        # Relative speed
        rel_v = v2 - v1
        rel_v2 = np.dot(rel_v, rel_v) 
        # Relative position
        displacement = r1 - r2; 
        
        # Midpoint time
        t_ca = (  np.dot(displacement, rel_v) + np.dot((v2*t2 - v1*t1), rel_v)  )/rel_v2;    
        
        # Midpoint position
        pos1 = r1 + v1*(t_ca - t1);
        pos2 = r2 + v2*(t_ca - t2);
        r_ca = (pos1 + pos2)*(0.5)
        
        # Closest distance
        line_distance = np.linalg.norm((pos1- pos2))

        return np.concatenate((r_ca, [t_ca])), line_distance



    @staticmethod
    def track_to_point_dist(track, point):
        v1 = np.array(track[3:6])
        r1 = np.array(track[:3])
        t1 = np.array(track[6])
        
        r0 = point[:3]
        t0 = point[3]
        pos = r1 + v1*(t1-t0)
        dist = np.linalg.norm(pos-r0)
        return dist
        
