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

sys.path.append("/project/rrg-mdiamond/owhgabri/MATHUSLA_JupyterNoteBooks/tracker")
sys.path.insert(1, "/project/rrg-mdiamond/owhgabri/MATHUSLA_JupyterNoteBooks/tracker")
os.chdir('/project/rrg-mdiamond/owhgabri/pyTracker')
print(os.getcwd())
print(joblib.__version__)

# ROOT
import ROOT as root

import pprint

import kalmanfilter as KF
import utilities as Util
import trackfinder as TF
import vertexfinder as VF
import datatypes
from datatypes import *

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
# y_layer_starts = [8550, 8631.6, 9730, 9811.6, 9893.2, 9974.8, 10056.4, 10138] # 6 layer 
# z_wall_starts = [10900, 10981.6, 11063.2, 11144.8, 11226.4, 11308] # 6 layer
# back_top =  y_layer_starts[4] # 6 layer
# -----4 layer--------------------
y_layer_starts = [8550, 8631.6, 9893.2, 9974.8, 10056.4, 10138] # 4 layer
z_wall_starts = [10900, 10981.6, 11063.2, 11144.8, 11226.4, 11308] # 4 layer
back_top =  y_layer_starts[2] # 4 layer

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
        if proj_point[1] < y_bottoms[0] or \
        proj_point[0] < x_lims[0] or proj_point[0] > x_lims[1] or \
        proj_point[2] < z_fronts[0]:
            return None
        proj_cov = FW_Veto.ProjectionCovariance(track, TrackPoints, proj_time)
        proj_cov += np.diag(hit_err)
        residual = proj_point - hit_pos 
        try:
            inv_cov = np.linalg.inv(proj_cov)
        except:
            return None
        return residual.T @ inv_cov @ residual

    def Chi2Veto(vertices, tracks, hits, cutoff):
        """
        Decide if the event is vetoed based on the chi2 cutoff
        Return True if vetoed. Return false if not vetoed
        """
        wf_hits = []
        for hit in hits:
            if hit.layer < 4:
                wf_hits.append(hit)
        if len(wf_hits) == 0: # Event passes
            return False
        min_chi2 = None
        for vertex in vertices:
            v_tracks = []
            for t_index in vertex.tracks:
                v_tracks.append(tracks[t_index])
            for track in v_tracks:
                for hit in wf_hits:
                    cur_chi2 = FW_Veto.GetChiSquared(hit,track)
                    if cur_chi2 is not None and (min_chi2 is None or cur_chi2 < min_chi2):
                        min_chi2 = cur_chi2
        if min_chi2 is None:
            return False
        chi2_red = min_chi2
        if chi2_red > 0 and chi2_red < cutoff:
            return True
        return False
 
    def DistanceVeto(vertices, tracks, hits, cutoff):
        """
        Decide if the event is vetoed based on the distance cutoff
        Return True if vetoed. Return false if not vetoed
        """        
        wf_hits = []
        for hit in hits:
            if hit.layer < 4:
                wf_hits.append(hit)
        if len(wf_hits) == 0: # Event passes 
            return False
        min_dist = None
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
                    if proj_point is not None:
                        dist = GetDistance(hit_pos, proj_point)
                    if (min_dist is None or dist < min_dist):
                        min_dist = dist
        if min_dist is None or min_dist < cutoff:
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

    def DistanceVeto(vertices, tracks, hits, cutoff):
        """
        Decide if the event is vetoed based on the distance cutoff
        Return True if vetoed. Return false if not vetoed
        """      
        wf_hits = []
        for hit in hits:
            if hit.layer < 4:
                wf_hits.append(hit)
        if len(wf_hits) == 0:
            return False
        min_dist = None
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
        if min_dist is None or min_dist < cutoff:
            return True
        return False

#----------------------------------------------------------------------------------------
class MaterialVeto:
    """
    This only vetoes if all vertices are within the cutoff range
    """

    def MakeLayers():
        """
        Returns a list of the layer modules of 
        the form
        ((x_min, x_max),(y_min,y_max),(z_min, z_max))
        """
        layers = []
        for i in range(2, len(y_layer_starts)):
            for j in range(len(x_mod_starts)):
                for k in range(len(z_mod_starts)):
                    y_min = y_layer_starts[i]; y_max = y_layer_starts[i] + thickness
                    x_min = x_mod_starts[j]; x_max = x_mod_starts[j] + mod_length
                    z_min = z_mod_starts[k]; z_max = z_mod_starts[k] + mod_length
                    layers.append(np.array(((x_min, x_max),(y_min,y_max),(z_min,z_max))))
        return layers
    
    def MakeBacks():
        """
        Returns a list of the back modules of 
        the form
        ((x_min, x_max),(y_min,y_max),(z_min, z_max))
        """
        backs = []
        y_min = back_top - mod_length; y_max = back_top
        for j in range(len(x_mod_starts)):
            for k in range(len(z_wall_starts)):
                x_min = x_mod_starts[j]; x_max = x_mod_starts[j] + mod_length
                z_min = z_wall_starts[k]; z_max = z_wall_starts[k] + thickness
                backs.append(np.array(((x_min, x_max),(y_min,y_max),(z_min,z_max))))
        return backs
    
    def MakeWalls():
        """
        Returns a list of the hermetic walls of the form
        ((x_min, x_max),(y_min,y_max),(z_min, z_max))
        """
        walls = []
        y_min = y_layer_starts[0]; y_max = back_top
        x_min = x_mod_starts[0]; x_max = x_mod_starts[-1] + mod_length
        for i in range(len(front_wall_starts)):
            z_min = front_wall_starts[i]; z_max = front_wall_starts[i] + thickness
            walls.append(np.array(((x_min,x_max),(y_min,y_max),(z_min,z_max))))
        return walls
    
    def MakeFloors():
        """
        Returns a list of the hermetic walls of the form
        ((x_min, x_max),(y_min,y_max),(z_min, z_max))
        """
        floors = []
        z_min = z_mod_starts[0]; z_max = z_mod_starts[-1] + mod_length
        x_min = x_mod_starts[0]; x_max = x_mod_starts[-1] + mod_length
        for i in range(n_floors):
            y_min = y_layer_starts[i]; y_max = y_layer_starts[i] + thickness
            floors.append(np.array(((x_min,x_max),(y_min,y_max),(z_min,z_max))))
        return floors
        
    def MakeBeams():
        """
        Returns a list with beams of the form 
        ((x_min, x_max),(y_min,y_max),(z_min, z_max))
        """
        beams = []
        for i in range(len(y_layer_starts) - 1):
            for j in range(len(x_mod_starts)):
                for k in range(len(z_mod_starts)):
                    y_min = y_layer_starts[i] + thickness; y_max = y_layer_starts[i+1]
                    
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

    def MakeGround():
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
    
    def MakeMaterials():
        """
        Returns all the materials concatenated together
        """
        beams = MaterialVeto.MakeBeams()
        layers = MaterialVeto.MakeLayers()
        backs = MaterialVeto.MakeBacks()
        walls = MaterialVeto.MakeWalls()
        floors = MaterialVeto.MakeFloors()
        ground = MaterialVeto.MakeGround()
        return np.concatenate((beams,layers,backs,walls,floors, ground))
    
    def GetResidual(vertex, rectangle):
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
    
    def GetChi2(vertex, rectangle):
        """
        Get the chi2 between the vertex and the rectangle
        """
        res = MaterialVeto.GetResidual(vertex, rectangle)
        trimmed_cov = np.delete(vertex.cov, 3, 0)
        trimmed_cov = np.delete(trimmed_cov, 3, 1)
        return res.T@np.linalg.inv(trimmed_cov)@res
    
    def MinDistance(vertex, materials):
        """
        Returns the minimum distance to a material
        """
        min_dist = None
        for material in materials:
            res = MaterialVeto.GetResidual(vertex, material)
            dist = np.sqrt(res[0]**2 + res[1]**2 + res[2]**2)
            if min_dist is None or dist < min_dist:
                min_dist = dist
        return min_dist
    
    def MinChi2(vertex):
        """
        Returns the minimum chi2 of the vertex
        """
        materials = MakeMaterials()
        min_chi2 = None
        for material in materials:
            chi2 = MaterialVeto.GetChi2(vertex, material)
            if min_chi2 is None or chi2 < min_chi2:
                min_chi2 = chi2
        return min_chi2

    def DistanceVeto(vertices, materials, cutoff, strict=False):
        """
        Vetoes an event if all vertices are within the cutoff distance from a material.
        Return True if vetoed. Return false if not vetoed
        if strict: veto if a single vertex is within the cutoff distance from a material
        """
        max_dist = None # Max of the minimum distances
        min_dist = None #
        for vertex in vertices:
            dist = MaterialVeto.MinDistance(vertex, materials)
            if max_dist is None or dist > max_dist:
                max_dist = dist
            if min_dist is None or dist < min_dist:
                min_dist = dist
        if not strict and max_dist is None or max_dist > cutoff:
            return False
        elif strict and min_dist is None or min_dist > cutoff:
            retrn False
        return True

#----------------------------------------------------------------------------------------
def CutoffVeto(vertices, strict=False):
    """
    Vetoes an event based on vertex being in decay volume
    Return True if vetoed. Return false if not vetoed
    if strict: veto if a single vertex is outside the cutoff
    """
    for vertex in vertices:
        if not (vertex.y0 < y_layer_starts[2] and \
        vertex.y0 > 8500 and vertex.x0 > -2000 and \
        vertex.x0 < 2000 and vertex.z0 > 6800 and \
        vertex.z0 < 11350): # One outside 
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
    if strict: # 
        return False

#----------------------------------------------------------------------------------------
def UpVertexVeto(vertices, tracks):
        """
        Vetoes an event based on whether there is a vertex with
        downward going tracks. Return True if vetoed. Return false if not vetoed
        """
        for vertex in vertices:
            for index in vertex.tracks:
                track = tracks[index]
                hits = track.hits_filtered
                point1 = hits[0]
                point2 = hits[-1]
                if point1[3] > point2[3]: # point 1 is later
                    ydiff = point1[1] - point2[1]
                    if ydiff < 1:
                        return True
                else: #Point 2 is layter
                    ydiff = point2[1] - point1[1]
                    if ydiff < 1:
                        return True
        return False


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



