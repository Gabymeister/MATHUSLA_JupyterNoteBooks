import numpy as np
import ROOT as root
from joblib import dump, load
import matplotlib.pyplot as plt
#from scipy.stats import chisquare
from scipy import stats
from array import array

import physics
import visualization
from detector import Detector
import util

class Event:
    """
    Process and plot simulation event
    
    ...
    
    Attributes
    ----------
    
    
    Methods
    ----------
    ExtractTruthPhysics():
        Read the truth of one simulated event
    
    
    """
    TRUTH_PARTICLE_E_THRESHOLD = 10 #MeV

    t0 = 0.0 #ns, Start time of one event
    tm = 0.0 #ns, End time of one event
    globalTime = 0.0 #ns 
    timeResolution = 1.0 #ns
    timeRangeCut = 0.1 #ns

    tracksAtGlobalTime = []

    recoTrackList = []
    recoVertexList = []
    digiHitList = []
    simHitList = []

    def __init__(self, filename, EventNumber, tree_name="integral_tree", ):
        self.EventNumber = EventNumber
        self.tfile = root.TFile.Open(filename)
        self.Tree_nm = filename
        self.Tree = self.tfile.Get(tree_name)
        self.truthTrackList = []
        self.writeDirectory = './'

        self.kalman = True
        self.unused = False
        self.used = False
        self.oldvertex = False
        self.vertex = False
        self.merged = False
        self.vert_vel = True


    def ExtractTruthPhysics(self):
        
        self.Tree.GetEntry(self.EventNumber)
        self.truthTrackList=[]
        particleSet = set()
        excluded_pdgs = set([]) # PDG particles to exclude

        for hitn in range(int(self.Tree.NumHits)):
            if self.Tree.Hit_particleEnergy.at(hitn) > self.TRUTH_PARTICLE_E_THRESHOLD:
                particleSet.add( physics.Particle( int(self.Tree.Hit_G4TrackId.at(hitn)), int(self.Tree.Hit_particlePdgId.at(hitn)), int(self.Tree.Hit_G4ParentTrackId.at(hitn)) ))

        for particle in particleSet:
            if not particle.pdgID in excluded_pdgs:
                currentTrack = physics.Track(particle)
                for hitn in range(int(self.Tree.NumHits)):
                    if ( int(self.Tree.Hit_G4TrackId.at(hitn)) == particle.trackID):
                        time = self.Tree.Hit_time.at(hitn)
                        location = physics.Vector(self.Tree.Hit_x.at(hitn), self.Tree.Hit_y.at(hitn), self.Tree.Hit_z.at(hitn))
                        energy = self.Tree.Hit_particleEnergy.at(hitn)
                        momentum = physics.Vector( self.Tree.Hit_particlePx.at(hitn), self.Tree.Hit_particlePy.at(hitn), self.Tree.Hit_particlePz.at(hitn) )
                        point = physics.TrackPoint(time, location, energy, momentum)
                        currentTrack.AddPoint(point)

            if (len(currentTrack.pointList)) <= 2:
                continue

            if currentTrack.TimeRange() < self.timeRangeCut:
                continue
            currentTrack.pointList.sort()

            self.truthTrackList.append(currentTrack)

        self.ResetTracks()
        self.truthTrackList.sort()
#        self.t0 = min([track.pointList[0].time for track in self.truthTrackList])
#        self.tm = max([track.pointList[len(track.pointList)-1].time for track in self.truthTrackList])

    def ExtractTruthPhysics_list(self):
        
        """
        Return a list of
        [x,y,z,t, PID, Energy, TRACK_ID]
        """

        self.Tree.GetEntry(self.EventNumber)
        self.truthTrackList_list=[]
        particleSet = set()

        track_ids = np.array(util.c2list(self.Tree.Hit_G4TrackId))
        hits_x = np.array(util.c2list(self.Tree.Hit_x))
        hits_y = np.array(util.c2list(self.Tree.Hit_y))
        hits_z = np.array(util.c2list(self.Tree.Hit_z))
        hits_t = np.array(util.c2list(self.Tree.Hit_time))
        hits_pid = np.array(util.c2list(self.Tree.Hit_particlePdgId))
        hits_energy = np.array(util.c2list(self.Tree.Hit_particleEnergy))

        track_ids_unique = np.unique(track_ids)
        hit_inds = np.arange(len(track_ids))

        for trackid in track_ids_unique:
            hit_inds_track = hit_inds[track_ids==trackid]

            # each hit is a list of [x, y, z, t, pid, energy]
            currentTrack = [hits_x[hit_inds_track].tolist(),hits_y[hit_inds_track].tolist(),hits_z[hit_inds_track].tolist(),hits_t[hit_inds_track].tolist(),hits_pid[hit_inds_track].tolist(),hits_energy[hit_inds_track].tolist(),track_ids[hit_inds_track]]

            # Require >2 hits
            if (len(currentTrack[0])) <= 2:
                continue

            # Require time range
            if currentTrack[3][-1]-currentTrack[3][0] < self.timeRangeCut:
                continue

            # sort hits by time
            current_track_sort  = util.sortbyrow(np.array(currentTrack), 3) 
            self.truthTrackList_list.append(current_track_sort)

        return self.truthTrackList_list
   
    def Print(self):

        self.Tree.SetBranchStatus("GenParticle_G4index",1)
        self.Tree.SetBranchStatus("GenParticle_x",1) # start position of particles
        self.Tree.SetBranchStatus("GenParticle_y",1)
        self.Tree.SetBranchStatus("GenParticle_z",1)

        visEngine_local = visualization.Visualizer()

        visEngine_local.writeDirectory = self.writeDirectory

        # try:
#         if 1:
#             used_gens_inds = np.where(np.array(self.Tree.GenParticle_G4index) != -1)[0]
#             #dump(used_gens_inds,"used_gens_inds.joblib")
#             #used_gens_inds = load("used_gens_inds.joblib")
  
#             if len(used_gens_inds) != 0:
#                 vert_truth = [self.Tree.GenParticle_y[int(used_gens_inds[0])] / 10,
#                                     self.Tree.GenParticle_x[int(used_gens_inds[0])] / 10,
#                                     self.Tree.GenParticle_z[int(used_gens_inds[0])] / 10] # only first since all used
  
#                         #truth vertex location
#                 visEngine_local.AddVertex( {'point':[vert_truth[0], vert_truth[1], vert_truth[2]], 'col':'tab:green', 'vert vel':[[],[],[]]} )

        # except:
        #     pass

        num_truth = 0
#        list_of_truth_hits = []
        for track in self.truthTrackList:
            x = []
            y = []
            z = []
            # print(track)
            for point in track.pointList:
                x.append(point.location.x)
                y.append(point.location.y)
                z.append(point.location.z)
#                temp.append([x, y, z] )

#            list_of_truth_hits.append(temp)
#            print(y)
            visEngine_local.TrackDisplayPoints(x, y, z, track.color(), track.LabelString())

            num_truth += 1

        self.num_truth = num_truth
        visEngine_local.Draw(self.Tree_nm.split('/')[-1].split('.')[0]+'_{}_{}.pdf'.format(self.EventNumber,'truth'))



    def ResetTracks(self):
        self.tracksAtGlobalTime = [ [] for x in range(len(self.truthTrackList))  ]
        self.globalTime = self.t0

    def StepParticles(self):
        self.globalTime += self.timeResolution

        for n, track in enumerate(self.truthTrackList):
            self.tracksAtGlobalTime[n].append(track.PointAtTime(self.globalTime))


    #returns location of particles at a given time t, called by visulaization engine
    def TruthAtTime(self, t):
        if ( np.absolute(self.globalTime - t) < self.timeResolution ):
            return self.tracksAtGlobalTime

        if (t < self.globalTime):
            self.globalTime = self.t0
            self.ResetTracks()

        while (self.globalTime < t):
            self.StepParticles()

        return self.tracksAtGlobalTime

    def find_with_bool(self, bool_func, Op=0):
        ''' find indices of events that satisfy a criteria given by bool_func (boolean valued function)'''

        inds = []

        for ev in range(int(self.Tree.GetEntries())):
            self.Tree.GetEntry(ev)

            if bool_func(self.Tree, op=Op):
                inds.append(ev)

        return inds


    def RecoLinVertex(self):
        ''' visualise tracks and vertices made by ML Tracker, and the digi hits and truth decay locations '''
        

        visEngine_local = visualization.Visualizer()

        visEngine_local.vert_vel = self.vert_vel

        visEngine_local.writeDirectory = self.writeDirectory

        det = Detector()

        self.Tree.SetBranchStatus("GenParticle_G4index",1)
        self.Tree.SetBranchStatus("GenParticle_x",1) # start position of particles
        self.Tree.SetBranchStatus("GenParticle_y",1)
        self.Tree.SetBranchStatus("GenParticle_z",1)

        self.Tree.SetBranchStatus("Track_x0", 1)
        self.Tree.SetBranchStatus("Track_y0", 1)
        self.Tree.SetBranchStatus("Track_z0", 1)
        self.Tree.SetBranchStatus("Track_velX", 1)
        self.Tree.SetBranchStatus("Track_velY", 1)
        self.Tree.SetBranchStatus("Track_velZ", 1)

        self.Tree.SetBranchStatus("Vertex_x", 1)
        self.Tree.SetBranchStatus("Vertex_y", 1)
        self.Tree.SetBranchStatus("Vertex_z", 1)


        self.Tree.GetEntry(self.EventNumber)

        digi_hit_inds = util.unzip(self.Tree.Track_hitIndices)

        used_inds = set(self.Tree.Track_hitIndices)

        numtracks = int(self.Tree.NumTracks)

        used_track_inds = set()

        vert_inds = util.unzip(self.Tree.Vertex_trackIndices)

        # get truth vertex info (decay location of parent particle used in the event)
        used_gens_inds = np.where(np.array(self.Tree.GenParticle_G4index) != -1)[0]
                #dump(used_gens_inds,"used_gens_inds.joblib")
                #used_gens_inds = load("used_gens_inds.joblib")

        if len(used_gens_inds) != 0:
                        # Vertex truth location (translated from Pythia coordinate system)
            vert_truth = [self.Tree.GenParticle_y[int(used_gens_inds[0])] / 10,
                                      self.Tree.GenParticle_x[int(used_gens_inds[0])] / 10,
                                      self.Tree.GenParticle_z[int(used_gens_inds[0])] / 10]

                        # pass truth vertex to visualiser
            visEngine_local.AddVertex( {'point':[vert_truth[0], vert_truth[1], vert_truth[2]], 'col':'tab:green', 'vert vel':[[],[],[]]} )

        colors = ["tab:"+col for col in ['blue','orange','red','purple','brown','pink','olive','cyan']]

        for n in range(len(self.Tree.Vertex_x)):
            c = n if n < len(colors) else n%len(colors)

            xv = self.Tree.Vertex_x[n]
            yv = self.Tree.Vertex_y[n]
            zv = self.Tree.Vertex_z[n]

            visEngine_local.AddVertex( {'point':[xv, yv, zv], 'col':colors[c], 'vert vel':[[]*3]} )

            for trk_ind in vert_inds[n]:
                if self.used:
                    for ind in digi_hit_inds[int(trk_ind)]:
                        x = self.Tree.Digi_x[ind]
                        y = self.Tree.Digi_y[ind]
                        z = self.Tree.Digi_z[ind]
                        visEngine_local.AddPoint( [[x, y, z], colors[c]] )


                x0, y0, z0 = self.Tree.Track_x0[int(trk_ind)], self.Tree.Track_y0[int(trk_ind)], self.Tree.Track_z0[int(trk_ind)]
                vx, vy, vz = self.Tree.Track_velX[int(trk_ind)], self.Tree.Track_velY[int(trk_ind)], self.Tree.Track_velZ[int(trk_ind)]

                [xi, yi, zi] = det.FindIntercept(x0, y0, z0, vx, vy, vz) # find intercept with boundary

                visEngine_local.TrackDisplayPoints([x0,xi], [y0,yi], [z0,zi], color=colors[c])

                used_track_inds.add(int(trk_ind))


        if self.unused:
            for n in range(len(self.Tree.Digi_numHits)):
                if not (n in used_inds):
                    x = self.Tree.Digi_x[n]
                    y = self.Tree.Digi_y[n]
                    z = self.Tree.Digi_z[n]
                    visEngine_local.AddHit( [x, y, z] )

            #print(used_track_inds)

            for n in range(int(self.Tree.NumTracks)):
                if not (n in used_track_inds):

                    x0, y0, z0 = self.Tree.Track_x0[n], self.Tree.Track_y0[n], self.Tree.Track_z0[n]
                    vx, vy, vz = self.Tree.Track_velX[n], self.Tree.Track_velY[n], self.Tree.Track_velZ[n]

                    [xi, yi, zi] = det.FindIntercept(x0, y0, z0, vx, vy, vz) # find intercept with boundary

                    visEngine_local.TrackDisplayPoints([x0,xi], [y0,yi], [z0,zi], color='k', opac=0.2)

                    for ind in digi_hit_inds[n]:
                        x = self.Tree.Digi_x[ind]
                        y = self.Tree.Digi_y[ind]
                        z = self.Tree.Digi_z[ind]
                        visEngine_local.AddPoint( [[x, y, z], 'tab:gray'] )




        visEngine_local.Draw(self.Tree_nm.split('/')[-1].split('.')[0]+'_{}_{}.pdf'.format(self.EventNumber,'_linear'))


    def RecoVertex(self):
        ''' visualise tracks and vertices made by kalman filter, and the digi hits and truth decay locations '''

        visEngine_local = visualization.Visualizer()

        visEngine_local.vert_vel = self.vert_vel

        visEngine_local.writeDirectory = self.writeDirectory

        self.Tree.SetBranchStatus("GenParticle_G4index",1)
        self.Tree.SetBranchStatus("GenParticle_x",1) # start position of particles
        self.Tree.SetBranchStatus("GenParticle_y",1)
        self.Tree.SetBranchStatus("GenParticle_z",1)

        self.Tree.SetBranchStatus("Vertex_k_m_x",1)
        self.Tree.SetBranchStatus("Vertex_k_m_y",1)
        self.Tree.SetBranchStatus("Vertex_k_m_z",1)
        self.Tree.SetBranchStatus("Vertex_k_m_trackIndices",1)

        #self.Tree.SetBranchStatus("king_move_inds", 1)
        self.Tree.SetBranchStatus("x_estimates_m", 1)
        self.Tree.SetBranchStatus("y_estimates_m", 1)
        self.Tree.SetBranchStatus("z_estimates_m", 1)
        self.Tree.SetBranchStatus("NumTracks_k_m", 1)

        #self.Tree.SetBranchStatus("vertex_vx_m", 1)
        #self.Tree.SetBranchStatus("vertex_vy_m", 1)
        #self.Tree.SetBranchStatus("vertex_vz_m", 1)

        self.Tree.SetBranchStatus("Digi_numHits", 1)
        self.Tree.SetBranchStatus("Digi_x", 1)
        self.Tree.SetBranchStatus("Digi_y", 1)
        self.Tree.SetBranchStatus("Digi_z", 1)

        self.Tree.GetEntry(self.EventNumber)

        # indices of used digi hits (organised by track and organised into a set)
        digi_hit_inds = util.unzip(self.Tree.Track_k_m_hitIndices)
        used_inds = set(self.Tree.Track_k_m_hitIndices)

        used_track_inds = set()

        # indices of tracks used in the made vertices
        vert_inds = util.unzip(self.Tree.Vertex_k_m_trackIndices)

#        print(vert_inds)

        # track best estimates
        x_s = util.unzip(self.Tree.x_estimates_m)
        y_s = util.unzip(self.Tree.y_estimates_m)
        z_s = util.unzip(self.Tree.z_estimates_m)

        # velocity best estimates at the vertex
        #vxv = util.unzip(self.Tree.vertex_vx_m)
        #vyv = util.unzip(self.Tree.vertex_vy_m)
        #vzv = util.unzip(self.Tree.vertex_vz_m)

        # get truth vertex info (decay location of parent particle used in the event)
        used_gens_inds = np.where(np.array(self.Tree.GenParticle_G4index) != -1)[0]
        #dump(used_gens_inds,"used_gens_inds.joblib")
        #used_gens_inds = load("used_gens_inds.joblib")

        if len(used_gens_inds) != 0:
            # Vertex truth location (translated from Pythia coordinate system)
            vert_truth = [self.Tree.GenParticle_y[int(used_gens_inds[0])] / 10,
                      self.Tree.GenParticle_x[int(used_gens_inds[0])] / 10,
                      self.Tree.GenParticle_z[int(used_gens_inds[0])] / 10]

            # pass truth vertex to visualiser
            visEngine_local.AddVertex( {'point':[vert_truth[0], vert_truth[1], vert_truth[2]], 'col':'tab:green', 'vert vel':[[],[],[]]} )

        # track and vertex colors
        colors = ["tab:"+col for col in ['blue','orange','red','purple','brown','pink','olive','cyan']]

        # pass vertices, tracks, and hits used in them to visualiser
        # (use one color per vertex and its tracks and their hits)
        for n in range(len(self.Tree.Vertex_k_m_x)):
            # choose index of a color to use for visualiser
            c = n if n < len(colors) else n%len(colors)

            xv = self.Tree.Vertex_k_m_x[n]
            yv = self.Tree.Vertex_k_m_y[n]
            zv = self.Tree.Vertex_k_m_z[n]

            # made vertex passed to visualiser
            #if len(vxv) != 0:
            #    visEngine_local.AddVertex( {'point':[xv, yv, zv], 'col':colors[c], 'vert vel':[vxv[n], vyv[n], vzv[n]]} )

            #else:
            if True:
                visEngine_local.AddVertex( {'point':[xv, yv, zv], 'col':colors[c], 'vert vel':[[]*3]} )

            # loop over tracks in the current vertex
            for trk_ind in vert_inds[n]:
                if self.used:
                    for ind in digi_hit_inds[int(trk_ind)]:
                        x = self.Tree.Digi_x[ind]
                        y = self.Tree.Digi_y[ind]
                        z = self.Tree.Digi_z[ind]

                        # pass each hit used in the current track to visualiser
                        visEngine_local.AddPoint( [[x, y, z], colors[c]] )

                # pass track to visualiser
                visEngine_local.TrackDisplayPoints(x_s[int(trk_ind)], y_s[int(trk_ind)], z_s[int(trk_ind)], color=colors[c])

                used_track_inds.add(int(trk_ind))


        if self.unused:
            # pass hits that weren't used in a track to visualiser
            for n in range(len(self.Tree.Digi_numHits)):
                if not (n in used_inds):
                    x = self.Tree.Digi_x[n]
                    y = self.Tree.Digi_y[n]
                    z = self.Tree.Digi_z[n]
                    visEngine_local.AddHit( [x, y, z] )

            #print(used_track_inds)

            # pass tracks not used in a vertex and their used hits to visualiser
            for n in range(self.Tree.NumTracks_k_m):
                if not (n in used_track_inds):
                    visEngine_local.TrackDisplayPoints(x_s[n], y_s[n], z_s[n], color='k', opac=0.2)

                    #print(n)

                    for ind in digi_hit_inds[n]:
                        x = self.Tree.Digi_x[ind]
                        y = self.Tree.Digi_y[ind]
                        z = self.Tree.Digi_z[ind]
                        visEngine_local.AddPoint( [[x, y, z], 'tab:gray'] )

            #king_move_inds = util.unzip(self.Tree.king_move_inds)

            #print('king move indices are ',king_move_inds)

            # show hits dropped from the hit pool (via king moves algorithm)
            '''
            for inds in king_move_inds:
                for ind in inds:
                    x = self.Tree.Digi_x[ind]
                    y = self.Tree.Digi_y[ind]
                    z = self.Tree.Digi_z[ind]
                    visEngine_local.AddPoint( [[x, y, z], 'm'] )
            '''
        # draw the plot and save (show) it
        visEngine_local.Draw(self.Tree_nm.split('/')[-1].split('.')[0]+'_{}_{}.pdf'.format(self.EventNumber,'merged')) #_linear'))



    def VertexAccuracy(self,id_str=''):
        ''' get differences in each coordinate between decay location (truth vertex) and the closest made vertex'''

        self.Tree.SetBranchStatus("GenParticle_G4index",1)
        self.Tree.SetBranchStatus("GenParticle_x",1) # start position of particles
        self.Tree.SetBranchStatus("GenParticle_y",1)
        self.Tree.SetBranchStatus("GenParticle_z",1)
        self.Tree.SetBranchStatus("GenParticle_time",1)

        self.Tree.SetBranchStatus("Vertex_k_m_x",1)
        self.Tree.SetBranchStatus("Vertex_k_m_y",1)
        self.Tree.SetBranchStatus("Vertex_k_m_z",1)
        self.Tree.SetBranchStatus("Vertex_k_m_t",1)
        self.Tree.SetBranchStatus("Vertex_k_m_ErrorX",1)
        self.Tree.SetBranchStatus("Vertex_k_m_ErrorY",1)
        self.Tree.SetBranchStatus("Vertex_k_m_ErrorZ",1)

        dx, dy, dz, dt = [], [], [], []
        ex, ey, ez = [], [], []
        px, py, pz = [], [], []
        chi = []

        long_res, trans_res = [], []

        max = 1e2

        res_info_txt = []

        for ev_num in range(self.Tree.GetEntries()):
        #for ev_num in [34, 38, 52]:

            print("event {}".format(ev_num)) if ev_num % 1000 == 0 else None

            self.Tree.GetEntry(ev_num)

            used_gens_inds = np.where(np.array(self.Tree.GenParticle_G4index) != -1)[0]
            #dump(used_gens_inds,"used_gens_inds.joblib")
            #used_gens_inds = load("used_gens_inds.joblib")

            # particles form one vertex only first since all used
            # Conversion from Gen (Pythia) info to CMS: all in mm, y <-> x, same origin
            # in CMS: x, y, z, t
            vert_truth = [self.Tree.GenParticle_y[int(used_gens_inds[0])] / 10,
                      self.Tree.GenParticle_x[int(used_gens_inds[0])] / 10,
                      self.Tree.GenParticle_z[int(used_gens_inds[0])] / 10,
                      self.Tree.GenParticle_time[int(used_gens_inds[0])] / 300]
            '''
            vert_truth = [self.Tree.GenParticle_y[int(used_gens_inds[0])],
                      self.Tree.GenParticle_x[int(used_gens_inds[0])],
                      self.Tree.GenParticle_z[int(used_gens_inds[0])],
                      self.Tree.GenParticle_time[int(used_gens_inds[0])]]
            '''

            x = self.Tree.Vertex_k_m_x
            y = self.Tree.Vertex_k_m_y
            z = self.Tree.Vertex_k_m_z
            t = self.Tree.Vertex_k_m_t
            _ex = self.Tree.Vertex_k_m_ErrorX
            _ey = self.Tree.Vertex_k_m_ErrorY
            _ez = self.Tree.Vertex_k_m_ErrorZ

            min_dist = 1e6
            min_index = -1

            for vert in range(len(x)):
                if self.Tree.NumTracks_k_m > 0 or id_str not in ['1','3']:
                        _dx = x[vert] - vert_truth[0]
                        _dy = y[vert] - vert_truth[1]
                        _dz = z[vert] - vert_truth[2]

                        dr = np.sqrt(_dx**2 + _dy**2 + _dz**2)

                        if dr < min_dist:
                            min_index = vert

            if min_index != -1:
                dx.append(x[min_index] - vert_truth[0])
                dy.append(y[min_index] - vert_truth[1])
                dz.append(z[min_index] - vert_truth[2])
                dt.append(t[min_index] - vert_truth[3])

                # vIP is given by vertex position
                # want to seperate out longitudinal and transverse
                # project res onto vIP -> longitudinal
                # subtract longitudinal from res -> transverse

                vIP = np.array(vert_truth[:3])
                res = np.array([dx[-1],dy[-1],dz[-1]])

                vIP_norm = np.sum(vIP * vIP)
                _long_res = np.sum(res * vIP) / np.sqrt(vIP_norm) # signed to indicated direction
                long_proj = _long_res / np.sqrt(vIP_norm) * vIP

                trans_proj = res - long_proj
                _trans_res = np.sqrt(np.sum(trans_proj * trans_proj))
                if trans_proj[0] < 0: # assign negative if x residual is negative (for gaussian fit)
                        _trans_res *= -1.0

                long_res.append(_long_res)
                trans_res.append(_trans_res)

                ex.append(_ex[min_index])
                ey.append(_ey[min_index])
                ez.append(_ez[min_index])
                px.append(dx[-1] / ex[-1])
                py.append(dy[-1] / ey[-1])
                pz.append(dz[-1] / ez[-1])
                #chi.append(dx[-1]**2 / ex[-1]**2 + dy[-1]**2 / ey[-1]**2 + dz[-1]**2 / ez[-1]**2)
                '''
                res_info_txt.append("\nev num is " + str(ev_num)+ '\n')
                res_info_txt.append("chi is " + str(chi[-1]) + '\n')
                res_info_txt.append("pull X is " + str(px[-1]) + '\n')
                res_info_txt.append("pull Y is " + str(py[-1]) + '\n')
                res_info_txt.append("pull Z is " + str(pz[-1]) + '\n')
                res_info_txt.append("error X is " + str(ex[-1]) + '\n')
                res_info_txt.append("error Y is " + str(ey[-1]) + '\n')
                res_info_txt.append("error Z is " + str(ez[-1]) + '\n')
                res_info_txt.append("diff X is " + str(dx[-1]) + '\n')
                res_info_txt.append("diff Y is " + str(dy[-1]) + '\n')
                res_info_txt.append("diff Z is " + str(dz[-1]) + '\n')
                '''
        #f = open("res_info.txt","w+")
        #f.writelines(res_info_txt)
        #f.close()
        '''
        vis_engine = visualization.Histogram(dx, rng=(-max,max), Title='Vertex X Resolution', \
                    xaxis='Vertex X - Decay Location X [cm]', fname='resolutionX.png')
        vis_engine = visualization.Histogram(dy, rng=(-max,max), Title='Vertex Y Resolution', \
                    xaxis='Vertex Y - Decay Location Y [cm]', fname='resolutionY.png')
        vis_engine = visualization.Histogram(dz, rng=(-max,max), Title='Vertex Z Resolution', \
                    xaxis='Vertex Z - Decay Location Z [cm]', fname='resolutionZ.png')
        vis_engine = visualization.Histogram(ex, rng=(0,max), Title='Vertex X Error', \
                    xaxis='Error X [cm]', fname='errorX.png')
        vis_engine = visualization.Histogram(ey, rng=(0,max), Title='Vertex Y Error', \
                    xaxis='Error Y [cm]', fname='errorY.png')
        vis_engine = visualization.Histogram(ez, rng=(0,max), Title='Vertex Z Error', \
                    xaxis='Error Z [cm]', fname='errorZ.png')
        vis_engine = visualization.Histogram(px, rng=(-100,100), Title='Vertex X Pull', \
                    xaxis='Pull X []', fname='pullX.png')
        vis_engine = visualization.Histogram(py, rng=(-100,100), Title='Vertex Y Pull', \
                    xaxis='Pull Y []', fname='pullY.png')
        vis_engine = visualization.Histogram(pz, rng=(-100,100), Title='Vertex Z Pull', \
                    xaxis='Pull Z []', fname='pullZ.png')
        vis_engine = visualization.Histogram(chi, rng=(0,200), Title='Vertex chi', \
                    xaxis='Chi []', fname='chi.png')
        '''
        '''
        vis_engine = visualization.root_Histogram(dx, rng=(-300,300), ft_rng=(-20,20), bins=50, Title='Vertex X Residual', \
                    xaxis='Vertex X - Decay Location X [cm]', fname='residualX.png')
        vis_engine = visualization.root_Histogram(dy, rng=(-700,700), ft_rng=(-50,50), bins=50, Title='Vertex Y Residual', \
                    xaxis='Vertex Y - Decay Location Y [cm]', fname='residualY.png')
        vis_engine = visualization.root_Histogram(dz, rng=(-700,700), ft_rng=(-50,50), bins=50, Title='Vertex Z Residual', \
                    xaxis='Vertex Z - Decay Location Z [cm]', fname='residualZ.png')
#        vis_engine = visualization.root_Histogram(dt, rng=(-max,max), Title='Vertex T Resolution', \
        vis_engine = visualization.root_Histogram(dt, rng=(-100,100), ft_rng=(-3,3), bins=100, Title='Vertex T Residual', \
                    xaxis='Vertex T - Decay Time T [ns]', fname='residualT.png')
        '''

        if id_str == '0' or id_str == '2':
            scale = 10

        elif id_str == '4' or id_str == '1':
            scale = 3

        else:
            scale = 5

        x_rng = np.std(dx)/scale
        y_rng = np.std(dy)/scale
        z_rng = np.std(dz)/scale
        t_rng = np.std(dt)/scale

        pull_x = np.abs(np.mean(dx) - dx) / (x_rng * scale)
        pull_x[pull_x > 2] = 0
        pull_y = np.abs(np.mean(dy) - dy) / (y_rng * scale)
        pull_y[pull_y > 2] = 0
        pull_z = np.abs(np.mean(dz) - dz) / (z_rng * scale)
        pull_z[pull_z > 2] = 0
        pull_t = np.abs(np.mean(dt) - dt) / (t_rng * scale)
        pull_t[pull_t > 2] = 0
        print(np.count_nonzero(pull_x))
        print(np.count_nonzero(pull_y))
        print(np.count_nonzero(pull_z))
        print(np.count_nonzero(pull_t))


    

        vis_engine = visualization.root_Histogram(dx, rng=(-x_rng*3,x_rng*3), ft_rng=(-x_rng/2,x_rng/2), bins=50, Title='Vertex X Residual', \
                    xaxis='Vertex X - Decay Location X [cm]', fname='residualX{}.png'.format(id_str))
        vis_engine = visualization.root_Histogram(dy, rng=(-y_rng*3,y_rng*3), ft_rng=(-y_rng/2,y_rng/2), bins=50, Title='Vertex Y Residual', \
                    xaxis='Vertex Y - Decay Location Y [cm]', fname='residualY{}.png'.format(id_str))
        vis_engine = visualization.root_Histogram(dz, rng=(-z_rng*3,z_rng*3), ft_rng=(-z_rng/2,z_rng/2), bins=50, Title='Vertex Z Residual', \
                    xaxis='Vertex Z - Decay Location Z [cm]', fname='residualZ{}.png'.format(id_str))
        vis_engine = visualization.root_Histogram(dt, rng=(-t_rng*3,t_rng*3), ft_rng=(-t_rng/2,t_rng/2), bins=100, Title='Vertex T Residual', \
                    xaxis='Vertex T - Decay Time T [ns]', fname='residualT{}.png'.format(id_str))

        t_perc_35_65 = (np.percentile(dt,17.5), np.percentile(dt,82.5))
        x_perc_35_65 = (np.percentile(dx,17.5), np.percentile(dx,82.5))
        y_perc_35_65 = (np.percentile(dy,17.5), np.percentile(dy,82.5))
        z_perc_35_65 = (np.percentile(dz,17.5), np.percentile(dz,82.5))
        t_perc_5_95 = (np.percentile(dt,2.5), np.percentile(dt,97.5))
        x_perc_5_95 = (np.percentile(dx,2.5), np.percentile(dx,97.5))
        y_perc_5_95 = (np.percentile(dt,2.5), np.percentile(dy,97.5))
        z_perc_5_95 = (np.percentile(dz,2.5), np.percentile(dz,97.5))
        print("t percentile 17.5 {}, 82.5 {}".format(*t_perc_35_65))
        print("x percentile 17.5 {}, 82.5 {}".format(*x_perc_35_65))
        print("y percentile 17.5 {}, 82.5 {}".format(*y_perc_35_65))
        print("z percentile 17.5 {}, 82.5 {}".format(*z_perc_35_65))
        print("t percentile 2.5 {}, 97.5 {}".format(*t_perc_5_95))
        print("x percentile 2.5 {}, 97.5 {}".format(*x_perc_5_95))
        print("y percentile 2.5 {}, 97.5 {}".format(*y_perc_5_95))
        print("z percentile 2.5 {}, 97.5 {}".format(*z_perc_5_95))

        long_rng = np.std(long_res)/scale
        trans_rng = np.std(trans_res)/scale

        vis_engine = visualization.root_Histogram(long_res, rng=(-long_rng*3,long_rng*3), ft_rng=(-long_rng/2,long_rng/2), bins=50, Title='Vertex Longitudinal Residual', \
                    xaxis='Longitudinal Comp. of Res. [cm]', fname='residualLong{}.png'.format(id_str))
        vis_engine = visualization.root_Histogram(trans_res, rng=(-trans_rng*3,trans_rng*3), ft_rng=(-trans_rng/2,trans_rng/2), bins=100, Title='Vertex Transverse Residual', \
                    xaxis='Transverse Comp. of Res. [cm]', fname='residualTrans{}.png'.format(id_str))

        long_perc_35_65 = (np.percentile(long_res,17.5), np.percentile(long_res,82.5))
        trans_perc_35_65 = (np.percentile(trans_res,17.5), np.percentile(trans_res,82.5))
        long_perc_5_95 = (np.percentile(long_res,2.5), np.percentile(long_res,97.5))
        trans_perc_5_95 = (np.percentile(trans_res,2.5), np.percentile(trans_res,97.5))
        print("long percentile 17.5 {}, 82.5 {}".format(*long_perc_35_65))
        print("trans percentile 17.5 {}, 82.5 {}".format(*trans_perc_35_65))
        print("long percentile 2.5 {}, 97.5 {}".format(*long_perc_5_95))
        print("trans percentile 2.5 {}, 97.5 {}".format(*trans_perc_5_95))


    def PlotMomentum(self, num=100, cut=0):
        ''' plot momentum distribution of sim hits '''

        self.Tree.SetBranchStatus("Hit_particlePx",1)
        self.Tree.SetBranchStatus("Hit_particlePy",1)
        self.Tree.SetBranchStatus("Hit_particlePz",1)
        self.Tree.SetBranchStatus("Hit_particlePdgId", 1)

            # plot momentum distribution
        px,py,pz = [], [], []

        for ev in range(num):
            self.Tree.GetEntry(ev)
            #print(np.array(self.Tree.Hit_particlePx))
            px.extend(np.array(self.Tree.Hit_particlePx))
            py.extend(np.array(self.Tree.Hit_particlePy))
            pz.extend(np.array(self.Tree.Hit_particlePz))
#                print(px)
#                print(np.where(px == -1))
#                for pi in [px, py, pz]:
#                    np.delete(pi, np.where(pi == -1))
        px = np.array(px)
        py = np.array(py)
        pz = np.array(pz)
        p = np.sqrt(px**2 + py**2 + pz**2)

        cap = 1e5

        p[p < cut] = 0
        p[p > cap] = 0

        p = np.delete(p,np.where(p==0))

        visualization.Histogram(p, Title="Sim Momenta (cut at {} MeV)".format(cut), xaxis="Momentum [MeV]")


    def PlotBeta(self):
        ''' plot vertex velocity best estimates (only if kalman vertexer is used)'''

        self.Tree.SetBranchStatus("vertex_vx_m", 1)
        self.Tree.SetBranchStatus("vertex_vy_m", 1)
        self.Tree.SetBranchStatus("vertex_vz_m", 1)

        betas = []

        for ev in range(self.Tree.GetEntries()): # event
            self.Tree.GetEntry(ev)

            vx = util.unzip(self.Tree.vertex_vx_m)
            vy = util.unzip(self.Tree.vertex_vy_m)
            vz = util.unzip(self.Tree.vertex_vz_m)

            for v in range(len(vx)): #vertex
                for t in range(len(vx[v])): #track
                    betas.append(np.sqrt(vx[v][t]**2 + vy[v][t]**2 + vz[v][t]**2) / physics.c)

        visualization.Histogram(betas, Title="Vertex Velocity Best Estimates", xaxis="beta", log=True, rng=(0,10))


    def PlotExpectedPositions(self):
        ''' plot expected positions of tracks in made vertices on floor layers if no hits are in the floor'''

        self.Tree.SetBranchStatus("Digi_y", 1)

        self.Tree.SetBranchStatus("Track_k_m_x0", 1)
        self.Tree.SetBranchStatus("Track_k_m_y0", 1)
        self.Tree.SetBranchStatus("Track_k_m_z0", 1)

        self.Tree.SetBranchStatus("Track_k_m_velX", 1)
        self.Tree.SetBranchStatus("Track_k_m_velY", 1)
        self.Tree.SetBranchStatus("Track_k_m_velZ", 1)

        self.Tree.SetBranchStatus("Vertex_k_m_trackIndices",1)

        expected_x = []
        expected_z = []

        floor_y = np.sum(Detector.LayerYLims[2]) / 2

        hits_in_floor = False

        print('starting now')

        for ev in range(self.Tree.GetEntries()): # event
            self.Tree.GetEntry(ev)

            if ev % 100 == 0:
                print(ev)

            for digi in self.Tree.Digi_y:
                if self.Tree.Digi_y in np.mean(Detector.LayerYLims[:3],axis=1):
                    hits_in_floor = True

            vertex_trackIndices = self.Tree.Vertex_k_m_trackIndices
            vertex_trackIndices = util.unzip(vertex_trackIndices)

            #print(vertex_trackIndices)

            if len(vertex_trackIndices) > 0 and not hits_in_floor:
                for vert in range(len(vertex_trackIndices)): # vertex
                    vert = int(vert)

                    for ind in vertex_trackIndices[vert]: # track
                        ind = int(ind)

                        x = self.Tree.Track_k_m_x0[ind]
                        y = self.Tree.Track_k_m_y0[ind]
                        z = self.Tree.Track_k_m_z0[ind]

                        vx = self.Tree.Track_k_m_velX[ind]
                        vy = self.Tree.Track_k_m_velY[ind]
                        vz = self.Tree.Track_k_m_velZ[ind]

                        del_t = (floor_y - y) / vy

                        expected_x.append(x + vx * del_t)
                        expected_z.append(z + vz * del_t)

            hits_in_floor = False


#        visualization.root_Histogram(expected_x,rng=[Detector.BoxLimits[0][0]-2000,Detector.BoxLimits[0][1]+2000],Title='Expected x Positions',xaxis='x expected on {} [cm]'.format(floor_y),fname='Exp_x.png')
#        visualization.root_Histogram(expected_z,rng=[Detector.BoxLimits[2][0]-5000,Detector.BoxLimits[2][1]+5000],Title='Expected z Positions',xaxis='z expected on {} [cm]'.format(floor_y),fname='Exp_z.png')

        canv = root.TCanvas("canv","newCanvas")

        xbins, zbins = 100, 100

        xlims = (np.amin(expected_x),np.max(expected_x))
        zlims = (np.amin(expected_z),np.max(expected_z))

        #xlims = [Detector.BoxLimits[0][0]-2000,Detector.BoxLimits[0][1]+2000]
        #zlims = [Detector.BoxLimits[2][0]-5000,Detector.BoxLimits[2][1]+5000]

        Title = "Expected Hit Positions"

        hist = root.TH2F("hist",Title,xbins,xlims[0],xlims[1],zbins,zlims[0],zlims[1])

        hist.SetStats(0)
        hist.GetXaxis().SetTitle('x expected on {} [cm]'.format(floor_y))
        hist.GetYaxis().SetTitle('z expected on {} [cm]'.format(floor_y))

        for i in range(len(expected_x)):
                hist.Fill(expected_x[i],expected_z[i])

        hist.Draw("colz")

        canv.Update()
        canv.Draw()

        canv.SaveAs("Expected.png")

        visualization.root_Histogram(expected_x,rng=xlims,Title=Title,fname='Exp_x.png')
        visualization.root_Histogram(expected_z,rng=zlims,Title=Title,fname='Exp_z.png')


    def PlotVertexBeta(self):
        ''' Plot beta required for a particle to move from bottom of a track to location of the vertex it belongs to'''

        self.Tree.SetBranchStatus("Vertex_x",1)
        self.Tree.SetBranchStatus("Vertex_y",1)
        self.Tree.SetBranchStatus("Vertex_z",1)
        self.Tree.SetBranchStatus("Vertex_t",1)

        self.Tree.SetBranchStatus("Track_x0", 1)
        self.Tree.SetBranchStatus("Track_y0", 1)
        self.Tree.SetBranchStatus("Track_z0", 1)
        self.Tree.SetBranchStatus("Track_t0", 1)

        self.Tree.SetBranchStatus("Vertex_trackIndices",1)

        betas = []

        for ev in range(self.Tree.GetEntries()): # event
            self.Tree.GetEntry(ev)

            vertex_trackIndices = self.Tree.Vertex_trackIndices
            vertex_trackIndices = util.unzip(vertex_trackIndices)

            #print(vertex_trackIndices)

            if len(vertex_trackIndices) > 0:
                for vert in range(len(vertex_trackIndices)): # vertex
                    vert = int(vert)

                    vx = self.Tree.Vertex_x[vert]
                    vy = self.Tree.Vertex_y[vert]
                    vz = self.Tree.Vertex_z[vert]
                    vt = self.Tree.Vertex_t[vert]

                    for ind in vertex_trackIndices[vert]: # track
                        x = self.Tree.Track_x0[ind]
                        y = self.Tree.Track_y0[ind]
                        z = self.Tree.Track_z0[ind]
                        t = self.Tree.Track_t0[ind]

                        beta = np.sqrt(((x-vx)**2 + (y-vy)**2 + (z-vz)**2) / (physics.c**2 * (t-vt)**2))

                        betas.append(beta)

        visualization.Histogram(betas)


    def TrackEfficiency(self):
        ''' plot distribution of reconstructed lowest hit location'''

        self.Tree.SetBranchStatus("Track_k_m_x0", 1)
        self.Tree.SetBranchStatus("Track_k_m_y0", 1)
        self.Tree.SetBranchStatus("Track_k_m_z0", 1)
        self.Tree.SetBranchStatus("NumTracks_k_m", 1)

        self.Tree.SetBranchStatus("GenParticle_y",1)
        self.Tree.SetBranchStatus("GenParticle_px",1)
        self.Tree.SetBranchStatus("GenParticle_py",1)
        self.Tree.SetBranchStatus("GenParticle_pz",1)
        self.Tree.SetBranchStatus("GenParticle_energy",1)

        xT, zT = [], []
        xG, zG = [], []

        has_a_track = 0

        y_layer = 8500 # [cm]

        for ev in range(self.Tree.GetEntries()): # event
            self.Tree.GetEntry(ev)

            for tr in range(self.Tree.NumTracks_k_m):

                xT.append(self.Tree.Track_k_m_x0[tr])
                zT.append(self.Tree.Track_k_m_z0[tr])
            '''
            for pt in range(len(GenParticle_px)):

                xG.append()
                zG.append()
            '''
            if self.Tree.NumTracks_k_m > 0:
                has_a_track += 1

        #plt.hist2d(trackx0, trackz0)
        '''
        hist, bin_x, bin_y = np.histogram2d(xT,zT,bins=150) # change this to ROOT 2D

        fig, ax = plt.subplots(figsize=(8,5))
        plt.imshow(hist,alpha=0.5)
        plt.colorbar(orientation='vertical')
        plt.savefig("Efficiency.png")
        '''
        canv = root.TCanvas("canv","newCanvas")

        xbins, zbins = 100, 100

        #zlims = (np.amin(zT), np.max(zT))
        #xlims = (np.amin(xT),np.max(xT))
        zlims = (1.5e4, 1e4)
        xlims = (-2000,3000)

        Title = "Lowest Hit Positions for made tracks"

        hist = root.TH2F("hist",Title,xbins,xlims[0],xlims[1],zbins,zlims[0],zlims[1])

        hist.SetStats(0)
        hist.GetXaxis().SetTitle("X [cm]")
        hist.GetYaxis().SetTitle("Z [cm]")

        for i in range(len(xT)):
            hist.Fill(xT[i],zT[i])

        hist.Draw("colz")

        canv.Update()
        canv.Draw()

        canv.SaveAs("Efficiency.png")

        efficiency = has_a_track / self.Tree.GetEntries()

        print("efficiency is ", efficiency)


    def TrackResolution(self):
        ''' plot distribution of truth vs reconstructed for angle and lowest hit location '''

        self.Tree.SetBranchStatus("Track_k_m_x0", 1)
        self.Tree.SetBranchStatus("Track_k_m_y0", 1)
        self.Tree.SetBranchStatus("Track_k_m_z0", 1)

        self.Tree.SetBranchStatus("Track_k_velX", 1)
        self.Tree.SetBranchStatus("Track_k_velY", 1) # should this be _k_m_?
        self.Tree.SetBranchStatus("Track_k_velZ", 1)
        self.Tree.SetBranchStatus("NumTracks_k_m", 1)

        self.Tree.SetBranchStatus("Hit_x", 1)
        self.Tree.SetBranchStatus("Hit_y", 1)
        self.Tree.SetBranchStatus("Hit_z", 1)

        self.Tree.SetBranchStatus("Hit_particlePx", 1)
        self.Tree.SetBranchStatus("Hit_particlePy", 1)
        self.Tree.SetBranchStatus("Hit_particlePz", 1)

        min_distances = []
        angle_difference = []

        for ev in range(self.Tree.GetEntries()): # event
            self.Tree.GetEntry(ev)

            for tr in range(self.Tree.NumTracks_k_m):
                min_sim_distance = 1e6
                min_sim_index = -1

                for sim in range(len(self.Tree.Hit_x)):
                    distance = (self.Tree.Track_k_m_x0[tr] - self.Tree.Hit_x[sim])**2 \
                         + (self.Tree.Track_k_m_y0[tr] - self.Tree.Hit_y[sim])**2 \
                         + (self.Tree.Track_k_m_z0[tr] - self.Tree.Hit_z[sim])**2

                    if np.sqrt(distance) < min_sim_distance:
                        min_sim_distance = np.sqrt(distance)
                        min_sim_index = sim

                min_distances.append(min_sim_distance)

                momentum = np.array([self.Tree.Hit_particlePx[min_sim_index],
                             self.Tree.Hit_particlePy[min_sim_index],
                             self.Tree.Hit_particlePz[min_sim_index]])

                track_vel = np.array([self.Tree.Track_k_velX[tr],
                              self.Tree.Track_k_velY[tr],
                              self.Tree.Track_k_velZ[tr]])

                angle = np.sum(momentum * track_vel) \
                    / np.sqrt(np.sum(track_vel**2) * np.sum(momentum**2))

                angle_difference.append(angle)

        visualization.Histogram(min_distances,Title='Truth vs Track distance',xaxis='Distance [cm]',fname='distance.png')
        visualization.Histogram(angle_difference,Title='Truth vs Track angle',xaxis='cos(angle) []',log=True,fname='angle.png')


    def GunStudies(self):

        self.Tree.SetBranchStatus("*",0)
        self.Tree.SetBranchStatus("Hit_particlePx", 1)
        self.Tree.SetBranchStatus("Hit_particlePy", 1)
        self.Tree.SetBranchStatus("Hit_particlePz", 1)
        self.Tree.SetBranchStatus("Hit_particleEnergy", 1)

        self.Tree.SetBranchStatus("Hit_x", 1)
        self.Tree.SetBranchStatus("Hit_y", 1)
        self.Tree.SetBranchStatus("Hit_z", 1)
        self.Tree.SetBranchStatus("Hit_time", 1)

        self.Tree.SetBranchStatus("Track_k_m_t0", 1)
        self.Tree.SetBranchStatus("Track_k_m_x0", 1)
        self.Tree.SetBranchStatus("Track_k_m_y0", 1)
        self.Tree.SetBranchStatus("Track_k_m_z0", 1)

        # t, x, y, z
        vertex_truth = [0, 0, 10547, 12000]

        dxs_y, dzs_y = [], []
        dxs_t, dys_t, dzs_t = [], [], []
        dang_yx, dang_yz = [], []

        for ev in range(self.Tree.GetEntries()): # event
            self.Tree.GetEntry(ev)

            print("event {}".format(ev)) if ev % 100 == 0 else None

            #if len(self.Tree.Track_k_m_y0) != 0:
            for tr in range(len(self.Tree.Track_k_m_x0)):
                t_0 = self.Tree.Track_k_m_t0[tr]
                x_0 = self.Tree.Track_k_m_x0[tr]
                y_0 = self.Tree.Track_k_m_y0[tr]
                z_0 = self.Tree.Track_k_m_z0[tr]

                dy = vertex_truth[2] - y_0

                #highest_hit_y_ind = np.where(self.Tree.Hit_y == np.max(self.Tree.Hit_y))[0]
                #highest_hit_y_ind = int(highest_hit_y_ind)

                hit_y = np.array(self.Tree.Hit_y)
                highest_hit_y_ind = np.where(hit_y == np.max(hit_y))[0][0]
                highest_hit_y_ind = int(highest_hit_y_ind)

                px = self.Tree.Hit_particlePx[highest_hit_y_ind]
                py = self.Tree.Hit_particlePy[highest_hit_y_ind]
                pz = self.Tree.Hit_particlePz[highest_hit_y_ind]

                if py != 0:
                    pred_x = vertex_truth[1] + dy * (px / py)
                    pred_z = vertex_truth[3] + dy * (pz / py)

                    dx = np.abs(x_0 - pred_x)
                    dz = np.abs(z_0 - pred_z)

                    dxs_y.append(dx)
                    dzs_y.append(dz)

                    yx_ang = np.abs(np.arctan(px / py))
                    yz_ang = np.abs(np.arctan(pz / py))

                    dang_yx.append(yx_ang)
                    dang_yz.append(yz_ang)

                dt = t_0 - vertex_truth[0]

                #p = np.sqrt(px**2 + py**2 + pz**2)
                E = self.Tree.Hit_particleEnergy[highest_hit_y_ind]

                if E != 0:
                    vx = physics.c * px / E
                    vz = physics.c * pz / E
                    vy = physics.c * py / E

                    pred_x = vertex_truth[1] + dt * vx
                    pred_y = vertex_truth[2] + dt * vy
                    pred_z = vertex_truth[3] + dt * vz

                    dx = np.abs(x_0 - pred_x)
                    dy = np.abs(y_0 - pred_y)
                    dz = np.abs(z_0 - pred_z)

                    dxs_t.append(dx)
                    dys_t.append(dy)
                    dzs_t.append(dz)

        '''
        '''
        visualization.Histogram(dxs_y, fname="dxs_y.png", Title="Distance between truth (predicted with y) and x0", xaxis="Distance [cm]")
        visualization.Histogram(dzs_y, fname="dzs_y.png", Title="Distance between truth (predicted with y) and z0", xaxis="Distance [cm]")

        visualization.Histogram(dxs_t, fname="dxs_t.png", Title="Distance between truth (predicted with t) and x0", xaxis="Distance [cm]", rng=(0,1000))
        visualization.Histogram(dys_t, fname="dys_t.png", Title="Distance between truth (predicted with t) and y0", xaxis="Distance [cm]", rng=(0,1000))
        visualization.Histogram(dzs_t, fname="dzs_t.png", Title="Distance between truth (predicted with t) and z0", xaxis="Distance [cm]", rng=(0,1000))

        visualization.Histogram(dang_yx, fname="dang_yx.png", Title="x-y plane momentum projection angle", xaxis="angle [rad]")
        visualization.Histogram(dang_yz, fname="dang_yz.png", Title="z-y plane momentum projection angle", xaxis="angle [rad]")



    def Compare_sim_and_digi(self):

        self.Tree.SetBranchStatus("Hit_x", 1)
        self.Tree.SetBranchStatus("Hit_y", 1)
        self.Tree.SetBranchStatus("Hit_z", 1)
        self.Tree.SetBranchStatus("Digi_x", 1)
        self.Tree.SetBranchStatus("Digi_y", 1)
        self.Tree.SetBranchStatus("Digi_z", 1)

        hit_x, hit_z = [], []
        digi_x, digi_z = [], []

        for ev in range(self.Tree.GetEntries()): # event
            self.Tree.GetEntry(ev)

            print("event {}".format(ev)) if ev % 100 == 0 else None

            for hit in range(len(self.Tree.Hit_x)):
                hit_x.append(self.Tree.Hit_x[hit])
                hit_z.append(self.Tree.Hit_z[hit]+9900)

            for digi in range(len(self.Tree.Digi_x)):
                digi_x.append(self.Tree.Digi_x[digi])
                digi_z.append(self.Tree.Digi_z[digi])


        visualization.root_2D_Histogram(hit_x, hit_z, [np.amin(hit_x),np.max(hit_x)], [np.amin(hit_z),np.max(hit_z)], Title='Sim x-z', xbins=100, zbins=100, xlabel='x', zlabel='z', fname='sims.png')
        visualization.root_2D_Histogram(digi_x, digi_z, [np.amin(digi_x),np.max(digi_x)], [np.amin(hit_z),np.max(hit_z)], Title='Digi x-z', xbins=100, zbins=100, xlabel='x', zlabel='z', fname='digis.png')


    def Chi_by_ndof(self):
        self.Tree.SetBranchStatus("local_chi_s",1)
        self.Tree.SetBranchStatus("Track_k_smoothchi",1)
        self.Tree.SetBranchStatus("Track_k_hitIndices",1)

        num_layers = 10

        chi_by_hits = [[] for i in range(num_layers)]
        chi_by_hits_all = [[] for i in range(num_layers)]

        for ev in range(self.Tree.GetEntries()): # event
            self.Tree.GetEntry(ev)

            chi_s = util.unzip(self.Tree.Track_k_smoothchi)
            chi_s_all = util.unzip(self.Tree.local_chi_s)

            hit_inds = util.unzip(self.Tree.Track_k_hitIndices)

            #if ev in set([i for i in range(1,20)]):
            #    print("chi sum for all tracks \n",chi_s_all)

            #print("chi smooth: ",chi_s)
            #print("hit inds ",hit_inds)

            for tr in range(len(chi_s)): # track

                #if len(hit_inds) == 0:
                #    continue

                num_hits = len(hit_inds[tr])

                #chi_by_hits[num_hits-1].append(np.sum(chi_s[tr]))
                #chi_by_hits[num_hits-1].append(chisquare(chi_s[tr], ddof = -(4 * num_hits - 6) - num_hits + 1)[1])
                chi_by_hits[num_hits-1].append(stats.chi2.cdf(np.sum(chi_s[tr]), (4 * num_hits - 6)))
                #for ht in range(len(chi_s[tr])):
                #    chi_by_hits[num_hits-1].append(stats.chi2.cdf(chi_s[tr][ht], (4 * num_hits - 6)))

            #print("chi smooth: ",chi_s_all)

            for tr in range(len(chi_s_all)): # track

                #chi_by_hits_all[len(chi_s_all[tr][:-1])].append(np.sum(chi_s_all[tr][:-1]))

                num_hits = len(chi_s_all[tr])

                #chi_by_hits_all[num_hits-1].append(stats.chi2.cdf(np.sum(chi_s_all[tr]), (4 * num_hits - 6))) # for p-value
                chi_by_hits_all[num_hits-1].append(np.sum(chi_s_all[tr])) # for chi sum

                #for ht in range(len(chi_s_all[tr])):
                #    chi_by_hits_all[num_hits-1].append(stats.chi2.cdf(chi_s_all[tr][ht], (4 * num_hits - 6)))



        for hits in range(num_layers):

            if len(chi_by_hits[hits]) == 0:
                continue

            canv = root.TCanvas("canv","newCanvas")

            #xlims = (0,1) #,200 * (hits + 1))
            xlims = (0,1) #,200 * (hits + 1))

            #Title = "Chi sum for {} hits of made tracks".format(hits+1)
            Title = "P-value for {} hits of made tracks".format(hits+1)

            #hist = root.TH1F("hist",Title,100,xlims[0],xlims[1])
            num_bins = 30
            hist = root.TH1F("hist",Title,num_bins-1,array('d',sorted([0]+[10**(-i) for i in range(num_bins-1)])))
            #hist = root.TH1F("hist",Title,num_bins-1,[10**(-i) for i in range(num_bins)])
            #hist = root.TH1F("hist",Title,10-1,np.array([10**(-i) for i in range(10)].reverse()))

            #hist.GetXaxis().SetTitle("Chi / ndof sum by track")
            hist.GetXaxis().SetTitle("Track 1 - p-values")

            for chi in chi_by_hits[hits]:
                hist.Fill(1 - chi)
                #hist.Fill(chi / (4 * (hits + 1) - 6))

            root.gStyle.SetOptStat(111111)

            canv.SetLogy()
            canv.SetLogx()

            hist.Draw()

            canv.Update()
            canv.Draw()

            #print("(Passed tracks) / (total tracks): ", len(chi_by_hits[hits])/len(chi_by_hits_all[hits]))

            canv.SaveAs("Chis_{}_hits.png".format(hits + 1))



        for hits in range(num_layers):

            if len(chi_by_hits[hits]) == 0:
                continue

            canv = root.TCanvas("canv","newCanvas")

            xlims = (0,1) #,200 * (hits + 1))

            #Title = "Chi sum for {} hits of made tracks".format(hits+1)
            Title = "P-value for {} hits of tracks at global chi cut".format(hits+1)

            #hist = root.TH1F("hist",Title,100,xlims[0],xlims[1])
            num_bins = 30
            hist = root.TH1F("hist",Title,num_bins-1,array('d',sorted([0]+[10**(-i) for i in range(num_bins-1)])))

            #hist.GetXaxis().SetTitle("Chi / ndof sum by track")
            hist.GetXaxis().SetTitle("Track 1 - p-values")

            for chi in chi_by_hits_all[hits]:
                hist.Fill(1-chi)
                #hist.Fill(chi / (4 * (hits + 1) - 6))

            root.gStyle.SetOptStat(111111)

            canv.SetLogy()
            canv.SetLogx()

            hist.Draw()

            canv.Update()
            canv.Draw()

            print("(Passed tracks) / (total tracks = {}): ".format(len(chi_by_hits_all[hits])), len(chi_by_hits[hits])/len(chi_by_hits_all[hits]))

            canv.SaveAs("Chis_all_{}_hits.png".format(hits + 1))

        '''
        '''
        for hits in range(num_layers):

            if len(chi_by_hits_all[hits]) == 0:
                continue

            canv = root.TCanvas("canv","newCanvas")

            xlims = (np.amin(chi_by_hits_all[hits])/ (4 * (hits + 1) - 6),np.max(chi_by_hits_all[hits])/ (4 * (hits + 1) - 6)) #,200 * (hits + 1))
            #xlims = (np.amin(chi_by_hits_all[hits]),np.max(chi_by_hits_all[hits])) #,200 * (hits + 1))

            #Title = "Chi sum for {} hits of made tracks".format(hits+1)
            Title = "Chi / ndof for {} hits of tracks that reach global chi cut".format(hits+1)

            hist = root.TH1F("hist",Title,100,xlims[0],xlims[1])
            #num_bins = 30

            #hist.GetXaxis().SetTitle("Chi / ndof sum by track")
            hist.GetXaxis().SetTitle("Track Chi / ndof at Global Cut")

            for chi in chi_by_hits_all[hits]:
                #hist.Fill(chi)
                hist.Fill(chi / (4 * (hits + 1) - 6))

            root.gStyle.SetOptStat(111111)

            #canv.SetLogy()
            #canv.SetLogx()

            hist.Draw()

            canv.Update()
            canv.Draw()

            #print("(Passed tracks) / (total tracks): ", len(chi_by_hits[hits])/len(chi_by_hits_all[hits]))

            canv.SaveAs("Chis_all_{}_hits_sum.png".format(hits + 1))
        '''
        [print("1 - p-value cut for {} dof: {:0.2g}".format(i,1-stats.chi2.cdf(100,(4 * i - 6)))) for i in range(3,11)]
        [print("1 - p-value cut for {} dof: {:0.2g}".format(i,1-root.Math.chisquared_cdf(100,(4 * i - 6)))) for i in range(3,11)]
        '''



    def Signal_Evaluation_Plots(self):

        passed_events = load('passed_events.joblib')

        self.Tree.SetBranchStatus("vertex_k_m_chi2",1)
        #self.Tree.SetBranchStatus("Track_k_m_chi2PerNdof",1)
        self.Tree.SetBranchStatus("Track_k_chi2PerNdof",1)
        self.Tree.SetBranchStatus("Track_k_smooth_chi_sum",1)
        vertex_chis, track_chis = [], []

        print(self.Tree_nm)

        for ev in passed_events[self.Tree_nm][-1]: # look at survivors
            ev = int(ev)
            self.Tree.GetEntry(ev)

            for vchi in self.Tree.vertex_k_m_chi2:
                vertex_chis.append(vchi)

#            for tchi in self.Tree.Track_k_chi2PerNdof:
            for tchi in self.Tree.Track_k_smooth_chi_sum:
                track_chis.append(tchi)

        visualization.root_Histogram(track_chis, fname="track_chis_{}.png".format(self.Tree_nm.split("_")[-2]),
            Title="Signal Track Chi's", xaxis="Track Chi per Ndof []")
        visualization.root_Histogram(vertex_chis, fname="vertex_chis_{}.png".format(self.Tree_nm.split("_")[-2]),
            Title="Signal Vertex Chi's", xaxis="Vertex Chi per Ndof []")

    def get_recon_kalman(self, Tree=None, EventNumber=-1):
        if Tree is None:
            Tree=self.Tree
            EventNumber=self.EventNumber
        Tree.GetEntry(EventNumber);

        event_vis={}
        event_vis["vertex_primary"]=[]
        event_vis["vertex"]=[]
        event_vis["track"]=[]
        event_vis["hit"]=[]
        event_vis["hit_unused"]=[]
        event_vis["track_nonvertex"]=[]        
        event_vis["hit_nonvertex"]=[]        


        # indices of used digi hits (organised by track and organised into a set)
        digi_hit_inds = util.unzip(Tree.Track_k_m_hitIndices)
        used_inds = set(Tree.Track_k_m_hitIndices)
        used_track_inds = set()

        # indices of tracks used in the made vertices
        vert_inds = util.unzip(Tree.Vertex_k_m_trackIndices)

        # track best estimates
        x_s = util.unzip(Tree.x_estimates_m)
        y_s = util.unzip(Tree.y_estimates_m)
        z_s = util.unzip(Tree.z_estimates_m)

        # velocity best estimates at the vertex
        #vxv = util.unzip(self.Tree.vertex_vx_m)
        #vyv = util.unzip(self.Tree.vertex_vy_m)
        #vzv = util.unzip(self.Tree.vertex_vz_m)

        # get truth vertex info (decay location of parent particle used in the event)
        try:
            used_gens_inds = np.where(np.array(Tree.GenParticle_G4index) != -1)[0]
            if len(used_gens_inds) != 0:
                # Vertex truth location (translated from Pythia coordinate system)
                vert_truth = [Tree.GenParticle_y[int(used_gens_inds[0])] / 10,
                          Tree.GenParticle_x[int(used_gens_inds[0])] / 10,
                          Tree.GenParticle_z[int(used_gens_inds[0])] / 10]

                # pass truth vertex to visualiser
                event_vis["vertex_primary"]=[vert_truth[0], vert_truth[1], vert_truth[2]]
        except:
            pass


        # Extract vertices, tracks, and hits
        for n in range(len(Tree.Vertex_k_m_x)):
            xv = Tree.Vertex_k_m_x[n]
            yv = Tree.Vertex_k_m_y[n]
            zv = Tree.Vertex_k_m_z[n]

            # save vertex
            event_vis["vertex"].append([xv, yv, zv])

            # loop over tracks in the current vertex
            for trk_ind in vert_inds[n]:
                # save track
                event_vis["track"].append([x_s[int(trk_ind)], y_s[int(trk_ind)], z_s[int(trk_ind)]])
                # save each hit used in the current track
                event_vis["hit"].append([])
                for ind in digi_hit_inds[int(trk_ind)]:
                    x = self.Tree.Digi_x[ind]
                    y = self.Tree.Digi_y[ind]
                    z = self.Tree.Digi_z[ind]
                    event_vis["hit"][-1].append([x, y, z])
                used_track_inds.add(int(trk_ind))


        # Extract hits that weren't used in a track
        for n in range(len(Tree.Digi_numHits)):
            if not (n in used_inds):
                x = Tree.Digi_x[n]
                y = Tree.Digi_y[n]
                z = Tree.Digi_z[n]
                event_vis["hit_unused"].append([x, y, z])

        # Extract tracks and hits that are not used in a vertex
        for n in range(Tree.NumTracks_k_m):
            if not (n in used_track_inds):
                event_vis["track_nonvertex"].append([x_s[n], y_s[n], z_s[n]])
                event_vis["hit_nonvertex"].append([])
                for ind in digi_hit_inds[n]:
                    x = Tree.Digi_x[ind]
                    y = Tree.Digi_y[ind]
                    z = Tree.Digi_z[ind]
                    event_vis["hit_nonvertex"][-1].append([x, y, z])

        #king_move_inds = util.unzip(self.Tree.king_move_inds)
        #print('king move indices are ',king_move_inds)
        # show hits dropped from the hit pool (via king moves algorithm)
        # for inds in king_move_inds:
        #     for ind in inds:
        #         x = self.Tree.Digi_x[ind]
        #         y = self.Tree.Digi_y[ind]
        #         z = self.Tree.Digi_z[ind]
        #         visEngine_local.AddPoint( [[x, y, z], 'm'] )

        return event_vis
    
    def get_truthtrack(self, event=None):
        if event is None:
            event=self
        event.ExtractTruthPhysics()
        tracks=[]
        for track in event.truthTrackList:
            x = [];    y = [];    z = [];  t = [];
            for point in track.pointList:
                x.append(point.location.x)
                y.append(point.location.y)
                z.append(point.location.z)
                t.append(point.time)
            tracks.append(np.array([x,y,z,t]))
        return tracks
    
    
    def get_hits_truth(self, event=None, keys_truth = ['Hit_energy', 'Hit_time', 'Hit_detId', 'Hit_particlePdgId', 'Hit_G4TrackId', 'Hit_G4ParentTrackId', 'Hit_x', 'Hit_y', 'Hit_z', 'Hit_particleEnergy', 'Hit_particlePx', 'Hit_particlePy', 'Hit_particlePz']):
        import numpy as np
        if event is None:
            event=self
        Tree = event.Tree
        Tree.GetEntry(event.EventNumber)

        n_hits = len(Tree.Hit_x)
        hits={}
        scope = locals()
        for key in keys_truth:
            exec(f'hits["{key}"]=np.array([Tree.{key}[i] for i in range(n_hits)])', scope)
        return hits   

    def get_hits_digi(self, event=None, keys_truth = ['Digi_time', 'Digi_x', 'Digi_y', 'Digi_z', 'Digi_energy', 'Digi_px', 'Digi_py', 'Digi_pz', 'Digi_particle_energy', 'Digi_pdg_id']):
        import numpy as np
        if event is None:
            event=self
        Tree = event.Tree
        Tree.GetEntry(event.EventNumber)

        n_hits = len(Tree.Digi_x)
        hits={}
        scope = locals()
        for key in keys_truth:
            exec(f'hits["{key}"]=np.array([Tree.{key}[i] for i in range(n_hits)])', scope)
        return hits    
    
    def projectMomentum(self, event=None):
        """Adds a second point one unit vector away based on the momentum"""
        if event==None:
            event=self
        tracks = []
        self.Tree.GetEntry(self.EventNumber)
        hitn=0
        x = [];    y = [];    z = [];
        x.append(self.Tree.Hit_x.at(hitn))
        y.append(self.Tree.Hit_y.at(hitn))
        z.append(self.Tree.Hit_z.at(hitn))

        norm = 1/500 * np.sqrt(self.Tree.Hit_particlePx.at(hitn)**2 + self.Tree.Hit_particlePy.at(hitn)**2 + self.Tree.Hit_particlePz.at(hitn)**2)
        if norm == 0:
            return []

        secondx = self.Tree.Hit_x.at(hitn) - self.Tree.Hit_particlePx.at(hitn)/norm
        secondy = self.Tree.Hit_y.at(hitn) - self.Tree.Hit_particlePy.at(hitn)/norm
        secondz = self.Tree.Hit_z.at(hitn) - self.Tree.Hit_particlePz.at(hitn)/norm
        x.append(secondx)
        y.append(secondy)
        z.append(secondz)
        tracks.append(np.array([x,y,z]))
        return tracks
    
    def projectInitialMomentum(self, event=None):
        """Adds a second point one unit vector away based on the momentum"""
        if event==None:
            event=self
        tracks = []
        self.Tree.GetEntry(self.EventNumber)
        hitn=0
        x = [];    y = [];    z = [];
        x.append(self.Tree.Hit_x.at(hitn))
        y.append(self.Tree.Hit_y.at(hitn))
        z.append(self.Tree.Hit_z.at(hitn))
        
        norm = 1/500 * np.sqrt(self.Tree.Hit_particlePx.at(hitn)**2 + self.Tree.Hit_particlePy.at(hitn)**2 + self.Tree.Hit_particlePz.at(hitn)**2)
        if norm == 0:
            return []

        norm = 1/500 * np.sqrt(self.Tree.GenParticle_px.at(hitn)**2 + self.Tree.GenParticle_py.at(hitn)**2 + self.Tree.GenParticle_pz.at(hitn)**2)
        if norm == 0:
            return []
        
        secondx = self.Tree.Hit_x.at(hitn) + self.Tree.GenParticle_px.at(hitn)/norm
        secondy = self.Tree.Hit_y.at(hitn) + self.Tree.GenParticle_py.at(hitn)/norm
        secondz = self.Tree.Hit_z.at(hitn) + self.Tree.GenParticle_pz.at(hitn)/norm
        x.append(secondx)
        y.append(secondy)
        z.append(secondz)
        tracks.append(np.array([x,y,z]))
        return tracks
