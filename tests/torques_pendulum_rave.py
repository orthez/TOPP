# -*- coding: utf-8 -*-
# Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
#
# This file is part of the Time-Optimal Path Parameterization (TOPP) library.
# TOPP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
sys.path.append('..')

import TOPPbindings
import TOPPpy
import TOPPopenravepy
import time
import string
from pylab import *
from numpy import *
from openravepy import *


ion()

########################### Robot ################################
env = Environment() # create openrave environment
#------------------------------------------#
robotfile = "../robots/twodof.robot.xml"
env.Load(robotfile)
robot=env.GetRobots()[0]
robot.SetTransform(array([[0,0,1,0],[0,1,0,0],[-1,0,0,0.3],[0,0,0,1]]))
#------------------------------------------#
grav=[0,0,-9.8]
n=robot.GetDOF()
dof_lim=robot.GetDOFLimits()
vel_lim=robot.GetDOFVelocityLimits()
robot.SetDOFLimits(-10*ones(n),10*ones(n))
robot.SetDOFVelocityLimits(100*vel_lim)


############################ Tunings ############################
discrtimestep = 0.001
integrationtimestep = discrtimestep
reparamtimestep = 0 #auto
passswitchpointnsteps = 10
tuningsstring = "%f %f %f %d"%(discrtimestep,integrationtimestep,reparamtimestep,passswitchpointnsteps)


############################ Trajectory ############################
#------------------------------------------#

# TODO: allure mechante !
trajectorystring = """1.000000
2
-0.0950348403623 -0.0410497641469 -4.82725068944 1.82174264036
-0.165204811614 0.135763930185 0.663416639512 -0.633975758083"""

# TODO: calcul supra long et resultat chelou
trajectorystring = """1.000000
2
0.0 3.63455628188 -16.353755721 9.57760678553
0.0 0.0 0.0 0.0"""

trajectorystring = """0.534643
2
-0.496005602127 -0.491268736192 -0.705298166892
-0.879406406486 -0.871008053258 0.701512327246"""

trajectorystring = """1.000000
2
-0.0950348403547 0.0346142777958 -9.82659842058 6.74542632955
-0.165204811682 -0.546977881311 1.62609934207 -0.913916649073"""

trajectorystring = """1.000000
2
-0.0539850762113 -0.0539850762059 0.131202643008 -0.0710088455766
-0.300968741813 -0.300968741783 0.480216768411 -0.231499819896"""
sdbeg_min, sdbeg_max = 0,0



#------------------------------------------#
traj0 = TOPPpy.PiecewisePolynomialTrajectory.FromString(trajectorystring)


############################ Constraints ############################
#------------------------------------------#
vmax = array([0,0])
taumin = array([-11,-7])
taumax = array([11,7])
t0 = time.time()
constraintstring = string.join([str(x) for x in taumin]) + "\n" + string.join([str(a) for a in taumax]) + "\n" + string.join([str(a) for a in vmax])
constraintstring += "\n" + robotfile
#------------------------------------------#


############################ Run TOPP ############################
t1 = time.time()
x = TOPPbindings.TOPPInstance("TorqueLimitsRave",constraintstring,trajectorystring,tuningsstring,robot)
t2 = time.time()
#ret = x.RunComputeProfiles(0,0)
ret = x.RunVIP(sdbeg_min, sdbeg_max)
print ret
t3 = time.time()

#print x.resduration
print "sdendmin =", x.sdendmin
print "sdendmax =", x.sdendmax

# if(ret == 1):
#     x.ReparameterizeTrajectory()

t4 = time.time()

################ Plotting the MVC and the profiles #################
x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
#axis([0, 1, 0, 100])
TOPPpy.PlotAlphaBeta(x)



##################### Plotting the trajectories #####################
# if(ret == 1):
#     x.WriteResultTrajectory()
#     traj1 = TOPPpy.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
#     dtplot = 0.01
#     TOPPpy.PlotKinematics(traj0,traj1,dtplot,vmax)
#     TOPPopenravepy.PlotTorques(robot,traj0,traj1,dtplot,taumin,taumax,3)


print "\n--------------"
print "Python preprocessing: ", t1-t0
print "Building TOPP Instance: ", t2-t1
print "Compute profiles: ", t3-t2
print "Reparameterize trajectory: ", t4-t3
print "Total: ", t4-t0
if(ret == 1):
    print "Trajectory duration (estimate): ", x.resduration

raw_input()