// -*- coding: utf-8 -*-
// Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
//
// This file is part of the Time-Optimal Path Parameterization (TOPP) library.
// TOPP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option, any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <algorithm>

#include <vector>
#include <list>
#include <cmath>
//#include <armadillo>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "Errors.h"
#include "Trajectory.h"

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#ifdef WITH_OPENRAVE
#include <openrave/openrave.h>
#endif

#ifdef _MSC_VER
#include <boost/typeof/std/string.hpp>
#include <boost/typeof/std/vector.hpp>
#include <boost/typeof/std/list.hpp>
#include <boost/typeof/std/map.hpp>
#include <boost/typeof/std/string.hpp>
#include <boost/typeof/typeof.hpp>

#define FOREACH(it, v) for(BOOST_TYPEOF(v) ::iterator it = (v).begin(); it != (v).end(); (it)++)
#define FOREACH_NOINC(it, v) for(BOOST_TYPEOF(v) ::iterator it = (v).begin(); it != (v).end(); )

#define FOREACHC(it, v) for(BOOST_TYPEOF(v) ::const_iterator it = (v).begin(); it != (v).end(); (it)++)
#define FOREACHC_NOINC(it, v) for(BOOST_TYPEOF(v) ::const_iterator it = (v).begin(); it != (v).end(); )

#else

#if __cplusplus > 199711L || defined(__GXX_EXPERIMENTAL_CXX0X__)
#define FOREACH(it, v) for(decltype((v).begin()) it = (v).begin(); it != (v).end(); (it)++)
#define FOREACH_NOINC(it, v) for(decltype((v).begin()) it = (v).begin(); it != (v).end(); )
#define FOREACHC FOREACH
#define FOREACHC_NOINC FOREACH_NOINC
#else
#define FOREACH(it, v) for(typeof((v).begin())it = (v).begin(); it != (v).end(); (it)++)
#define FOREACH_NOINC(it, v) for(typeof((v).begin())it = (v).begin(); it != (v).end(); )
#define FOREACHC FOREACH
#define FOREACHC_NOINC FOREACH_NOINC
#endif

#endif

namespace TOPP {

#ifdef WITH_OPENRAVE
typedef OpenRAVE::dReal dReal;
#else
typedef double dReal;
#endif

#define TINY 1e-7
#define TINY2 1e-5
#define INF 1.0e15
#define MAXSD 200
#define BOBROWEXCLUDENOTDEFINED -1

/// \brief Exception that all OpenRAVE internal methods throw; the error codes are held in \ref OpenRAVEErrorCode.
class TOPPException : public std::exception
{
public:
    TOPPException() : std::exception(), _s("unknown exception"), _errorcode(0) {
    }
    TOPPException(const std::string& s, int errorcode=0) : std::exception() {
        _errorcode = errorcode;
        _s = "topp (";
        _s += boost::lexical_cast<std::string>(_errorcode);
        _s += "): ";
        _s += s;
    }
    virtual ~TOPPException() throw() {
    }
    char const* what() const throw() {
        return _s.c_str();
    }
    const std::string& message() const {
        return _s;
    }
    int GetCode() const {
        return _errorcode;
    }
private:
    std::string _s;
    int _errorcode;
};

#define TOPP_EXCEPTION_FORMAT0(s, errorcode) TOPP::TOPPException(boost::str(boost::format("[%s:%d] " s)%(__PRETTY_FUNCTION__)%(__LINE__)),errorcode)

/// adds the function name and line number to an TOPP exception
#define TOPP_EXCEPTION_FORMAT(s, args,errorcode) TOPP::TOPPException(boost::str(boost::format("[%s:%d] " s)%(__PRETTY_FUNCTION__)%(__LINE__)%args),errorcode)

////////////////////////////////////////////////////////////////////
/////////////////////////// Switch Point ///////////////////////////
////////////////////////////////////////////////////////////////////

enum SwitchPointType {
    SP_TANGENT = 0,
    SP_SINGULAR = 1,
    SP_DISCONTINUOUS = 2,
    SP_ZLAJPAH = 3
};

class SwitchPoint{
public:
    SwitchPoint(dReal s0, dReal sd0, int switchpointtype0){
        s = s0;
        sd = sd0;
        switchpointtype = switchpointtype0;
    }
    dReal s, sd; ///< the position of the singularity
    int switchpointtype; ///< the type of singularity SwitchPointType
    std::vector<dReal> slopesvector;  ///< the slope of the profiles near a dynamic singularity
};


////////////////////////////////////////////////////////////////////
/////////////////////////// Profile ////////////////////////////////
////////////////////////////////////////////////////////////////////

// \brief Velocity profile from integration
class Profile {
public:
    /// if !forward, then stores svect as the reverse of the inputs
    Profile(const std::list<dReal>&slist, const std::list<dReal>&sdlist, const std::list<dReal>&sddlist, dReal integrationtimestep, bool forward);
    /// if !forward, then stores svect as the reverse of the inputs
    Profile(const std::vector<dReal>&svect, const std::vector<dReal>&sdvect, const std::vector<dReal>&sddvect, dReal integrationtimestep, bool forward);
    Profile(){
    }

    inline void Reset() {
        svect.resize(0);
        sdvect.resize(0);
        sddvect.resize(0);
        nsteps=0;
    }
    std::vector<dReal> svect, sdvect, sddvect; ///< the points on the profile when integrating. svect is always in increasing order.
    dReal integrationtimestep; ///< the integration delta timestep
    dReal duration; ///< duration of the profile (svect.size()-1)*integrationtimestep
    int nsteps; ///< svect.size()-1
    bool forward; ///< if 1 then forward integrate in time, if 0 then backward integrate

    bool FindTimestepIndex(dReal t, int &index, dReal& remainder) const;
    // Find t such that Eval(t) = s
    // Return false if no such t
    bool Invert(dReal s,  dReal& t) const;
    dReal Eval(dReal t) const;
    dReal Evald(dReal t) const;
    dReal Evaldd(dReal t) const;
    void Print() const;
    void Write(std::stringstream& ss, dReal dt=0.001) const;

};




////////////////////////////////////////////////////////////////////
//////////////////////// Integration ///////////////////////////////
////////////////////////////////////////////////////////////////////

// Return type for the integration methods
enum IntegrationReturnType {
    INT_MAXSTEPS = 0,     // reached maximum number of steps
    INT_END = 1,          // reached s = send (FW) or s = 0 (BW)
    INT_BOTTOM = 2,       // crossed sd = 0
    INT_MVC = 3,          // crossed the MVC
    INT_PROFILE = 4       // crossed a profile
};

// Return type for the CLC computation
enum CLCReturnType {
    CLC_OK = 0,        // no problem
    CLC_SWITCH = 1,    // could not pass a switchpoint
    CLC_BOTTOM = 2     // crossed sd = 0
};

// Integrate forward from (sstart,sdstart)
/// \param[in] dt the integration timestep
/// \param[in] maxsteps the maximum steps to integrate for
int IntegrateForward(Constraints& constraints, dReal sstart, dReal sdstart, dReal dt, Profile& resprofile, int maxsteps=1e6, bool testaboveexistingprofiles=true, bool testmvc=true, bool zlajpah=false);

// Integrate backward from (sstart,sdstart)
/// \param[in] dt the integration timestep
/// \param[in] maxsteps the maximum steps to integrate for
int IntegrateBackward(Constraints& constraints, dReal sstart, dReal sdstart, dReal dt, Profile& resprofile, int maxsteps=1e6, bool testaboveexistingprofiles=true, bool testmvc=true, bool zlajpah=false);

// Compute the CLC
int ComputeLimitingCurves(Constraints& constraints);

// Compute all the profiles (CLC, forward from 0, backward from send)
int ComputeProfiles(Constraints& constraints, dReal sdbeg, dReal sdend);

// Velocity Interval Propagation
int VIP(Constraints& constraints, dReal sdbegmin, dReal sdbegmax, dReal& sdendmin, dReal& sdendmax);

// Velocity Interval Propagation from sdend backwards
int VIPBackward(Constraints& constraints, dReal& sdbegmin, dReal& sdbegmax, dReal sdendmin, dReal sdendmax);

// Emergency stopping (integrate forward the minimum acceleration)
dReal  EmergencyStop(Constraints& constraints, dReal sdbeg, Trajectory& restrajectory);

//////// ////////////////////////////////////////////////////////////
///////////////////////// Utilities ////////////////////////////////
////////////////////////////////////////////////////////////////////

// Various operations on Vectors
dReal VectorMin(const std::vector<dReal>& v);
dReal VectorMax(const std::vector<dReal>& v);
dReal VectorArgMin(const std::vector<dReal>& v);
dReal VectorArgMax(const std::vector<dReal>&v);
void VectorAdd(const std::vector<dReal>&a, const std::vector<dReal>&b,  std::vector<dReal>& res, dReal coefa=1, dReal coefb=1);
void VectorMultScalar(const std::vector<dReal>&a, std::vector<dReal>& res, dReal scalar);
dReal VectorNorm(const std::vector<dReal>&v);
void PrintVector(const std::vector<dReal>& v);

// Read a vector of dReal from a space-separated string
void VectorFromString(std::string& s,std::vector<dReal>&resvect);

/// \brief read N items from a stream and put them into vector
void ReadVectorFromStream(std::istream& s, size_t N, std::vector<dReal>& resvect);

// Solve a0 + a1*x + a2*x^2 = 0 in the interval [lowerbound,upperbound]
// Return false if no solution in [lowerbound,upperbound], true otherwise
// If more than one solution, choose the smallest
/// \brief epsilon Used to threshold the coefficients and results to avoid dividing by 0
bool SolveQuadraticEquation(dReal a0, dReal a1, dReal a2, dReal& sol, dReal lowerbound=-INF, dReal upperbound=INF, dReal epsilon=TINY);

// Check whether the point (s,sd) is above at least one profile in the list
bool IsAboveProfilesList(dReal s, dReal sd, const std::list<Profile>&resprofileslist, bool searchbackward=false, dReal softborder=TINY2);

// Find the lowest profile (in terms of sd) at s
// Return false if no profile covers s, true otherwise
bool FindLowestProfile(dReal s, Profile& profile, dReal& tres, const std::list<Profile>& resprofileslist);

/// \brief describes a location in the profile
class ProfileSample
{
public:
    ProfileSample() : s(0), sd(0), sdd(0), sindex(0), t(0) {
    }

    std::list<Profile>::const_iterator itprofile; ///< iterator into the Profile. if == constraints.resprofileslist.end(), then invalid
    dReal s, sd, sdd; ///< samples os svect, sdvect, and sddvector
    int sindex; ///< index into svect such that svect[sindex] <= s <= svect[sindex+1]
    dReal t; ///< time from svect[index] to s
};

/// \brief Find the lowest profile (in terms of sd) at s and return the iterator of the profile
///
/// If there's no profile that covers s, return resprofileslist.end()
/// \param scur the s where to look for a low sd
/// \param sdmax need to look for sd that are <= this value
ProfileSample FindLowestProfileFast(dReal scur, dReal sdmax, const std::list<Profile>& resprofileslist);

/// \brief finds if any profiles in resprofileslist intersect with ramp [s, sd, sdd] during times [0, tmax].
///
/// \param tintersect The time to the intersection point computed as sstart + tintersect*(sdstart + tintersect*0.5*sddstart)
/// \param itprofileexclude profile to exclude for intersection, usually the profile that [sstart, sdstart, sddstart] come from.
/// \return the intersection sample point
ProfileSample FindEarliestProfileIntersection(dReal sstart, dReal sdstart, dReal sddstart, dReal tmax, const std::list<Profile>& resprofileslist, std::list<Profile>::const_iterator itprofileexclude, dReal& tintersect);

void CheckInsert(std::list<std::pair<dReal,dReal> >& reslist, std::pair<dReal,dReal> e, bool inverse = false);
void FindMaxima(const std::list<std::pair<dReal,dReal> >& origlist, std::list<std::pair<dReal,dReal> >& reslist, bool inverse = false);

Profile StraightProfile(dReal sbackward,dReal sforward,dReal sdbackward,dReal sdforward);

};
