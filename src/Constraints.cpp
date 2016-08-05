#include "Constraints.h"

namespace TOPP{
bool Constraints::Preprocess() {
    switchpointslist.clear();
    resprofileslist.resize(0);
    // Change discrtimestep so as it becomes a divisor of trajectory duration
    int ndiscrsteps = int((trajectory.duration+1e-10)/discrtimestep);
    if(ndiscrsteps<1) {
        return false;
    }
    discrtimestep = trajectory.duration/ndiscrsteps;
    Discretize();

    std::chrono::time_point<std::chrono::system_clock> t0,t1,t2,t3;
    std::chrono::duration<double> d1,d2,d3;

    t0 = std::chrono::system_clock::now();
    ComputeMVCBobrow();
    t1 = std::chrono::system_clock::now();
    ComputeMVCCombined();
    t2 = std::chrono::system_clock::now();
    FindSwitchPoints();
    t3 = std::chrono::system_clock::now();

    //d1 = t1-t0;
    //d2 = t2-t1;
    //d3 = t3-t2;

    //std::cout << d1.count() <<  " " << d2.count() << " " << d3.count() << "\n";

    if(passswitchpointnsteps == 0) {
        passswitchpointnsteps = 5;
    }

    // Set integration timestep automatically if it is initially set to 0
    dReal meanmvc = 0;
    if(integrationtimestep == 0) {
        for(size_t i=0; i< mvccombined.size(); i++) {
            meanmvc += std::min(mvccombined[i],10.);
        }
        meanmvc /= mvccombined.size();
        meanmvc = std::min(1.,meanmvc);
        integrationtimestep = discrtimestep/meanmvc;
        //std::cout << "\n--------------\nIntegration timestep: " << integrationtimestep << "\n";
    }

    return true;
}


void Constraints::Discretize() {
    ndiscrsteps = int((trajectory.duration+1e-5)/discrtimestep);
    ndiscrsteps++;
    discrsvect.resize(ndiscrsteps);
    for(int i=0; i<ndiscrsteps; i++) {
        discrsvect[i] = i*discrtimestep;
    }
}


void Constraints::ComputeMVCBobrow() {
    mvcbobrow.resize(ndiscrsteps);

    for(int i=0; i<ndiscrsteps; i++) {
        mvcbobrow[i] = SdLimitBobrowInit(discrsvect[i]);
    }
}


void Constraints::ComputeMVCCombined()
{
    mvccombined.resize(ndiscrsteps);
    for(int i=0; i<ndiscrsteps; i++) {
        mvccombined[i] = SdLimitCombinedInit(discrsvect[i]);
    }
}


dReal Constraints::Interpolate1D(dReal s, const std::vector<dReal>& v) {
    assert(s>=-TINY && s<=trajectory.duration+TINY);
    if(s<0) {
        s=0;
    }
    if(s>=trajectory.duration) {
        int n = ndiscrsteps-1;
        return v[n];
    }
    int n = int(s/discrtimestep);
    dReal coef = (s-n*discrtimestep)/discrtimestep;
    return (1-coef)*v[n] + coef*v[n+1];
}


dReal Constraints::SdLimitCombinedInit(dReal s){
    dReal res = SdLimitBobrow(s);
    std::vector<dReal> qd(trajectory.dimension);
    if(hasvelocitylimits) {
        trajectory.Evald(s, qd);
        for(int i=0; i<trajectory.dimension; i++) {
            if(std::abs(qd[i])>TINY && std::abs(vmax[i])>0) {
                res = std::min(res,vmax[i]/std::abs(qd[i]));
            }
        }
    }
    return res;
}


dReal Constraints::SdLimitCombined(dReal s) {
    return Interpolate1D(s,mvccombined);
}


dReal Constraints::SdLimitBobrow(dReal s) {
    return Interpolate1D(s,mvcbobrow);
}


void Constraints::WriteMVCBobrow(std::stringstream& ss, dReal dt){
    dReal duration = trajectory.duration;
    ss << duration << " " << dt << "\n";
    for(dReal t=0; t<=duration; t+=dt) {
        ss << t << " ";
    }
    ss << "\n";
    for(dReal t=0; t<=duration; t+=dt) {
        ss << SdLimitBobrow(t) << " ";
    }
}


void Constraints::WriteMVCDirect(std::stringstream& ss, dReal dt){
    std::vector<dReal> qd(trajectory.dimension);
    dReal duration = trajectory.duration;
    ss << duration << " " << dt << "\n";
    for(dReal t=0; t<=duration; t+=dt) {
        ss << t << " ";
    }
    ss << "\n";
    for(dReal t=0; t<=duration; t+=dt) {
        dReal res = INF;
        trajectory.Evald(t, qd);
        for(int i=0; i<trajectory.dimension; i++) {
            if(std::abs(qd[i])>TINY && std::abs(vmax[i])>0) {
                res = std::min(res,vmax[i]/std::abs(qd[i]));
            }
        }
        ss << res << " ";
    }
}


void Constraints::FindSwitchPoints()
{
    switchpointslist.clear();
    FindSingularSwitchPoints();
    FindTangentSwitchPoints();
    FindDiscontinuousSwitchPoints();
    TrimSwitchPoints();
}

void Constraints::AddSwitchPoint(int i, int switchpointtype, dReal sd){
    dReal s = discrsvect[i];
    // If no sd is specified, then take the value of the mvc
    // (The case when sd is specified corresponds to a singular switchpoint in some cases)
    if(sd<0) {
        sd = mvcbobrow[i];
    }
    if(sd > MAXSD) {
        return;
    }
    std::list<SwitchPoint>::iterator it = switchpointslist.begin();
    while(it!=switchpointslist.end()) {
        if(s == it->s) {
            return;
        }
        if(s<=it->s) {
            break;
        }
        it++;
    }
    SwitchPoint sw(s,sd,switchpointtype);
    if(switchpointtype == SP_SINGULAR) {
        sw.slopesvector.resize(0);
        ComputeSlopeDynamicSingularity(s,sd,sw.slopesvector);
    }
    switchpointslist.insert(it,sw);
}

bool CompareSwitchPoint(const SwitchPoint& sw0, const SwitchPoint& sw1)
{
    return sw0.s < sw1.s;
}

void Constraints::AddSwitchPoint2(dReal s, dReal sd, int switchpointtype)
{
    // If no sd is specified, then take the value of the mvc
    // (The case when sd is specified corresponds to a singular switchpoint in some cases)
    if(sd > MAXSD) {
        return;
    }
    SwitchPoint sw(s,sd,switchpointtype);
    std::list<SwitchPoint>::iterator it = std::lower_bound(switchpointslist.begin(), switchpointslist.end(), sw, CompareSwitchPoint);
    if( it != switchpointslist.end() ) {
        if( s >= it->s+TINY ) {
            std::cout << "switch point already exists, type=" << it->switchpointtype;
            return;
        }
    }
    if(switchpointtype == SP_SINGULAR) {
        sw.slopesvector.resize(0);
        ComputeSlopeDynamicSingularity(s,sd,sw.slopesvector);
    }
    switchpointslist.insert(it,sw);
}

void Constraints::FindTangentSwitchPoints(){
    if(ndiscrsteps<3)
        return;
    int i = 1;
    dReal s,sd,snext,sdnext,alpha,diff,diffprev,tangent,prevtangent;
    std::pair<dReal,dReal> sddlimits;

    s = discrsvect[i];
    snext = discrsvect[i+1];
    sd = SdLimitBobrow(s);
    sdnext = SdLimitBobrow(snext);
    tangent = (sdnext-sd)/discrtimestep;
    prevtangent = (sd - SdLimitBobrow(discrsvect[i-1]))/discrtimestep;
    sddlimits = SddLimits(s,sd);
    alpha = sddlimits.first;
    //beta = sddlimits.second;
    diffprev = alpha/sd - tangent;

    for(int i=2; i<ndiscrsteps-1; i++) {
        s = discrsvect[i];
        snext = discrsvect[i+1];
        sd = SdLimitBobrow(s);
        sdnext = SdLimitBobrow(snext);
        sddlimits = SddLimits(s,sd);
        alpha = sddlimits.first;
        if(std::abs(prevtangent-tangent)>2 && prevtangent < 0 && tangent >0) {
            AddSwitchPoint2(s,sd,SP_DISCONTINUOUS);
        }
        prevtangent = tangent;
        tangent = (sdnext-sd)/discrtimestep;
        //if(std::abs(tangent-prevtangent)>1.) {
        //    continue;
        //}
        //beta = sddlimits.second;
        diff = alpha/sd - tangent;
        if(diffprev*diff<0 && std::abs(diff)<1) {
            AddSwitchPoint2(s,sd,SP_TANGENT);
        }
        diffprev = diff;
    }
}

void Constraints::FindDiscontinuousSwitchPoints() {
    if(ndiscrsteps<3)
        return;
    int i = 0;
    dReal sd, sdn, sdnn;
    sd = SdLimitBobrow(discrsvect[i]);
    sdn = SdLimitBobrow(discrsvect[i+1]);
    // also look for the start of the chunks for the trajectory
    //std::list<dReal>::const_iterator itchunkstart = trajectory.chunkcumulateddurationslist.begin();
    //int nLastAddedSwitchIndex = -1;
    for(int i=0; i<ndiscrsteps-2; i++) {
        sdnn = SdLimitBobrow(discrsvect[i+2]);
        if(std::abs(sdnn-sdn)>100*std::abs(sdn-sd)) {
            if(sdn<sdnn) {
                AddSwitchPoint2(discrsvect[i+1],mvcbobrow[i+1],SP_DISCONTINUOUS);
            }
            else{
                AddSwitchPoint2(discrsvect[i+2],mvcbobrow[i+2],SP_DISCONTINUOUS);
            }
        }
        // if( trajectory.degree <= 3 ) {
        //     // if the trajectory degree is <= 3, then the accelerations will not be differentiable at the trajectory chunk edges.
        //     // therefore add those discontinuity points.
        //     // perhaps there's a better way to compute this, but the above threshold doesn't catch it.
        //     if( itchunkstart != trajectory.chunkcumulateddurationslist.end() && *itchunkstart <= discrsvect[i+2]+TINY ) {
        //         if( nLastAddedSwitchIndex < i+1 ) {
        //             AddSwitchPoint2(discrsvect[i+1],mvcbobrow[i+1],SP_DISCONTINUOUS);
        //             nLastAddedSwitchIndex = i+1;
        //         }
        //         ++itchunkstart;
        //     }
        // }
        sd = sdn;
        sdn = sdnn;
    }
}


void InsertInSpList(std::list<SwitchPoint>& splist, SwitchPoint sp){
    if(splist.size()==0) {
        splist.push_back(sp);
        return;
    }
    std::list<SwitchPoint>::iterator it = splist.begin();
    while(it!=splist.end()) {
        if(sp.switchpointtype == SP_SINGULAR) {
            if((it->switchpointtype == SP_SINGULAR && it->sd >= sp.sd) ||  it->switchpointtype!=SP_SINGULAR) {
                splist.insert(it,sp);
                return;
            }
        }
        else if(it->switchpointtype!=SP_SINGULAR && it->sd >= sp.sd ) {
            splist.insert(it,sp);
            return;
        }
        it++;

    }
    splist.push_back(sp);
}


void Constraints::TrimSwitchPoints() {
    dReal radius = discrtimestep*2.1;
    dReal scur = -1, snext, sdcur = -1, sdnext;
    int stypenext;
    std::list<SwitchPoint>::iterator itcur;
    std::vector<dReal> slopesvector;
    std::list<SwitchPoint>::iterator it = switchpointslist.begin();

    // Merge singular points
    // Find the first singular point, if any
    while(it!=switchpointslist.end()) {
        if(it->switchpointtype == SP_SINGULAR) {
            scur = it->s;
            sdcur = it->sd;
            itcur = it;
            break;
        }
        it++;
    }
    // Merge consecutive singular points that are in a small radius
    if(scur>=0) {
        it++;
        while(it!=switchpointslist.end()) {
            snext = it->s;
            sdnext = it->sd;
            stypenext = it->switchpointtype;
            if(stypenext == SP_SINGULAR) {
                if(snext-scur<radius) {
                    if(sdcur < sdnext) {
                        slopesvector = it->slopesvector;
                        it =  switchpointslist.erase(it);
                    }
                    else{
                        slopesvector = itcur->slopesvector;
                        switchpointslist.erase(itcur);
                        scur = snext;
                        sdcur = sdnext;
                        itcur = it;
                        it++;
                    }
                    for(int j=0; j<int(slopesvector.size()); j++) {
                        itcur->slopesvector.push_back(slopesvector[j]);
                    }
                }
                else{
                    scur = snext;
                    sdcur = sdnext;
                    itcur = it;
                    it++;
                }
            }
            else{
                it++;
            }

        }
    }

    // Merge non-singular switchpoints
    // Find the first non singular switchpoint, if any
    it = switchpointslist.begin();
    while(it!=switchpointslist.end()) {
        if(it->switchpointtype != SP_SINGULAR) {
            scur = it->s;
            sdcur = it->sd;
            itcur = it;
            break;
        }
        it++;
    }
    // Merge consecutive non-singular switchpoints that are in a small radius
    if(scur>=0) {
        it++;
        while(it!=switchpointslist.end()) {
            snext = it->s;
            sdnext = it->sd;
            stypenext = it->switchpointtype;
            if(stypenext != SP_SINGULAR) {
                if(snext-scur<radius) {
                    if(sdcur < sdnext) {
                        slopesvector = it->slopesvector;
                        it =  switchpointslist.erase(it);
                    }
                    else{
                        slopesvector = itcur->slopesvector;
                        switchpointslist.erase(itcur);
                        // don't set scur/sdcur since then radius will be sliding
                        //scur = snext;
                        //sdcur = sdnext;
                        itcur = it;
                        it++;
                    }
                    for(int j=0; j<int(slopesvector.size()); j++) {
                        itcur->slopesvector.push_back(slopesvector[j]);
                    }
                }
                else{
                    scur = snext;
                    sdcur = sdnext;
                    itcur = it;
                    it++;
                }
            }
            else{
                it++;
            }

        }
    }

    // Sort by sd
    std::list<SwitchPoint> splist;
    it = switchpointslist.begin();
    while(it!=switchpointslist.end()) {
        InsertInSpList(splist,*it);
        it++;
    }
    switchpointslist = splist;
}

}

