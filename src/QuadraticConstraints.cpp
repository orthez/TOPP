#include "QuadraticConstraints.h"

namespace TOPP {

QuadraticConstraints::QuadraticConstraints(const std::string& constraintsstring) {
    std::vector<dReal> tmpvect;
    std::string buff;
    std::istringstream iss(constraintsstring);
    getline(iss, buff, '\n');
    discrtimestep = atof(buff.c_str());
    getline(iss, buff, '\n');
    VectorFromString(buff, vmax);
    while(iss.good()) {
        getline(iss, buff, '\n');
        VectorFromString(buff, tmpvect);
        avect.push_back(tmpvect);
        getline(iss, buff, '\n');
        VectorFromString(buff,tmpvect);
        bvect.push_back(tmpvect);
        getline(iss, buff, '\n');
        VectorFromString(buff,tmpvect);
        cvect.push_back(tmpvect);
    }

    nconstraints = int(avect.front().size());
    hasvelocitylimits =  VectorMax(vmax) > TINY;
}


void QuadraticConstraints::WriteConstraints(std::stringstream& ss){
    ss << discrtimestep << "\n";
    for(int i=0; i<int(vmax.size()); i++) {
        ss << vmax[i] << " ";
    }
    ss << "\n";
    for(int i=0; i<int(avect.size()); i++) {
        for(int j=0; j<int(avect[0].size()); j++) {
            ss << avect[i][j] << " ";
        }
        ss << "\n";
        for(int j=0; j<int(avect[0].size()); j++) {
            ss << bvect[i][j] << " ";
        }
        ss << "\n";
        for(int j=0; j<int(avect[0].size()); j++) {
            ss << cvect[i][j] << " ";
        }
        if(i<int(avect.size())-1) {
            ss << "\n";
        }
    }
}

void QuadraticConstraints::CheckInput() {
    if ((int)vmax.size() != trajectory.dimension) {
        std::ostringstream msg;
        msg << "vmax has dimension " << vmax.size()
            << " but trajectory has dimension " << trajectory.dimension << ".";
        std::cout << "[TOPP] " << msg.str() << std::endl;
        throw TOPPException(msg.str());
    }
}


void QuadraticConstraints::InterpolateDynamics(dReal s, std::vector<dReal>& a, std::vector<dReal>& b, std::vector<dReal>& c) {
    a.resize(nconstraints);
    b.resize(nconstraints);
    c.resize(nconstraints);
    BOOST_ASSERT(s>=-TINY && s<=trajectory.duration+TINY);

    if(s < 0){
        s = 0;
    }
    if(s >= trajectory.duration-TINY) {
        int n = ndiscrsteps-1;
        for(int i = 0; i < nconstraints; i++) {
            a[i]= avect[n][i];
            b[i]= bvect[n][i];
            c[i]= cvect[n][i];
        }
        return;
    }

    int n = int(s/discrtimestep);
    dReal coef = (s-n*discrtimestep);
    if( std::abs(coef) <= TINY ) {
        a = avect[n];
        b = bvect[n];
        c = cvect[n];
    }
    else {
        coef /= discrtimestep;
        for (int i = 0; i < nconstraints; i++) {
            a[i] = (1-coef)*avect[n][i] + coef*avect[n+1][i];
            b[i] = (1-coef)*bvect[n][i] + coef*bvect[n+1][i];
            c[i] = (1-coef)*cvect[n][i] + coef*cvect[n+1][i];
        }
    }
}

void QuadraticConstraints::ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector) {
    dReal delta = TINY2, s2, ap, bp, cp, slope;
    std::vector<dReal> a, b, c, a2, b2, c2;
    if(s>delta) {
        delta = -delta;
    }
    s2 = s + delta;
    InterpolateDynamics(s,a,b,c);
    InterpolateDynamics(s2,a2,b2,c2);
    dReal idelta=1/delta;
    slopesvector.resize(a.size());
    for(size_t i=0; i< a.size(); i++) {
        ap = (a2[i]-a[i])*idelta;
        bp = (b2[i]-b[i])*idelta;
        cp = (c2[i]-c[i])*idelta;
        if(std::abs((2*b[i]+ap)*sd)>TINY) {
            slope = (-bp*sd*sd-cp)/((2*b[i]+ap)*sd);
        }
        else{
            slope = 0;
        }
        slopesvector[i] = slope;
    }
}

std::pair<dReal,dReal> QuadraticConstraints::SddLimits(dReal s, dReal sd){
    dReal dtsq = integrationtimestep;
    dtsq = dtsq*dtsq;
    dReal alpha = -INF;
    dReal beta = INF;
    dReal sdsq = sd*sd;
    std::vector<dReal> a, b, c;

    dReal alpha_i, beta_i;
    InterpolateDynamics(s,a,b,c);

    for(int i=0; i<nconstraints; i++) {
        if(std::abs(a[i])<TINY) {
            if(b[i]*sdsq+c[i]>0) {
                // Constraint not satisfied
                beta = -INF;
                alpha = INF;
            }
            continue;
        }
        if(a[i]>0) {
            beta_i = (-sdsq*b[i]-c[i])/a[i];
            beta = std::min(beta,beta_i);
        }
        else{
            alpha_i = (-sdsq*b[i]-c[i])/a[i];
            alpha = std::max(alpha,alpha_i);
        }
    }
    std::pair<dReal,dReal> result(alpha,beta);
    return result;
}


dReal QuadraticConstraints::SdLimitBobrowInit(dReal s){
    std::vector<dReal> a, b, c;
    InterpolateDynamics(s,a,b,c);
    if(VectorNorm(a)<TINY) {
        if(s<1e-2) {
            s+=1e-3;
        }
        else{
            s-=1e-3;
        }
        InterpolateDynamics(s,a,b,c);
    }
    std::pair<dReal,dReal> sddlimits = SddLimits(s,0);
    if(sddlimits.first > sddlimits.second) {
        return 0;
    }

    dReal sdmin = INF;
    for(int k=0; k<nconstraints; k++) {
        for(int m=k+1; m<nconstraints; m++) {
            dReal num, denum, r;
            // If we have a pair of alpha and beta bounds, then determine the sdot for which the two bounds are equal
            if(a[k]*a[m]<0) {
                num = a[k]*c[m]-a[m]*c[k];
                denum = -a[k]*b[m]+a[m]*b[k];
                if(std::abs(denum)>TINY) {
                    r = num/denum;
                    if(r>=0) {
                        sdmin = std::min(sdmin,sqrt(r));
                    }
                }
            }
        }
    }
    return sdmin;
}

// Compute the SdLimitBobrow after removing one inequality (sdot^\dag in Pham 2014)
dReal QuadraticConstraints::SdLimitBobrowExclude(dReal s, int iexclude){
    std::vector<dReal> a, b, c;
    InterpolateDynamics(s,a,b,c);
    if(VectorNorm(a)<TINY) {
        if(s<1e-2) {
            s+=1e-3;
        }
        else{
            s-=1e-3;
        }
        InterpolateDynamics(s,a,b,c);
    }
    std::pair<dReal,dReal> sddlimits = SddLimits(s,0);
    if(sddlimits.first > sddlimits.second) {
        return 0;
    }

    dReal sdmin = INF;
    for(int k=0; k<nconstraints; k++) {
        for(int m=k+1; m<nconstraints; m++) {
            if (k==iexclude || m==iexclude) {
                continue;
            }
            dReal num, denum, r;
            // If we have a pair of alpha and beta bounds, then determine the sdot for which the two bounds are equal
            if(a[k]*a[m]<0) {
                num = a[k]*c[m]-a[m]*c[k];
                denum = -a[k]*b[m]+a[m]*b[k];
                if(std::abs(denum)>TINY) {
                    r = num/denum;
                    if(r>=0) {
                        sdmin = std::min(sdmin,sqrt(r));
                    }
                }
            }
        }
    }
    return sdmin;
}


void QuadraticConstraints::FindSingularSwitchPoints() {
    if(ndiscrsteps<3) {
        return;
    }
    int i = 0;
    std::vector<dReal> a,b,c, aprev, bprev, cprev;

    InterpolateDynamics(discrsvect[i],aprev,bprev,cprev);

    for(int i=1; i<ndiscrsteps-1; i++) {
        InterpolateDynamics(discrsvect[i],a,b,c);
        dReal minsd = INF;
        dReal mins = INF;
        bool found = false;
        for(int j=0; j< int(a.size()); j++) {
            if(a[j]*aprev[j]<=0) {
                dReal adiff = a[j] - aprev[j];
                dReal ccur=c[j], bcur=b[j], scur=discrsvect[i];
                if( fabs(adiff) > TINY ) {
                    // compute the zero-crossing and linearly interpolate dynamics
                    dReal interp=-aprev[j]/adiff;
                    scur = discrsvect[i-1] + interp*(discrsvect[i]-discrsvect[i-1]);
                    bcur = bprev[j] + interp*(b[j]-bprev[j]);
                    ccur = cprev[j] + interp*(c[j]-cprev[j]);
                }

                dReal f = ccur/bcur;
                if(f<0) {
                    dReal sdstar = sqrt(-f);
                    dReal sdplus = SdLimitBobrowExclude(scur,j);
                    if(sdplus >0 && sdplus < sdstar) {
                        continue;
                    }
                    if( !found || sdstar < minsd ) {
                        found = true;
                        minsd = sdstar;
                        mins = scur;
                    }
                }
            }
        }
        if(found) {
            AddSwitchPoint2(mins, minsd, SP_SINGULAR);
        }
        aprev.swap(a);
        bprev.swap(b);
        cprev.swap(c);
    }
}

}//namespace TOPP

