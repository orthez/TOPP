#pragma once
#include "TOPP.h"
#include "Constraints.h"

namespace TOPP {

// Constraints of the form a(s)sdd + b(s)sd^2 + c(s) <= 0
// See our article http://arxiv.org/abs/1312.6533 for more details
class QuadraticConstraints : public Constraints {
public:
    QuadraticConstraints() : Constraints(){
    }
    QuadraticConstraints(const std::string& constraintsstring);

    //////////////// Overloaded methods //////////////////////
    std::pair<dReal,dReal> SddLimits(dReal s, dReal sd);
    dReal SdLimitBobrowInit(dReal s);
    dReal SdLimitBobrowExclude(dReal s, int iexclude);
    void CheckInput();
    void FindSingularSwitchPoints();
    void ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector);
    void WriteConstraints(std::stringstream& ss);

    //////////////// Specific members and methods //////////////////////
    int nconstraints;  // Number of constraints
    std::vector<std::vector<dReal> > avect, bvect, cvect;  // Dynamics coefficients. avect[i], bvect[i], cvect[i] are vectors of length 2*ndof where the first ndof are the upper limit, the next ndof are for the lower limit. These incorporate any upper/lower limits.

    void InterpolateDynamics(dReal s, std::vector<dReal>& a, std::vector<dReal>& b, std::vector<dReal>& c);   // Linearly interpolate the dynamics coefficients a,b,c

};
}
