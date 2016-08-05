#pragma once

#include "TOPP.h"

namespace TOPP {
class Constraints {
public:

    Trajectory trajectory;

    // Tuning parameters
    dReal discrtimestep; // Time step to discretize the trajectory when computing the MVC
    dReal integrationtimestep; // Time step to integrate the profiles
    dReal reparamtimestep; // Time step to reparameterize the trajectory based on the optimal profile
    int passswitchpointnsteps; // Number of steps to integrate around a switch point
    int extrareps; // Number of reps to lower the integrationtimestep when integration fails
    dReal bisectionprecision; // Precision when determining the sd for tangent switch point
    dReal loweringcoef; // Lowerbound that multiplies sd when searching sd for tangent switch point

    // Maximum Velocity Curves
    int ndiscrsteps; // Number of discretization steps, around trajectory.duration/discrtimestep
    std::vector<dReal> discrsvect; // Discretization points on the s axis
    std::vector<dReal> mvcbobrow;
    std::vector<dReal> mvccombined;
    bool hasvelocitylimits;
    std::vector<dReal> vmax;

    std::list<SwitchPoint> switchpointslist; // list of switch points, ordered by s
    std::list<std::pair<dReal,dReal> > zlajpahlist; // list of zlajpah points
    std::list<Profile> resprofileslist; // resulting profiles

    dReal resduration;

    dReal stepthresh; /// threshold for amount of s to step around a singularity. larger values stabilizes the graph, but can make trajectory slower.

    //////////////////////// General ///////////////////////////

    Constraints(){
        _busingcache = false;
        stepthresh = 0.01;
    }

    // Check input after this->trajectory has been set (from TOPPbindings)
    virtual void CheckInput() {
    }

    // Compute the MVCs and the switchpoints and other initializations
    virtual bool Preprocess();

    // Discretize the time interval
    virtual void Discretize();

    // Compute the MVC given by acceleration constraints
    virtual void ComputeMVCBobrow();

    // Compute the combined MVC (incorporating pure velocity constraints)
    virtual void ComputeMVCCombined();

    // Write the MVC to stringstreams
    virtual void WriteMVCBobrow(std::stringstream& ss, dReal dt=0.01);
    virtual void WriteMVCDirect(std::stringstream& ss, dReal dt=0.01);
    virtual void WriteExtra(std::stringstream& ss){
        return;
    }
    virtual void WriteConstraints(std::stringstream& ss){
        return;
    }

    // Linear interpolation
    virtual dReal Interpolate1D(dReal s, const std::vector<dReal>& v);


    //////////////////////// Limits ///////////////////////////

    // Upper limit on sd given by acceleration constraints (Bobrow)
    virtual dReal SdLimitBobrow(dReal s);
    dReal SdLimitBobrowExclude(dReal s, int iexclude){
        return BOBROWEXCLUDENOTDEFINED;
    }

    // Compute the maximum velocity curve due to dynamics at s
    // Called at initialization
    virtual dReal SdLimitBobrowInit(dReal s){
        std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    }

    // Upper limit on sd after incorporating pure velocity constraints
    virtual dReal SdLimitCombined(dReal s);
    virtual dReal SdLimitCombinedInit(dReal s);

    // Pair of (lower,upper) limits on sdd
    virtual std::pair<dReal,dReal> SddLimits(dReal s, dReal sd){
        std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    }
    virtual dReal SddLimitAlpha(dReal s, dReal sd){
        return SddLimits(s, sd).first;
    }
    virtual dReal SddLimitBeta(dReal s, dReal sd){
        return SddLimits(s, sd).second;
    }


    ///////////////////////// Switch Points ///////////////////////

    int ntangenttreated; //number of tangent switchpoints that are treated (for stats purpose)
    int nsingulartreated;  //number of singular switchpoints that are treated (for stats purpose)

    // Find all switch points, add them to switchpointslist
    virtual void FindSwitchPoints();
    virtual void FindTangentSwitchPoints();
    virtual void FindDiscontinuousSwitchPoints();

    // Switch points that are close to each other will be replaced by a single swtich point
    // Modifies switchpointslist
    virtual void TrimSwitchPoints();

    // Compute the slope of the profiles near a dynamic singularity
    virtual void ComputeSlopeDynamicSingularity(dReal s, dReal sd, std::vector<dReal>& slopesvector){
        std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    }


    // Find all the singular switch points
    // Add them to switchpointslist
    virtual void FindSingularSwitchPoints(){
        std::cout << "Virtual method not implemented\n";
        throw TOPPException("Virtual method not implemented");
    };

    // critical point of a non-valid velocity profile
    int critical_point;
    dReal critical_point_value;
    
    void SetCriticalPoint(int s, dReal sval){
            this->critical_point = s;
            this->critical_point_value = sval;
    }
    int GetCriticalPoint(){
            return this->critical_point;
    }
    dReal GetCriticalPointValue(){
            return this->critical_point_value;
    }

    // Add a switch point to switchpointslist
    virtual void AddSwitchPoint(int i, int switchpointtype, dReal sd = -1);
    virtual void AddSwitchPoint2(dReal s, dReal sd, int switchpointtype);

    std::vector<dReal> _svectcache, _sdvectcache, _sddvectcache; ///< cache
    bool _busingcache;
};
}

