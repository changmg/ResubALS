#pragma once


#include "header.h"
#include "my_abc.h"
#include "simulator.h"
#include "error.h"
#include "lac.h"


struct ALSOpt {
    bool isSign = false;
    bool enableFastErrEst = false;
    unsigned sourceSeed = 0;
    LAC_TYPE lacType = LAC_TYPE::RESUB;
    DISTR_TYPE distrType = DISTR_TYPE::UNIF;
    METR_TYPE metrType = METR_TYPE::ER;
    int nFrame = 102400;
    int nFrame4ResubGen = 32;
    int maxCandResub = 50000;
    int nThread = 32;
    // int maxLevelDiff = INT_MAX;
    double errUppBound = 0.05;
    std::string outpPath = "./tmp";
    std::string standCellPath = "./input/standard-cell/nangate_45nm_typ.lib";
    abc::Abc_Ntk_t * pNtk = nullptr;

    void Print();
};


class ALSMan {
private:
    bool isSign;
    bool enableFastErrEst;
    unsigned sourceSeed;
    unsigned seed;
    LAC_TYPE lacType;
    DISTR_TYPE distrType;
    METR_TYPE metrType;
    int nFrame;
    int nFrame4ResubGen;
    int maxCandResub;
    int nThread;
    int maxLevelDiff;
    int round;
    double errUppBound;
    // double errVariationTolerance;
    double maxDelay;
    NetMan accNet;
    std::string standCellPath;
    std::string outpPath;
    boost::mt19937 randGen;

    ALSMan(const ALSMan &);
    ALSMan(ALSMan &&);
    ALSMan & operator = (const ALSMan &);
    ALSMan & operator = (ALSMan &&);

public:
    explicit ALSMan(ALSOpt & opt);
    ~ALSMan() = default;
    void Run(); 
    void RunMultipleSelection();
    double ApplyTheBestLAC(NetMan & net);
    double ApplyMultipleLACs(NetMan & net);
    void FindGoodLACs(NetMan& net, const IntVect& targetNodes, LACMan& lacMan, BigInt& errorMargin, RETURN_VAR LACPtrVect& goodLACs);
    void FindGoodLACsSmallMemory(NetMan& net, const IntVect& targetNodes, LACMan& lacMan, BigInt& errorMargin, RETURN_VAR LACPtrVect& goodLACs);
    unsigned NewSeed();
    void ApplyLacPro(NetMan & net, std::shared_ptr<LAC> pLac, double backErr);
    void ApplyLacs(NetMan & net, LACPtrVect& lacs);
    void ApplyLacsConsideringErrors(NetMan& net, LACPtrVect& lacs, Simulator& accSmlt);
    void ExactSimpl(NetMan & net);
    double Eval(NetMan & net, double err, bool useYosys = false, bool isInitialCircuit = false);
    NetManPtr BuildErrorRateMiter(NetMan& accNet, NetMan& appNet, RETURN_VAR IntVect& miterId2AppId, RETURN_VAR IntVect& appId2MiterId);
    NetManPtr BuildXorOrCircuit(int nBits);
    NetManPtr BuildErrorDistanceMiter(NetMan& accNet, NetMan& appNet, RETURN_VAR IntVect& miterId2AppId, RETURN_VAR IntVect& appId2MiterId);
    NetManPtr BuildAbsoluteDifferenceCircuit(int nBits);
    NetManPtr BuildMiterWithYosys(NetMan& accNet, NetMan& appNet, RETURN_VAR IntVect& miterId2AppId, RETURN_VAR IntVect& appId2MiterId);
    NetManPtr BuildDeviationCircuit(int nBits);
    double ComputeError(Simulator& accSmlt, NetMan& net);
    double ComputeError(Simulator& accSmlt, Simulator& appSmlt);
    bool CheckError(NetMan& net, double err, int seedChange = 888);
};


Int2DVect TempApplyLacs(NetMan& net, LACPtrVect& lacs, LAC_TYPE lacType, bool isVerb = false);
void Recov(NetMan& net, Int2DVect& replTraces, bool isVerb = false);