#pragma once


#include "simulator.h"
#include "my_abc.h"
#include "lac.h"


enum class METR_TYPE{
    ER, MED, MSE, MHD, SNR, MAXED, SELF
};


enum class DISTR_TYPE {
    UNIF, ENUM, MIX, SELF
};


static inline std::ostream & operator << (std::ostream & os, const METR_TYPE metrType) {
    const std::string strs[7] = {"ER", "MED", "MSE", "MHD", "SNR", "MAXED", "SELF"};
    os << strs[static_cast <ll> (metrType)];
    return os;
}


static inline std::ostream & operator << (std::ostream & os, const DISTR_TYPE distrType) {
    const std::string strs[4] = {"UNIF", "ENUM", "MIX", "SELF"};
    os << strs[static_cast <ll> (distrType)];
    return os;
}


class ErrMan {
private:
    NetMan & net0;
    NetMan & net1;
    std::shared_ptr <Simulator> pSmlt0;
    std::shared_ptr <Simulator> pSmlt1;
    unsigned seed;
    ll nFrame;
    DISTR_TYPE distrType;

public:
    ErrMan(NetMan & netMan0, NetMan & netMan1, unsigned _seed, ll n_frame, DISTR_TYPE distr_type);
    ~ErrMan() = default;
    ErrMan(const ErrMan &) = delete;
    ErrMan(ErrMan &&) = delete;
    ErrMan & operator = (const ErrMan &) = delete;
    ErrMan & operator = (ErrMan &&) = delete;

    void InitForStatErr();
    double CalcErrRate();
    double CalcMeanErrDist(bool isSign);
    double CalcMeanSquareErr(bool isSign);
    double CalcMeanHammDist();
    double CalcSigNoiseRat(bool isSign);
    ull CalcMaxErrDist(bool isSign);
    ull GET_MEM(abc::Abc_Ntk_t * pNtk1, abc::Abc_Ntk_t * pNtk2);
    ull NtkMiterComp(abc::Abc_Ntk_t * pNtk1, abc::Abc_Ntk_t * pNtk2);
    ull NtkMiterFinalize( abc::Abc_Ntk_t * pNtk1, abc::Abc_Ntk_t * pNtk2, abc::Abc_Ntk_t * pNtkMiter, ll fComb, ll nPartSize, ll fImplic, ll fMulti );
    abc::Abc_Obj_t ** X_subtract_Y_abs(abc::Abc_Ntk_t * pNtk, abc::Abc_Obj_t * X[], abc::Abc_Obj_t * Y[], ll n);
    ull GETMEM(abc::Abc_Ntk_t * pNtk, abc::Abc_Obj_t *R[], ll n);
    abc::Abc_Obj_t * X_lt_Y(abc::Abc_Ntk_t * pNtk, abc::Abc_Obj_t * X[], abc::Abc_Obj_t * Y[], ll n);
    bool SATSolver(abc::Abc_Ntk_t * pNtk);
};


// double CalcErrPro(NetMan& net0, NetMan& net1, bool isSign, unsigned seed, ll nFrame, METR_TYPE metrType, DISTR_TYPE distrType);
double CalcErr(NetMan & netMan0, NetMan & netMan1, bool isSign, unsigned seed, ll nFrame, METR_TYPE metrType, DISTR_TYPE distrType);
// double GetMSEFromSNR(NetMan & net, bool isSign, unsigned seed, ll nFrame, DISTR_TYPE distrType, double snr);


class VECBEEMan {
private:
    bool isSign;
    unsigned seed;
    int nFrame;
    METR_TYPE metrType;
    LAC_TYPE lacType;
    DISTR_TYPE distrType;
    const int nThread;
    std::vector< std::vector<BitVect> > bdPo2Nodes; // bdPo2Node[poId][nodeId], the boolean difference of poId in terms of nodeId
    std::vector< std::vector<BitVect> > bdCut2Nodes; // bdCut2Node[nodeId][cutId], the boolean difference of nodeId in terms of cutId
    std::vector<AbcObjList> disjCuts;
    std::vector<AbcObjVect> cutNtks;
    std::vector<BitVect> poMarks;
    std::vector<ll> topoIds;

public:
    VECBEEMan() = default;
    VECBEEMan(bool is_sign, unsigned _seed, int n_frame, METR_TYPE metr_type, LAC_TYPE lac_type, DISTR_TYPE distr_type, int n_thread):
        isSign(is_sign), seed(_seed), nFrame(n_frame), metrType(metr_type), lacType(lac_type), distrType(distr_type), nThread(n_thread) {}
    ~VECBEEMan() = default;
    VECBEEMan(const VECBEEMan &) = delete;
    VECBEEMan(VECBEEMan &&) = delete;
    VECBEEMan & operator = (const VECBEEMan &) = delete;
    VECBEEMan & operator = (VECBEEMan &&) = delete;

    void ComputeLocalErrorRate(Simulator& accSmlt, Simulator& appSmlt, LACMan & lacMan);
    void EstimateErrorBoundForEachLAC(Simulator& miterSmlt, LACMan& lacMan, const BigInt& upperBound, bool enableFastErrEst, const IntVect& miterId2AppId, const IntVect& appId2MiterId);
    void EstimateErrorBoundForEachLAC(Simulator& accSmlt, Simulator& appSmlt, LACMan& lacMan, const BigInt& upperBound, bool enableFastErrEst);
    void ComputeBooleanDifferenceOfPos2Nodes(Simulator& appSmlt, AbcObjVect& topoNodes, bool enableFastErrEst);
    void EstimateERBound(Simulator& accSmlt, Simulator& appSmlt, LACMan& lacMan);
    void EstimateMEDBound(Simulator& accSmlt, Simulator& appSmlt, LACMan& lacMan);
    void BatchErrEstPro(NetMan& accNet, NetMan& appNet, LACMan& lacMan, const BigInt& upperBound, bool enableFastErrEst, BigInt& backErrInt);
    void FindDisjCut(NetMan& net, AbcObjVect& topoNodes);
    void FindAppDisjCut(NetMan& net);
    // void FindAppDisjCutNew(NetMan& net, AbcObjVect& topoNodes);
    void FindDisjCutOfNode(abc::Abc_Obj_t* pObj, AbcObjList& disjCut);
    void ExpandCut(abc::Abc_Obj_t* pObj, AbcObjList& disjCut);
    abc::Abc_Obj_t* ExpandWhich(AbcObjList& disjCut);
    void CalcBoolDiffCut2Node(Simulator& appSmlt, AbcObjVect & topoNodes);
    void CalcBoolDiffCut2NodeParallelly(Simulator& appSmlt, AbcObjVect & topoNodes);
    void CalcBoolDiffPo2Node(Simulator& appSmlt, AbcObjVect& topoNodes);
    void CalcBoolDiffPo2NodeParallelly(Simulator& appSmlt, AbcObjVect& topoNodes);
    void CalcLACNonERErrs(Simulator& accSmlt, Simulator& appSmlt, LACMan& lacMan, const BigInt& upperBound);
    void CalcLACNonERErrsNew(Simulator& accSmlt, Simulator& appSmlt, LACMan& lacMan, const BigInt& upperBound, BigInt& backErrInt);
    void CalcLACERErrs(Simulator& accSmlt, Simulator& appSmlt, LACMan& lacMan);

    inline std::vector< std::vector<BitVect> >& GetBdPo2Nodes() {return bdPo2Nodes;}
};