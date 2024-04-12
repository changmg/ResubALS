#pragma once


#include "my_abc.h"
#include "header.h"
#include "simulator.h"
// #include "espresso_api.h"


enum class LAC_TYPE {
    CONS, RAC, RESUB
};


static inline std::ostream & operator << (std::ostream & os, const LAC_TYPE lacType) {
    const std::string strs[3] = {"CONST", "RAC", "RESUB"};
    os << strs[static_cast <ll> (lacType)];
    return os;
}


class LAC {
private:
    BigInt errBigInt;
    int targId;

public:
    explicit LAC(): errBigInt(LL_MAX), targId(-1) {}
    explicit LAC(int targ_node_id): errBigInt(LL_MAX), targId(targ_node_id) {}
    virtual ~LAC() = default;
    LAC(const LAC & oth_lac) = default;
    LAC(LAC &&) = default;
    LAC & operator = (const LAC & oth_lac) = default;
    LAC & operator = (LAC &&) = default;

    inline BigInt GetErrPro() const {return errBigInt;}
    inline void SetErrPro(BigInt& err_big_int) {errBigInt = err_big_int;}
    inline void SetErrPro(int err) {errBigInt = (BigInt)(err);}
    inline int GetTargId() const {return targId;}
    inline void SetTargId(ll targ_node_id) {targId = targ_node_id;}
    inline void Print(bool isNewLine) const {std::cout << "node " << targId << ", " << "error = " << errBigInt; if (isNewLine) std::cout << std::endl;}
};
using LACPtr = std::shared_ptr<LAC>;
using LACPtrVect = std::vector<LACPtr>;


class ResubLAC: public LAC {
public:
    int sizeGain;
    IntVect divs;
    std::string sop;

    // explicit ResubLAC() = default;
    explicit ResubLAC(int targ_node_id, int gain, const IntVect& _divs, const std::string & _sop): LAC(targ_node_id), sizeGain(gain), divs(_divs), sop(_sop) {}
    ~ResubLAC() = default;
    ResubLAC(const ResubLAC & oth_lac) = default;
    ResubLAC(ResubLAC && oth_lac) = default;
    ResubLAC & operator = (const ResubLAC & oth_lac) = default;
    ResubLAC & operator = (ResubLAC && oth_lac) = default;

    bool Check();
    ll GetAddedNodeNum();

    inline IntVect GetDivIds() const {return divs;}
    inline std::string GetSop() const {return sop;}
    inline void Print(bool isNewLine = true) const {LAC::Print(false); std::cout << ", sizeGain = " << sizeGain << ", divisors = "; PrintVect(divs, ", "); auto _sop = sop; replace(_sop.begin(), _sop.end(), '\n', ';'); std::cout << "sop = " << _sop; if (isNewLine) std::cout << std::endl;}
    inline std::string GetReprStr() const {std::ostringstream oss(""); oss << GetTargId() << ","; for (auto div: divs) oss << div << ","; oss << sop; return oss.str();}
    inline int GetSizeGain() const {return sizeGain;}
};
using ResubLACPtr = std::shared_ptr<ResubLAC>;


class LACMan {
private:
    LAC_TYPE lacType;
    int nThread;
    std::vector<LACPtr> pLacs;
    std::vector<std::vector<LACPtr>> pPromisingLacs;
    std::vector<std::vector<BitVect>> pPromisingLacPatterns;

public:
    explicit LACMan(LAC_TYPE lac_type, int n_thread): lacType(lac_type), nThread(n_thread) {}
    ~LACMan() = default;
    LACMan(const LACMan &) = delete;
    LACMan(LACMan &&) = delete;
    LACMan & operator = (const LACMan &) = delete;
    LACMan & operator = (LACMan &&) = delete;
    void Gen012ResubLACsPro(NetMan& net, IntVect& nodeIds, unsigned seed, int maxLevelDiff, int nFrame4ResubGen, int maxCandResub);
    int CollectPromisingLACs(Simulator& accSmlt, Simulator& appSmlt, ll errorMargin);
    int CollectPromisingLACs(ll nObjs, BigInt& errorMargin);
    std::shared_ptr <LAC> GetBestLac() const;
    std::shared_ptr <LAC> GetBestResubLacConsiderRealSize(NetMan& net) const;
    void GetDivs(abc::Abc_Obj_t* pNode, int nLevDivMax, RETURN_VAR IntVect& divs);
    // void GetDivsParallelly(NetMan& net, IntVect& targIds, RETURN_VAR Int2DVect& divs4Nodes);

    inline int GetLacNum() const {return static_cast<int> (pLacs.size());}
    inline std::shared_ptr <LAC> GetLac(int i) const {return pLacs[i];}
    inline ResubLACPtr GetResubLac(int i) const {return dynamic_pointer_cast<ResubLAC>(pLacs[i]);}
    inline int GetPromisingLacSize() const {return (int)(pPromisingLacs.size());}
    inline int GetPromisingLacNum(int nodeId) const {return (int)(pPromisingLacs[nodeId].size());}
    inline LACPtr GetPromisingLac(int nodeId, int lacId) const {return pPromisingLacs[nodeId][lacId];}
    inline ResubLACPtr GetPromisingResubLac(int nodeId, int lacId) const {return std::dynamic_pointer_cast<ResubLAC>(pPromisingLacs[nodeId][lacId]);}
    inline BitVect* GetPromisingLacPattern(int nodeId, int lacId) {return &pPromisingLacPatterns[nodeId][lacId];}
};


static bool IsSopXor(std::string & sop) {
    return (sop == std::string("01 1\n10 1\n")) || (sop == std::string("01 0\n10 0\n")) 
        || (sop == std::string("10 1\n01 1\n")) || (sop == std::string("10 0\n01 0\n"))
        || (sop == std::string("00 1\n11 1\n")) || (sop == std::string("00 0\n11 0\n"))
        || (sop == std::string("11 1\n00 1\n")) || (sop == std::string("11 0\n00 0\n"));
}