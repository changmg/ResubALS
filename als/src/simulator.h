#pragma once


#include "header.h"
// #include "exprtkWrap.h"
#include "my_abc.h"


class Simulator: public NetMan {
private:
    unsigned seed;
    ll nFrame;
    std::vector < BitVect > dat;
    std::vector < BitVect > tempDat;

public:
    explicit Simulator(NetMan & net_man, unsigned _seed, ll n_frame);
    ~Simulator() = default;
    Simulator(const Simulator &) = delete;
    Simulator(Simulator &&) = delete;
    Simulator & operator = (const Simulator &) = delete;
    Simulator & operator = (Simulator &&) = delete;

    void InpUnif();
    void InpUnifFast();
    void InpEnum();
    void InpMix();
    void InpSelf(const std::string & fileName);
    void Sim();
    void UpdAigNode(abc::Abc_Obj_t * pObj);
    void UpdSopNode(abc::Abc_Obj_t * pObj);
    void UpdGateNode(abc::Abc_Obj_t * pObj);
    void UpdSop(abc::Abc_Obj_t * pObj, char * pSop);
    BigInt GetInp(ll iPatt, ll lsb, ll msb, bool isSign) const;
    void PrintInpStream(ll iPatt, bool isRev) const;
    BigInt GetOutp(ll iPatt) const;
    BigInt GetOutpPro(ll iPatt, bool isSign) const;
    BigInt GetTempOutpPro(ll iPatt, bool isSign) const;
    ll GetOutpFast(ll iPatt, bool isSign) const;
    ll GetTempOutpFast(ll iPatt, bool isSign) const;
    void PrintOutpStream(ll iPatt) const;
    double GetSignalProb(ll objId) const;
    void PrintSignalProb() const;
    bool IsPIOSame(const Simulator & oth_smlt) const;
    double GetErrRate(const Simulator & oth_smlt, bool isCheck = false) const;
    double GetMeanErrDist(const Simulator & oth_smlt, bool isSign, bool isCheck = false) const;
    double GetMeanSquareErr(const Simulator & oth_smlt, bool isSign, bool isCheck = false) const;
    double GetMeanHammDist(const Simulator & oth_smlt, bool isCheck = false) const;
    double GetSigNoiseRat(const Simulator & oth_smlt, bool isSign, bool isCheck = false) const;
    double GetError() const;
    void CalcLocBoolDiff(abc::Abc_Obj_t * pObj, std::list <abc::Abc_Obj_t *> & disjCut, std::vector <abc::Abc_Obj_t *> & cutNtk, std::vector < BitVect > & bdCut2Node);
    void CalcLocPartDiff(abc::Abc_Obj_t * pObj, std::list <abc::Abc_Obj_t *> & disjCut, std::vector <abc::Abc_Obj_t *> & cutNtk, std::vector < std::vector <int8_t> > & pdCut2Node);
    void UpdSopNodeForBoolAndPartDiff(abc::Abc_Obj_t * pObj);
    void UpdGateNodeForBoolAndPartDiff(abc::Abc_Obj_t * pObj);
    void UpdSopForBoolAndPartDiff(abc::Abc_Obj_t * pObj, char * pSop);

    inline ll GetFrameNumb() {return nFrame;}
    inline BitVect * GetDat(abc::Abc_Obj_t * pObj) {return &dat[GetId(pObj)];}
    inline BitVect * GetDat(ll id) {return &dat[id];}
    inline bool GetDat(abc::Abc_Obj_t * pObj, ll frameId) {return dat[GetId(pObj)][frameId];}
    inline bool GetDat(ll id, ll frameId) {return dat[id][frameId];}
    inline BitVect * GetTempDat(abc::Abc_Obj_t * pObj) {return &tempDat[GetId(pObj)];}
    inline BitVect * GetTempDat(ll id) {return &tempDat[id];}
    inline void SetTempDat(ll id, BitVect & newDat) {tempDat[id] = newDat;}
    inline void SetTempDat(ll id, BitVect && newDat) {tempDat[id] = newDat;}
    inline int CountNumbOfOnes(int id) {return dat[id].count();}
};

bool IsPIOSame(Simulator & smlt0, Simulator & smlt1);
bool IsPIOSame(NetMan & net0, NetMan & net1);
void GetNewValue(Simulator & smlt, const IntVect& faninIds, const std::string & sop, RETURN_VAR BitVect & value);