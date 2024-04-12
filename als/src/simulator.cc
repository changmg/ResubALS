#include "simulator.h"


using namespace abc;
using namespace std;
using namespace boost;
using namespace random;


Simulator::Simulator(NetMan & net_man, unsigned _seed, ll n_frame): NetMan(net_man.GetNet(), false), seed(_seed), nFrame(n_frame) {
    auto type = NetMan::GetNetType();
    assert(type == NET_TYPE::AIG || type == NET_TYPE::GATE || type == NET_TYPE::SOP);
    dat.resize(NetMan::GetIdMaxPlus1(), BitVect (nFrame, 0));
}


void Simulator::InpUnif() {
    uniform_int <> unif01(0, 1);
    random::mt19937 eng(seed);
    variate_generator < random::mt19937, uniform_int <> > rand01(eng, unif01);

    for (ll i = 0; i < NetMan::GetPiNum(); ++i) {
        auto piId = NetMan::GetPiId(i);
        dat[piId].reset();
        for (ll j = 0; j < nFrame; ++j) {
            if (rand01())
                dat[piId].set(j);
        }
    }
    
    for (ll i = 0; i < NetMan::GetIdMaxPlus1(); ++i) {
        if (NetMan::IsConst0(i))
            dat[i].reset();
        else if (NetMan::IsConst1(i))
            dat[i].set();
    }
}


void Simulator::InpUnifFast() {
    random::uniform_int_distribution <ull> unifUll;
    random::mt19937 eng(seed);
    variate_generator <random::mt19937, random::uniform_int_distribution <ull> > randUll(eng, unifUll);
    const ll unitLength = 64;
    assert((nFrame & (unitLength - 1)) == 0);
    ll nUnit = nFrame / unitLength;

    for (ll i = 0; i < NetMan::GetPiNum(); ++i) {
        auto piId = NetMan::GetPiId(i);
        dat[piId].resize(0);
        for (ll j = 0; j < nUnit; ++j) {
            ull numb = randUll();
            dat[piId].append(numb);
        }
    }

    for (ll i = 0; i < NetMan::GetIdMaxPlus1(); ++i) {
        if (NetMan::IsConst0(i))
            dat[i].reset();
        else if (NetMan::IsConst1(i))
            dat[i].set();
    }
}


void Simulator::InpEnum() {
    #ifdef DEBUG
    assert(GetPiNum() < 30);
    assert(1ll << GetPiNum() == nFrame);
    #endif
    for (ll i = 0; i < GetPiNum(); ++i) {
        bool phase = 1;
        auto piId = GetPiId(i);
        dat[piId].reset();
        for (ll j = 0; j < nFrame; ++j) {
            if (j % (1 << i) == 0)
                phase = !phase;
            if (phase)
                dat[piId].set(j);
        }
    }

    for (ll i = 0; i < NetMan::GetIdMaxPlus1(); ++i) {
        if (NetMan::IsConst0(i))
            dat[i].reset();
        else if (NetMan::IsConst1(i))
            dat[i].set();
    }
}


void Simulator::InpMix() {
    // only for 9x8 multiplier
    const ll inWidth = 9;
    const ll wWidth = 8;
    const double inMean = 0.0;
    const double inStandDer = 85.0;

    assert(inWidth >= 1 && inWidth <= 16);
    assert(wWidth >= 1 && wWidth <= 16);
    const double inMax = (1ll << (inWidth - 1)) - 1;
    const ll wMax = (1ll << (wWidth - 1)) - 1;
    assert(GetPiNum() == inWidth + wWidth);
    assert(inMean >= -inMax && inMean <= inMax);

    // "in" follows normal distribution
    boost::random::mt19937 engine0(seed);
    boost::random::normal_distribution <double> gauss(inMean, inStandDer);
    boost::variate_generator <boost::random::mt19937, boost::random::normal_distribution <double> > randIn(engine0, gauss);

    // "w" follows uniform distribution
    boost::random::mt19937 engine1(seed);
    boost::random::uniform_int_distribution <ll> unif(-wMax, wMax);
    boost::variate_generator <boost::random::mt19937, boost::random::uniform_int_distribution <ll> > randW(engine1, unif);

    // primary inputs
    abc::Abc_Obj_t * pObj = nullptr;
    ll k = 0;
    for (ll i = 0; i < nFrame; ++i) {
        ll in = static_cast <ll> (std::min(std::max(round(randIn()), -inMax), inMax));
        ll w = randW();
        auto in2Comp = bitset <16> ((in >= 0)? in: (1ll << inWidth) + in);
        auto w2Comp = bitset <16> ((w >= 0)? w: (1ll << wWidth) + w);
        for (k = 0; k < inWidth; ++k) {
            pObj = GetPi(k);
            dat[pObj->Id][i] = in2Comp[k];
        }
        for (k = inWidth; k < inWidth + wWidth; ++k) {
            pObj = GetPi(k);
            dat[pObj->Id][i] = w2Comp[k - inWidth];
        }
        // cout << in << "," << GetInput(i, 0, inWidth - 1, true) << endl;
        // cout << w << "," << GetInput(i, inWidth, inWidth + wWidth - 1, true) << endl;
        assert(in == GetInp(i, 0, inWidth - 1, true) && w == GetInp(i, inWidth, inWidth + wWidth - 1, true));
    }

    // constant nodes
    for (ll i = 0; i < NetMan::GetIdMaxPlus1(); ++i) {
        if (NetMan::IsConst0(i))
            dat[i].reset();
        else if (NetMan::IsConst1(i))
            dat[i].set();
    }
}


void Simulator::InpSelf(const string & fileName) {
    // primary inputs
    FILE * fp = fopen(fileName.c_str(), "r");
    assert(fp != nullptr);
    const ll maxPiNumb = 1000;
    assert(GetPiNum() <= maxPiNumb);
    char buf[maxPiNumb];
    ll cnt = 0;
    while (fgets(buf, sizeof(buf), fp) != nullptr) {
        assert(static_cast <ll>(strlen(buf)) == GetPiNum() + 1);
        for (ll i = 0; i < GetPiNum(); ++i) {
            auto pObj = GetPi(i);
            dat[pObj->Id].set(cnt, buf[i] == '1');
        }
        ++cnt;
        assert(cnt <= nFrame);
    }
    assert(cnt == nFrame);
    fclose(fp);

    // constant nodes
    for (ll i = 0; i < NetMan::GetIdMaxPlus1(); ++i) {
        if (NetMan::IsConst0(i))
            dat[i].reset();
        else if (NetMan::IsConst1(i))
            dat[i].set();
    }
}


void Simulator::Sim() {
    auto type = GetNetType();
    auto nodes = TopoSort();
    for (const auto & pObj: nodes) {
        if (type == NET_TYPE::AIG)
            UpdAigNode(pObj);
        else if (type == NET_TYPE::SOP)
            UpdSopNode(pObj);
        else if (type == NET_TYPE::GATE)
            UpdGateNode(pObj);
        else
            assert(0);
    }
    for (ll i = 0; i < GetPoNum(); ++i) {
        auto pPo = GetPo(i);
        auto drivId = GetFaninId(pPo, 0);
        #ifdef DEBUG
        assert(!Abc_ObjIsComplement(pPo));
        #endif
        dat[GetId(pPo)] = dat[drivId];
    }
}


void Simulator::UpdAigNode(Abc_Obj_t * pObj) {
    #ifdef DEBUG
    assert(Abc_ObjIsNode(pObj));
    #endif
    auto pNtk = NetMan::GetNet();
    auto pMan = static_cast <Hop_Man_t *> (pNtk->pManFunc);
    auto pRoot = static_cast <Hop_Obj_t *> (pObj->pData);
    auto pRootR = Hop_Regular(pRoot);

    // skip constant node
    if (Hop_ObjIsConst1(pRootR))
        return;

    // get topological order of subnetwork in aig
    Vec_Ptr_t * vHopNodes = Hop_ManDfsNode(pMan, pRootR);

    // init internal hop nodes
    ll maxHopId = -1;
    ll i = 0;
    Hop_Obj_t * pHopObj = nullptr;
    Vec_PtrForEachEntry(Hop_Obj_t *, vHopNodes, pHopObj, i)
        maxHopId = max(maxHopId, static_cast <ll> (pHopObj->Id));
    Vec_PtrForEachEntry( Hop_Obj_t *, pMan->vPis, pHopObj, i )
        maxHopId = max(maxHopId, static_cast <ll> (pHopObj->Id));
    vector < BitVect > interData(maxHopId + 1, BitVect (nFrame, 0));
    unordered_map <ll, BitVect *> hop2Data;
    Abc_Obj_t * pFanin = nullptr;
    Abc_ObjForEachFanin(pObj, pFanin, i)
        hop2Data[Hop_ManPi(pMan, i)->Id] = &dat[pFanin->Id];

    // special case for inverter or buffer
    if (pRootR->Type == AIG_PI) {
        pFanin = Abc_ObjFanin0(pObj);
        dat[pObj->Id] = dat[pFanin->Id];
    }

    // simulate
    Vec_PtrForEachEntry(Hop_Obj_t *, vHopNodes, pHopObj, i) {
        assert(Hop_ObjIsAnd(pHopObj));
        auto pHopFanin0 = Hop_ObjFanin0(pHopObj);
        auto pHopFanin1 = Hop_ObjFanin1(pHopObj);
        #ifdef DEBUG
        assert(!Hop_ObjIsConst1(pHopFanin0));
        assert(!Hop_ObjIsConst1(pHopFanin1));
        #endif
        BitVect & data0 = Hop_ObjIsPi(pHopFanin0) ? *hop2Data[pHopFanin0->Id] : interData[pHopFanin0->Id];
        BitVect & data1 = Hop_ObjIsPi(pHopFanin1) ? *hop2Data[pHopFanin1->Id] : interData[pHopFanin1->Id];
        BitVect & out = (pHopObj == pRootR) ? dat[pObj->Id] : interData[pHopObj->Id];
        bool isFanin0C = Hop_ObjFaninC0(pHopObj);
        bool isFanin1C = Hop_ObjFaninC1(pHopObj);
        if (!isFanin0C && !isFanin1C)
            out = data0 & data1;
        else if (!isFanin0C && isFanin1C)
            out = data0 & ~data1;
        else if (isFanin0C && !isFanin1C)
            out = ~data0 & data1;
        else if (isFanin0C && isFanin1C)
            out = ~(data0 | data1);
    }

    // complement
    if (Hop_IsComplement(pRoot))
        dat[pObj->Id].flip();

    // recycle memory
    Vec_PtrFree(vHopNodes); 
}


void Simulator::UpdSopNode(Abc_Obj_t * pObj) {
    #ifdef DEBUG
    assert(Abc_ObjIsNode(pObj));
    #endif
    // skip constant node
    if (Abc_NodeIsConst(pObj))
        return;
    // update sop
    char * pSop = static_cast <char *> (pObj->pData);
    UpdSop(pObj, pSop);
}


void Simulator::UpdGateNode(Abc_Obj_t * pObj) {
    #ifdef DEBUG
    assert(Abc_ObjIsNode(pObj));
    #endif
    // skip constant node
    if (Abc_NodeIsConst(pObj))
        return;
    // update sop
    char * pSop = static_cast <char *> ((static_cast <Mio_Gate_t *> (pObj->pData))->pSop);
    UpdSop(pObj, pSop);
}


void Simulator::UpdSop(Abc_Obj_t * pObj, char * pSop) {
    ll nVars = Abc_SopGetVarNum(pSop);
    BitVect product(nFrame, 0);
    for (char * pCube = pSop; *pCube; pCube += nVars + 3) {
        bool isFirst = true;
        for (ll i = 0; pCube[i] != ' '; i++) {
            Abc_Obj_t * pFanin = Abc_ObjFanin(pObj, i);
            switch (pCube[i]) {
                case '-':
                    continue;
                    break;
                case '0':
                    if (isFirst) {
                        isFirst = false;
                        product = ~dat[pFanin->Id];
                    }
                    else
                        product &= ~dat[pFanin->Id];
                    break;
                case '1':
                    if (isFirst) {
                        isFirst = false;
                        product = dat[pFanin->Id];
                    }
                    else
                        product &= dat[pFanin->Id];
                    break;
                default:
                    assert(0);
            }
        }
        if (isFirst) {
            isFirst = false;
            product.set();
        }
        #ifdef DEBUG
        assert(!isFirst);
        #endif
        if (pCube == pSop)
            dat[pObj->Id] = product;
        else
            dat[pObj->Id] |= product;
    }

    // complement
    if (Abc_SopIsComplement(pSop))
        dat[pObj->Id].flip();
}


BigInt Simulator::GetInp(ll iPatt, ll lsb, ll msb, bool isSign) const {
    #ifdef DEBUG
    assert(lsb >= 0 && msb < NetMan::GetPiNum());
    assert(iPatt < nFrame);
    assert(lsb <= msb && msb - lsb < 512);
    #endif
    BigInt ret(0);
    for (ll k = msb; k >= lsb; --k) {
        ret <<= 1;
        if (dat[NetMan::GetPiId(k)][iPatt])
            ++ret;
    }
    if (isSign && ret >= (static_cast <BigInt> (1) << (msb - lsb)))
        ret = -((static_cast <BigInt> (1) << (msb - lsb + 1)) - ret);
    return ret;
}


void Simulator::PrintInpStream(ll iPatt, bool isRev) const {
    #ifdef DEBUG
    assert(iPatt < nFrame);
    #endif
    if (isRev) {
        for (ll k = GetPiNum() - 1; k >= 0; --k)
            cout << dat[GetPiId(k)][iPatt];
    }
    else {
        for (ll k = 0; k < GetPiNum(); ++k)
            cout << dat[GetPiId(k)][iPatt];
    }
    cout << endl;
}


BigInt Simulator::GetOutp(ll iPatt) const {
    ll lsb = 0;
    ll msb = NetMan::GetPoNum() - 1;
    #ifdef DEBUG
    assert(iPatt < nFrame);
    assert(msb < 512);
    #endif
    BigInt ret(0);
    for (ll k = msb; k >= lsb; --k) {
        ret <<= 1;
        if (dat[NetMan::GetPoId(k)][iPatt])
            ++ret;
    }
    return ret;
}


BigInt Simulator::GetOutpPro(ll iPatt, bool isSign) const {
    ll lsb = 0;
    ll msb = NetMan::GetPoNum() - 1;
    ll shift = msb - lsb;
    assert(iPatt < nFrame);
    assert(msb < 200);
    BigInt ret(0);
    for (ll k = msb; k >= lsb; --k) {
        ret <<= 1;
        if (dat[NetMan::GetPoId(k)][iPatt])
            ++ret;
    }
    if (isSign && ret >= (BigInt(1) << shift))
        ret = -((BigInt(1) << (shift + 1)) - ret);
    return ret;
}


BigInt Simulator::GetTempOutpPro(ll iPatt, bool isSign) const {
    ll lsb = 0;
    ll msb = NetMan::GetPoNum() - 1;
    ll shift = msb - lsb;
    assert(iPatt < nFrame);
    assert(msb < 200);
    BigInt ret(0);
    for (ll k = msb; k >= lsb; --k) {
        ret <<= 1;
        if (tempDat[NetMan::GetPoId(k)][iPatt])
            ++ret;
    }
    if (isSign && ret >= (BigInt(1) << shift))
        ret = -((BigInt(1) << (shift + 1)) - ret);
    return ret;
}


ll Simulator::GetOutpFast(ll iPatt, bool isSign) const {
    ll lsb = 0;
    ll msb = NetMan::GetPoNum() - 1;
    assert(iPatt < nFrame);
    assert(msb < 60);
    ll ret = 0;
    for (ll k = msb; k >= lsb; --k) {
        ret <<= 1;
        if (dat[NetMan::GetPoId(k)][iPatt])
            ++ret;
    }
    if (isSign && ret >= (1ll << (msb - lsb)))
        ret = -((1ll << (msb - lsb + 1)) - ret);
    return ret;
}


ll Simulator::GetTempOutpFast(ll iPatt, bool isSign) const {
    ll lsb = 0;
    ll msb = NetMan::GetPoNum() - 1;
    #ifdef DEBUG
    assert(iPatt < nFrame);
    assert(msb < 60);
    #endif
    ll ret = 0;
    for (ll k = msb; k >= lsb; --k) {
        ret <<= 1;
        if (tempDat[NetMan::GetPoId(k)][iPatt])
            ++ret;
    }
    if (isSign && ret >= (1ll << (msb - lsb)))
        ret = -((1ll << (msb - lsb + 1)) - ret);
    return ret;
}


void Simulator::PrintOutpStream(ll iPatt) const {
    #ifdef DEBUG
    assert(iPatt < nFrame);
    #endif
    for (ll k = GetPoNum() - 1; k >= 0; --k)
        cout << dat[GetPoId(k)][iPatt];
    cout << endl;
}


double Simulator::GetSignalProb(ll objId) const {
    #ifdef DEBUG
    assert(objId < NetMan::GetIdMaxPlus1());
    #endif
    return dat[objId].count() / static_cast <double> (nFrame);
}


void Simulator::PrintSignalProb() const {
    for (ll i = 0; i < NetMan::GetPoNum(); ++i) {
        cout << NetMan::GetName(NetMan::GetPo(i)) << " " << GetSignalProb(NetMan::GetPoId(i)) << endl;
    }
}


bool Simulator::IsPIOSame(const Simulator & oth_smlt) const {
    if (this->GetPiNum() != oth_smlt.GetPiNum())
        return false;
    for (ll i = 0; i < this->GetPiNum(); ++i) {
        // if (strcmp(Abc_ObjName(Abc_NtkPo(pNtk1, i)), Abc_ObjName(Abc_NtkPo(pNtk2, i))) != 0)
        if (this->GetPiName(i) != oth_smlt.GetPiName(i))
            return false;
    }
    if (this->GetPoNum() != oth_smlt.GetPoNum())
        return false;
    for (ll i = 0; i < this->GetPoNum(); ++i) {
        // if (strcmp(Abc_ObjName(Abc_NtkPo(pNtk1, i)), Abc_ObjName(Abc_NtkPo(pNtk2, i))) != 0)
        if (this->GetPoName(i) != oth_smlt.GetPoName(i))
            return false;
    }
    return true;
}


double Simulator::GetErrRate(const Simulator & oth_smlt, bool isCheck) const {
    if (isCheck)
        assert(IsPIOSame(oth_smlt));
    BitVect temp(nFrame, 0);
    for (ll i = 0; i < NetMan::GetPoNum(); ++i)
        temp |= (this->dat[NetMan::GetPoId(i)] ^ oth_smlt.dat[oth_smlt.GetPoId(i)]);
    return temp.count() / static_cast <double> (nFrame);
}


double Simulator::GetMeanErrDist(const Simulator & oth_smlt, bool isSign, bool isCheck) const {
    if (isCheck) {
        assert(IsPIOSame(oth_smlt));
        // assert(GetPoNum() <= 60);
    }
    // const ll warnPoint = numeric_limits <ll>::max() * 0.8;
    BigInt sed(0);
    for (int i = 0; i < nFrame; ++i) {
        auto accOut = GetOutpPro(i, isSign);
        auto appOut = oth_smlt.GetOutpPro(i, isSign);
        sed += abs(accOut - appOut);
        // assert(sed < warnPoint);
    }
    auto med = BigFlt(sed) / BigFlt(nFrame);
    assert(med < DBL_MAX);
    return (double)(med);
}


double Simulator::GetMeanSquareErr(const Simulator & oth_smlt, bool isSign, bool isCheck) const {
    if (isCheck) {
        assert(IsPIOSame(oth_smlt));
        assert(GetPoNum() < 200);
    }
    BigInt sse = 0;
    for (ll i = 0; i < nFrame; ++i) {
        BigInt accOut = GetOutpPro(i, isSign);
        BigInt appOut = oth_smlt.GetOutpPro(i, isSign);
        sse += (accOut - appOut) * (accOut - appOut);
    }
    return double(BigFlt(sse) / BigFlt(nFrame));
}


double Simulator::GetMeanHammDist(const Simulator & oth_smlt, bool isCheck) const {
    if (isCheck)
        assert(IsPIOSame(oth_smlt));
    ll sumHammDist = 0;
    for (int i = 0; i < NetMan::GetPoNum(); ++i) {
        auto diff = dat[NetMan::GetPoId(i)] ^ oth_smlt.dat[oth_smlt.GetPoId(i)];
        sumHammDist += diff.count();
    }
    return sumHammDist / static_cast<double>(nFrame);
}


double Simulator::GetSigNoiseRat(const Simulator & oth_smlt, bool isSign, bool isCheck) const {
    if (isCheck) {
        assert(IsPIOSame(oth_smlt));
        assert(GetPoNum() < 200);
    }
    BigInt sse = 0;
    BigInt sumAcc2 = 0;
    for (ll i = 0; i < nFrame; ++i) {
        BigInt accOut = GetOutpPro(i, isSign);
        BigInt appOut = oth_smlt.GetOutpPro(i, isSign);
        sse += (accOut - appOut) * (accOut - appOut);
        sumAcc2 += accOut * accOut;
    }
    if (sse != 0) {
        auto rat = BigFlt(sumAcc2) / BigFlt(sse);
        return double(BigFlt(10) * log10(rat));
    }
    return numeric_limits <double>::max();
}


double Simulator::GetError() const {
    BigInt res = 0;
    for (int i = GetPoNum() - 1; i >= 0; --i) {
        res <<= 1;
        res += (ll)(dat[GetPoId(i)].count());
    }
    auto ret = static_cast<double>(static_cast<BigFlt>(res) / static_cast<BigFlt>(nFrame));
    return ret;
}


void Simulator::CalcLocBoolDiff(Abc_Obj_t * pObj, list <Abc_Obj_t *> & disjCut, vector <Abc_Obj_t *> & cutNtk, vector < BitVect > & bdCut2Node) {
    #ifdef DEBUG
    assert(pObj->pNtk == GetNet());
    #endif
    if (tempDat.size() != dat.size())
        tempDat.resize(dat.size(), BitVect(nFrame, 0));
    // flip the node
    tempDat[pObj->Id] = ~dat[pObj->Id];
    // simulate
    Abc_NtkIncrementTravId(GetNet());
    Abc_NodeSetTravIdCurrent(pObj);
    for (auto & pInner: cutNtk)
        Abc_NodeSetTravIdCurrent(pInner);
    auto type = NetMan::GetNetType();
    if (type == NET_TYPE::AIG)
        assert(0);
        // UpdAigNodeForBoolAndPartDiff(pObj);
    else if (type == NET_TYPE::SOP) {
        for (auto & pInner: cutNtk)
            UpdSopNodeForBoolAndPartDiff(pInner);
    }
    else if (type == NET_TYPE::GATE) {
        for (auto & pInner: cutNtk)
            UpdGateNodeForBoolAndPartDiff(pInner);
    }
    else
        assert(0);
    // get boolean difference from the node to its disjoint cuts
    bdCut2Node.resize(disjCut.size(), BitVect(nFrame, 0));
    ll i = 0;
    for (auto & pCut: disjCut) {
        bdCut2Node[i] = dat[pCut->Id] ^ tempDat[pCut->Id];
        ++i;
    }
}


void Simulator::CalcLocPartDiff(Abc_Obj_t * pObj, list <Abc_Obj_t *> & disjCut, vector <Abc_Obj_t *> & cutNtk, vector < vector <int8_t> > & pdCut2Node) {
    #ifdef DEBUG
    assert(pObj->pNtk == GetNet());
    #endif
    if (tempDat.size() != dat.size())
        tempDat.resize(dat.size(), BitVect(nFrame, 0));
    // flip the node
    tempDat[pObj->Id] = ~dat[pObj->Id];
    // simulate
    Abc_NtkIncrementTravId(GetNet());
    Abc_NodeSetTravIdCurrent(pObj);
    for (auto & pInner: cutNtk)
        Abc_NodeSetTravIdCurrent(pInner);
    auto type = NetMan::GetNetType();
    if (type == NET_TYPE::AIG)
        assert(0);
    else if (type == NET_TYPE::SOP) {
        for (auto & pInner: cutNtk)
            UpdSopNodeForBoolAndPartDiff(pInner);
    }
    else if (type == NET_TYPE::GATE) {
        for (auto & pInner: cutNtk)
            UpdGateNodeForBoolAndPartDiff(pInner);
    }
    else
        assert(0);
    // get boolean difference from the node to its disjoint cuts
    pdCut2Node.resize(disjCut.size(), vector <int8_t> (nFrame, 0));
    ll i = 0;
    // parallel acceleration
    for (auto & pCut: disjCut) {
        // if (useMP) {
        //     omp_set_num_threads(numOfThread);
        //     #pragma omp parallel for schedule(dynamic)
        //     for (ll j = 0; j < nFrame; ++j) {
        //         pdCut2Node[i][j] = ((int8_t)1 - ((int8_t)dat[pObj->Id][j] << int8_t(1))) * ((int8_t)tempDat[pCut->Id][j] - (int8_t)dat[pCut->Id][j]);
        //     }
        // }
        // else 
        {
            for (ll j = 0; j < nFrame; ++j) {
                pdCut2Node[i][j] = ((int8_t)1 - ((int8_t)dat[pObj->Id][j] << int8_t(1))) * ((int8_t)tempDat[pCut->Id][j] - (int8_t)dat[pCut->Id][j]);
            }
        }
        ++i;
    }
}


void Simulator::UpdSopNodeForBoolAndPartDiff(Abc_Obj_t * pObj) {
    #ifdef DEBUG
    assert(!Abc_ObjIsPi(pObj));
    assert(!Abc_NodeIsConst(pObj));
    #endif
    if (Abc_ObjIsPo(pObj)) {
        #ifdef DEBUG
        assert(!Abc_ObjIsComplement(pObj));
        #endif
        Abc_Obj_t * pDriver = Abc_ObjFanin0(pObj);
        if (Abc_NodeIsTravIdCurrent(pDriver))
            tempDat[pObj->Id] = tempDat[pDriver->Id];
        else
            tempDat[pObj->Id] = dat[pDriver->Id];
        return;
    }
    // update sop
    char * pSop = static_cast <char *> (pObj->pData);
    UpdSopForBoolAndPartDiff(pObj, pSop);
}


void Simulator::UpdGateNodeForBoolAndPartDiff(Abc_Obj_t * pObj) {
    #ifdef DEBUG
    assert(!Abc_ObjIsPi(pObj));
    assert(!Abc_NodeIsConst(pObj));
    #endif
    if (Abc_ObjIsPo(pObj)) {
        Abc_Obj_t * pDriver = Abc_ObjFanin0(pObj);
        if (Abc_NodeIsTravIdCurrent(pDriver))
            tempDat[pObj->Id] = tempDat[pDriver->Id];
        else
            tempDat[pObj->Id] = dat[pDriver->Id];
        return;
    }
    // update sop
    char * pSop = static_cast <char *> ((static_cast <Mio_Gate_t *> (pObj->pData))->pSop);
    UpdSopForBoolAndPartDiff(pObj, pSop);
}


void Simulator::UpdSopForBoolAndPartDiff(Abc_Obj_t * pObj, char * pSop) {
    ll nVars = Abc_SopGetVarNum(pSop);
    BitVect product(nFrame, 0);
    for (char * pCube = pSop; *pCube; pCube += nVars + 3) {
        bool isFirst = true;
        for (ll i = 0; pCube[i] != ' '; i++) {
            Abc_Obj_t * pFanin = Abc_ObjFanin(pObj, i);
            BitVect & datFi = Abc_NodeIsTravIdCurrent(pFanin)? tempDat[pFanin->Id]: dat[pFanin->Id];
            switch (pCube[i]) {
                case '-':
                    continue;
                    break;
                case '0':
                    if (isFirst) {
                        isFirst = false;
                        product = ~datFi;
                    }
                    else
                        product &= ~datFi;
                    break;
                case '1':
                    if (isFirst) {
                        isFirst = false;
                        product = datFi;
                    }
                    else
                        product &= datFi;
                    break;
                default:
                    assert(0);
            }
        }
        if (isFirst) {
            isFirst = false;
            product.set();
        }
        #ifdef DEBUG
        assert(!isFirst);
        #endif
        if (pCube == pSop) {
            tempDat[pObj->Id] = product;
        }
        else
            tempDat[pObj->Id] |= product;
    }

    // complement
    if (Abc_SopIsComplement(pSop))
        tempDat[pObj->Id].flip();
}


bool IsPIOSame(Simulator & smlt0, Simulator & smlt1) {
    if (smlt0.GetPiNum() != smlt1.GetPiNum())
        return false;
    for (ll i = 0; i < smlt0.GetPiNum(); ++i) {
        // if (strcmp(Abc_ObjName(Abc_NtkPo(pNtk1, i)), Abc_ObjName(Abc_NtkPo(pNtk2, i))) != 0)
        if (smlt0.GetPiName(i) != smlt1.GetPiName(i))
            return false;
    }
    if (smlt0.GetPoNum() != smlt1.GetPoNum())
        return false;
    for (ll i = 0; i < smlt0.GetPoNum(); ++i) {
        // if (strcmp(Abc_ObjName(Abc_NtkPo(pNtk1, i)), Abc_ObjName(Abc_NtkPo(pNtk2, i))) != 0)
        if (smlt0.GetPoName(i) != smlt1.GetPoName(i))
            return false;
    }
    return true;
}


bool IsPIOSame(NetMan & net0, NetMan & net1) {
    if (net0.GetPiNum() != net1.GetPiNum())
        return false;
    for (ll i = 0; i < net0.GetPiNum(); ++i) {
        // if (strcmp(Abc_ObjName(Abc_NtkPo(pNtk1, i)), Abc_ObjName(Abc_NtkPo(pNtk2, i))) != 0)
        if (net0.GetPiName(i) != net1.GetPiName(i))
            return false;
    }
    if (net0.GetPoNum() != net1.GetPoNum())
        return false;
    for (ll i = 0; i < net0.GetPoNum(); ++i) {
        // if (strcmp(Abc_ObjName(Abc_NtkPo(pNtk1, i)), Abc_ObjName(Abc_NtkPo(pNtk2, i))) != 0)
        if (net0.GetPoName(i) != net1.GetPoName(i))
            return false;
    }
    return true;
}


void GetNewValue(Simulator & smlt, const IntVect& faninIds, const std::string & sop, RETURN_VAR BitVect & value) {
    if (sop == " 0\n") {
        value.reset();
        return;
    }
    if (sop == " 1\n") {
        value.set();
        return;
    }
    char * pSop = const_cast <char *> (sop.c_str());
    ll nVars = Abc_SopGetVarNum(pSop);
    assert(nVars == faninIds.size());

    ll nFrame = smlt.GetFrameNumb();
    BitVect product(nFrame, 0);
    for (char * pCube = pSop; *pCube; pCube += nVars + 3) {
        bool isFirst = true;
        for (ll i = 0; pCube[i] != ' '; i++) {
            ll faninId = faninIds[i];
            switch (pCube[i]) {
                case '-':
                    continue;
                    break;
                case '0':
                    if (isFirst) {
                        isFirst = false;
                        product = ~(*smlt.GetDat(faninId));
                    }
                    else
                        product &= ~(*smlt.GetDat(faninId));
                    break;
                case '1':
                    if (isFirst) {
                        isFirst = false;
                        product = (*smlt.GetDat(faninId));
                    }
                    else
                        product &= (*smlt.GetDat(faninId));
                    break;
                default:
                    assert(0);
            }
        }
        if (isFirst) {
            isFirst = false;
            product.set();
        }
        assert(!isFirst);
        if (pCube == pSop)
            value = product;
        else
            value |= product;
    }

    // complement
    if (Abc_SopIsComplement(pSop))
        value.flip();
}