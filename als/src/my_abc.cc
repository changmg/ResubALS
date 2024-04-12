#include "my_abc.h"


using namespace abc;
using namespace std;


std::ostream & operator << (std::ostream & os, const NET_TYPE netwType) {
    const std::string strs[4] = {"AIG", "GATE", "SOP", "STRASH"};
    os << strs[static_cast <ll> (netwType)];
    return os;
}


std::ostream & operator << (std::ostream & os, const ORIENT orient) {
    const std::string strs[2] = {"AREA", "DELAY"};
    os << strs[static_cast <ll> (orient)];
    return os;
}


std::ostream & operator << (std::ostream & os, const MAP_TYPE cell) {
    const std::string strs[2] = {"LUT", "SCL"};
    os << strs[static_cast <ll> (cell)];
    return os;
}


AbcMan::AbcMan() {
    #ifdef DEBUG
    assert(Abc_FrameGetGlobalFrame() != nullptr);
    #endif
}


void AbcMan::Comm(const string & cmd, bool isVerb) {
    if (isVerb)
        cout << "Execute abc command: " << cmd << endl;
    if (Cmd_CommandExecute(GetAbcFame(), cmd.c_str())) {
        cout << "Execuation failed." << endl;
        assert(0);
    }
}


void AbcMan::ReadNet(const std::string & fileName, bool inpMapVerilog) {
    #ifdef DEBUG
    assert(IsPathExist(fileName));
    #endif
    if (inpMapVerilog)
        Comm("r -m " + fileName);
    else
        Comm("r " + fileName);
}


void AbcMan::WriteNet(const std::string & fileName, bool isVerb) {
    Comm("w " + fileName, isVerb);
}


void AbcMan::ReadStandCell(const std::string & fileName) {
    #ifdef DEBUG
    assert(IsPathExist(fileName));
    #endif
    Comm("r " + fileName);
}


void AbcMan::ConvToAig() {
    Comm("aig");
}


void AbcMan::ConvToGate() {
    Map(MAP_TYPE::SCL, ORIENT::AREA);
}


void AbcMan::ConvToSop() {
    if (GetNetType() == NET_TYPE::STRASH)
        Comm("logic;");
    Comm("sop");
}


void AbcMan::ConvToStrash() {
    Comm("st");
}


void AbcMan::PrintStat() {
    if (GetNetType() == NET_TYPE::GATE && GetAbcFame()->pLibScl != nullptr) {
        TopoSort();
        StatTimeAnal();
        // Comm("ps");
    }
    else {
        Comm("ps");
    }
}


void AbcMan::TopoSort() {
    assert(GetNetType() == NET_TYPE::AIG || GetNetType() == NET_TYPE::SOP || GetNetType() == NET_TYPE::GATE);
    Comm("topo");

    // fix twin nodes
    auto pNtk = GetNet();
    if (Abc_NtkHasMapping(pNtk)) {
        Abc_Ntk_t * pNtkNew; 
        Abc_Obj_t * pObj, * pFanin;
        int i, k;
        assert(pNtk != nullptr);
        // start the network
        pNtkNew = Abc_NtkStartFrom( pNtk, pNtk->ntkType, pNtk->ntkFunc );
        // copy the internal nodes
        assert(!Abc_NtkIsStrash(pNtk));
        // duplicate the nets and nodes (CIs/COs/latches already dupped)
        set <ll> skip;
        Abc_NtkForEachObj( pNtk, pObj, i ) {
            if ( pObj->pCopy == NULL && skip.count(pObj->Id) == 0 ) {
                Abc_NtkDupObj(pNtkNew, pObj, Abc_NtkHasBlackbox(pNtk) && Abc_ObjIsNet(pObj));
                auto pTwin = GetTwinNode(pObj);
                if (pTwin != nullptr) {
                    Abc_NtkDupObj(pNtkNew, pTwin, Abc_NtkHasBlackbox(pNtk) && Abc_ObjIsNet(pTwin));
                    skip.insert(pTwin->Id);
                }
            }
        }
        // reconnect all objects (no need to transfer attributes on edges)
        Abc_NtkForEachObj( pNtk, pObj, i )
            if ( !Abc_ObjIsBox(pObj) && !Abc_ObjIsBo(pObj) )
                Abc_ObjForEachFanin( pObj, pFanin, k )
                    Abc_ObjAddFanin( pObj->pCopy, pFanin->pCopy );
        // duplicate the EXDC Ntk
        if ( pNtk->pExdc )
            pNtkNew->pExdc = Abc_NtkDup( pNtk->pExdc );
        if ( pNtk->pExcare )
            pNtkNew->pExcare = Abc_NtkDup( (Abc_Ntk_t *)pNtk->pExcare );
        // duplicate timing manager
        if ( pNtk->pManTime )
            Abc_NtkTimeInitialize( pNtkNew, pNtk );
        if ( pNtk->vPhases )
            Abc_NtkTransferPhases( pNtkNew, pNtk );
        if ( pNtk->pWLoadUsed )
            pNtkNew->pWLoadUsed = Abc_UtilStrsav( pNtk->pWLoadUsed );
        // check correctness
        if ( !Abc_NtkCheck( pNtkNew ) )
            fprintf( stdout, "Abc_NtkDup(): Network check has failed.\n" );
        pNtk->pCopy = pNtkNew;
        // return pNtkNew;
        SetMainNetw(pNtkNew);
    }
}


void AbcMan::StatTimeAnal() {
    #ifdef DEBUG
    assert(GetNetType() == NET_TYPE::GATE);
    assert(GetAbcFame()->pLibScl != nullptr);
    #endif
    TopoSort();
    Comm("stime");
}


void AbcMan::Synth(ORIENT orient, bool isVerb) {
    #ifdef DEBUG
    assert(Abc_NtkIsStrash(GetNet()));
    #endif
    if (isVerb)
        cout << orient << "-oriented synthesis" << endl;
    const ll commSize = 2;
    const string areaComm[commSize] = {"st; drwsat", "st; compress2rs"};
    const string delayComm[commSize] = {"st; dc2", "st; resyn2"};
    double oldArea = GetArea();
    double oldDelay = GetDelay();
    bool isCont = true;
    while (isCont) {
        isCont = false;
        for (ll i = 0; i < commSize; ++i) {
            auto oldNtk = Abc_NtkDup(GetNet());
            if (orient == ORIENT::AREA)
                Comm(areaComm[i], isVerb);
            else if (orient == ORIENT::DELAY)
                Comm(delayComm[i], isVerb);
            else
                assert(0);
            auto res = make_pair <double, double> (GetArea(), GetDelay());
            if (isVerb)
                PrintStat();
            double newArea = res.first;
            double newDelay = res.second;
            IMPR impr = UpdNetw(oldArea, oldDelay, oldNtk, newArea, newDelay, orient);
            if (impr == IMPR::GOOD) {
                oldArea = newArea;
                oldDelay = newDelay;
                isCont = true;
            }
            if (isVerb)
                cout << (impr == IMPR::GOOD? "accept": "cancel") << endl;
        }
    }
    if (isVerb)
        PrintStat();
}


void AbcMan::SynthWithResyn2Comm() {
    Comm("st; resyn2; logic; sop; ps;");
}


void AbcMan::SynthAndMap(double maxDelay, bool isVerb) {
    bool cont = true;
    if (isVerb)
        cout << "maxDelay = " << maxDelay << endl;
    TopoSort();
    while (cont) {
        double oldArea = numeric_limits <double>::max(), oldDelay = numeric_limits <double>::max();
        if (GetNetType(GetNet()) == NET_TYPE::GATE)
            oldArea = GetArea(), oldDelay = GetDelay();
        auto pOldNtk = Abc_NtkDup(GetNet());
        if (isVerb)
            cout << "oldArea = " << oldArea << ", " << "oldDelay = " << oldDelay << endl;
        Comm("st; resyn2; dch; amap;", isVerb);
        TopoSort();
        double newArea = GetArea(), newDelay = GetDelay();
        if (isVerb)
            cout << "newArea = " << newArea << ", " << "newDelay = " << newDelay << endl;
        if (newDelay <= maxDelay) {
            auto impr = UpdNetw(oldArea, oldDelay, pOldNtk, newArea, newDelay, ORIENT::AREA);
            if (impr != IMPR::GOOD) {
                cont = false;
                if (isVerb)
                    cout << "reject" << endl;
            }
            else {
                if (isVerb)
                    cout << "accept" << endl;
            }
        }
        else {
            SetMainNetw(pOldNtk);
                cont = false;
            if (isVerb)
                cout << "reject" << endl;
        }
    }
    PrintStat();
}


void AbcMan::Sweep() {
    #ifdef DEBUG
    assert(GetNetType() == NET_TYPE::SOP);
    #endif
    Comm("sweep; sop;");
}


pair <double, double> AbcMan::Map(MAP_TYPE cell, ORIENT orient, bool isVerb) {
    double oldArea = numeric_limits <double>::max();
    double oldDelay = numeric_limits <double>::max();
    ostringstream LutInpStr("");
    LutInpStr << LutInp;
    if ((cell == MAP_TYPE::SCL && GetNetType() == NET_TYPE::GATE) ||
        (cell == MAP_TYPE::LUT && IsLutNetw())) {
        oldArea = GetArea();
        oldDelay = GetDelay();
    }
    bool isFirst = true;
    bool isCont = true;
    while (isCont) {
        auto oldNtk = Abc_NtkDup(GetNet());
        if (isFirst) {
            Comm("st; dch;", isVerb);
            isFirst = false;
        }
        else
            Comm("st; b;", isVerb);
        if (cell == MAP_TYPE::SCL) {
            if (orient == ORIENT::AREA)
                Comm("amap", isVerb);
            else if (orient == ORIENT::DELAY)
                Comm("map", isVerb);
            else
                assert(0);
        }
        else if (cell == MAP_TYPE::LUT) {
            if (orient == ORIENT::AREA)
                Comm("if -a -K " + LutInpStr.str(), isVerb);
            else if (orient == ORIENT::DELAY)
                Comm("if -K " + LutInpStr.str(), isVerb);
            else
                assert(0);
        }
        else
            assert(0);
        double newArea = GetArea();
        double newDelay = GetDelay();
        IMPR impr = UpdNetw(oldArea, oldDelay, oldNtk, newArea, newDelay, orient);
        if (impr == IMPR::GOOD) {
            oldArea = newArea;
            oldDelay = newDelay;
        }
        else
            isCont = false;
        PrintStat();
    }
    return make_pair(oldArea, oldDelay);
}


pair <double, double> AbcMan::Map2(double maxDelay, bool isVerb) {
    double oldArea = numeric_limits <double>::max();
    double oldDelay = numeric_limits <double>::max();
    ostringstream LutInpStr("");
    LutInpStr << LutInp;
    assert(GetNetType() == NET_TYPE::STRASH);
    bool isFirst = true;
    bool isCont = true;
    while (isCont) {
        auto oldNtk = Abc_NtkDup(GetNet());
        if (isFirst) {
            Comm("st; dch;", isVerb);
            isFirst = false;
        }
        else
            Comm("st; b;", isVerb);
        ostringstream oss("");
        oss << "map -D " << maxDelay;
        Comm(oss.str(), isVerb);
        double newArea = GetArea();
        double newDelay = GetDelay();
        IMPR impr = UpdNetw(oldArea, oldDelay, oldNtk, newArea, newDelay, ORIENT::AREA);
        if (impr == IMPR::GOOD) {
            oldArea = newArea;
            oldDelay = newDelay;
        }
        else
            isCont = false;
        // PrintStat();
    }
    return make_pair(oldArea, oldDelay);
}


IMPR AbcMan::UpdNetw(double oldArea, double oldDelay, Abc_Ntk_t * oldNtk, double newArea, double newDelay, ORIENT orient) {
    IMPR impr = IMPR::SAME;
    if (orient == ORIENT::AREA) {
        if (DoubleGreat(newArea, oldArea) || (DoubleEqual(newArea, oldArea) && DoubleGreat(newDelay, oldDelay)))
            impr = IMPR::BAD;
        else if (DoubleEqual(newArea, oldArea) && DoubleEqual(newDelay, oldDelay))
            impr = IMPR::SAME;
        else
            impr = IMPR::GOOD;
    }
    else if (orient == ORIENT::DELAY) {
        if (DoubleGreat(newDelay, oldDelay) || (DoubleEqual(newDelay, oldDelay) && DoubleGreat(newArea, oldArea)))
            impr = IMPR::BAD;
        else if (DoubleEqual(newDelay, oldDelay) && DoubleEqual(newArea, oldArea))
            impr = IMPR::SAME;
        else
            impr = IMPR::GOOD;
    }
    else
        assert(0);
    if (impr == IMPR::BAD) {
        #ifdef DEBUG
        assert(oldArea != numeric_limits <double>::max() && oldDelay != numeric_limits <double>::max());
        assert(oldNtk != nullptr);
        // cout << "Cancel the last abc command" << endl;
        #endif
        SetMainNetw(oldNtk);
    }
    else
        Abc_NtkDelete(oldNtk);
    return impr;
}


NET_TYPE AbcMan::GetNetType(Abc_Ntk_t * pNtk) const {
    if (Abc_NtkIsAigLogic(pNtk))
        return NET_TYPE::AIG;
    else if (Abc_NtkIsMappedLogic(pNtk))
        return NET_TYPE::GATE;
    else if (Abc_NtkIsSopLogic(pNtk))
        return NET_TYPE::SOP;
    else if (Abc_NtkIsStrash(pNtk))
        return NET_TYPE::STRASH;
    else {
        cout << pNtk << endl;
        cout << "invalid network type" << endl;
        assert(0);
        return NET_TYPE::AIG;
    }
}


double AbcMan::GetArea(Abc_Ntk_t * pNtk) const {
    auto type = GetNetType(pNtk);
    if (type == NET_TYPE::AIG || type == NET_TYPE::STRASH)
        return Abc_NtkNodeNum(pNtk);
    else if (type == NET_TYPE::SOP) {
        Abc_Obj_t * pObj = nullptr;
        ll i = 0;
        ll ret = Abc_NtkNodeNum(pNtk);
        Abc_NtkForEachNode(pNtk, pObj, i) {
            if (Abc_NodeIsConst(pObj))
                --ret;
        }
        return ret;
    }
    else if (type == NET_TYPE::GATE) {
        auto pLibScl = static_cast <SC_Lib *> (GetAbcFame()->pLibScl);
        if (pLibScl == nullptr)
            return Abc_NtkGetMappedArea(pNtk);
        else {
            return Abc_NtkGetMappedArea(pNtk);
            // #ifdef DEBUG
            // assert(pNtk->nBarBufs2 == 0);
            // assert(CheckSCLNet(pNtk));
            // #endif
            // SC_Man * p = Abc_SclManStart(pLibScl, pNtk, 0, 1, 0, 0);
            // double area = Abc_SclGetTotalArea(p->pNtk);
            // Abc_SclManFree(p);
            // return area;
        }
    }
    else
        assert(0);
    return 0.0;
}


double AbcMan::GetDelay(Abc_Ntk_t * pNtk) const {
    auto type = GetNetType(pNtk);
    if (type == NET_TYPE::AIG || type == NET_TYPE::SOP || type == NET_TYPE::STRASH)
        return Abc_NtkLevel(pNtk);
    else if (type == NET_TYPE::GATE) {
        auto pLibScl = static_cast <SC_Lib *> (GetAbcFame()->pLibScl);
        if (pLibScl == nullptr) {
            return Abc_NtkDelayTrace(pNtk, nullptr, nullptr, 0);
        }
        else {
            assert(pNtk->nBarBufs2 == 0);
            assert(CheckSCLNet(pNtk));
            SC_Man * p = Abc_SclManStart(pLibScl, pNtk, 0, 1, 0, 0);
            int fRise = 0;
            Abc_Obj_t * pPivot = Abc_SclFindCriticalCo(p, &fRise); 
            double delay = Abc_SclObjTimeOne(p, pPivot, fRise);
            Abc_Obj_t * pObj = nullptr;
            ll i = 0;
            Abc_NtkForEachObj(pNtk, pObj, i)
                pObj->dTemp = Abc_SclObjTimeMax(p, pObj);
            Abc_SclManFree(p);
            return delay;
        }
    }
    else
        assert(0);
    return 0.0;
}


bool AbcMan::CheckSCLNet(abc::Abc_Ntk_t * pNtk) const {
    Abc_Obj_t * pObj, * pFanin;
    int i, k, fFlag = 1;
    Abc_NtkIncrementTravId( pNtk );        
    Abc_NtkForEachCi( pNtk, pObj, i )
        Abc_NodeSetTravIdCurrent( pObj );
    Abc_NtkForEachNode( pNtk, pObj, i )
    {
        Abc_ObjForEachFanin( pObj, pFanin, k )
            if ( !Abc_NodeIsTravIdCurrent( pFanin ) )
                printf( "obj %d and its fanin %d are not in the topo order\n", Abc_ObjId(pObj), Abc_ObjId(pFanin) ), fFlag = 0;
        Abc_NodeSetTravIdCurrent( pObj );
        if ( Abc_ObjIsBarBuf(pObj) )
            continue;
        // if ( Abc_ObjFanoutNum(pObj) == 0 )
        //     printf( "node %d has no fanout\n", Abc_ObjId(pObj) ), fFlag = 0;
        if ( !fFlag )
            break;
    }
    // if ( fFlag && fVerbose )
    //     printf( "The network is in topo order and no dangling nodes.\n" );
    return fFlag;
}


Abc_Obj_t * AbcMan::GetTwinNode( Abc_Obj_t * pNode ) {
    assert( Abc_NtkHasMapping(pNode->pNtk) );
    Mio_Gate_t * pGate = (Mio_Gate_t *)pNode->pData;
    if ( pGate == nullptr || Mio_GateReadTwin(pGate) == nullptr )
        return nullptr;
    Abc_Obj_t * pNode2 = nullptr;
    ll id = 0;
    Abc_Obj_t * pTwin = nullptr;
    ll count = 0;
    Abc_NtkForEachNode(pNode->pNtk, pNode2, id) {
        if ( Abc_ObjFaninNum(pNode) != Abc_ObjFaninNum(pNode2) )
            continue;
        bool sameFanin = true;
        for (ll faninId = 0; faninId < Abc_ObjFaninNum(pNode); ++faninId) {
            if (Abc_ObjFanin(pNode, faninId) != Abc_ObjFanin(pNode2, faninId)) {
                sameFanin = false;
                break;
            }
        }
        if (!sameFanin)
            continue;
        if ( Mio_GateReadTwin(pGate) != (Mio_Gate_t *)pNode2->pData )
            continue;
        pTwin = pNode2;
        ++count;
        if (count > 1)
            assert(0);
    }
    return pTwin;
}


void AbcMan::LoadAlias() {
    Comm("alias hi history", false);
    Comm("alias b balance", false);
    Comm("alias cg clockgate", false);
    Comm("alias cl cleanup", false);
    Comm("alias clp collapse", false);
    Comm("alias cs care_set", false);
    Comm("alias el eliminate", false);
    Comm("alias esd ext_seq_dcs", false);
    Comm("alias f fraig", false);
    Comm("alias fs fraig_sweep", false);
    Comm("alias fsto fraig_store", false);
    Comm("alias fres fraig_restore", false);
    Comm("alias fr fretime", false);
    Comm("alias ft fraig_trust", false);
    Comm("alias ic indcut", false);
    Comm("alias lp lutpack", false);
    Comm("alias pcon print_cone", false);
    Comm("alias pd print_dsd", false);
    Comm("alias pex print_exdc -d", false);
    Comm("alias pf print_factor", false);
    Comm("alias pfan print_fanio", false);
    Comm("alias pg print_gates", false);
    Comm("alias pl print_level", false);
    Comm("alias plat print_latch", false);
    Comm("alias pio print_io", false);
    Comm("alias pk print_kmap", false);
    Comm("alias pm print_miter", false);
    Comm("alias ps print_stats ", false);
    Comm("alias psb print_stats -b", false);
    Comm("alias psu print_supp", false);
    Comm("alias psy print_symm", false);
    Comm("alias pun print_unate", false);
    Comm("alias q quit", false);
    Comm("alias r read", false);
    Comm("alias ra read_aiger", false);
    Comm("alias r3 retime -M 3", false);
    Comm("alias r3f retime -M 3 -f", false);
    Comm("alias r3b retime -M 3 -b", false);
    Comm("alias ren renode", false);
    Comm("alias rh read_hie", false);
    Comm("alias ri read_init", false);
    Comm("alias rl read_blif", false);
    Comm("alias rb read_bench", false);
    Comm("alias ret retime", false);
    Comm("alias dret dretime", false);
    Comm("alias rp read_pla", false);
    Comm("alias rt read_truth", false);
    Comm("alias rv read_verilog", false);
    Comm("alias rvl read_verlib", false);
    Comm("alias rsup read_super mcnc5_old.super", false);
    Comm("alias rlib read_library", false);
    Comm("alias rlibc read_library cadence.genlib", false);
    Comm("alias rty read_liberty", false);
    Comm("alias rlut read_lut", false);
    Comm("alias rw rewrite", false);
    Comm("alias rwz rewrite -z", false);
    Comm("alias rf refactor", false);
    Comm("alias rfz refactor -z", false);
    Comm("alias re restructure", false);
    Comm("alias rez restructure -z", false);
    Comm("alias rs resub", false);
    Comm("alias rsz resub -z", false);
    Comm("alias sa set autoexec ps", false);
    Comm("alias scl scleanup", false);
    Comm("alias sif if -s", false);
    Comm("alias so source -x", false);
    Comm("alias st strash", false);
    Comm("alias sw sweep", false);
    Comm("alias ssw ssweep", false);
    Comm("alias tr0 trace_start", false);
    Comm("alias tr1 trace_check", false);
    Comm("alias trt \"r c.blif; st; tr0; b; tr1\"", false);
    Comm("alias u undo", false);
    Comm("alias w write", false);
    Comm("alias wa write_aiger", false);
    Comm("alias wb write_bench", false);
    Comm("alias wc write_cnf", false);
    Comm("alias wh write_hie", false);
    Comm("alias wl write_blif", false);
    Comm("alias wp write_pla", false);
    Comm("alias wv write_verilog", false);
    Comm("alias resyn       \"b; rw; rwz; b; rwz; b\"", false);
    Comm("alias resyn2      \"b; rw; rf; b; rw; rwz; b; rfz; rwz; b\"", false);
    Comm("alias resyn2a     \"b; rw; b; rw; rwz; b; rwz; b\"", false);
    Comm("alias resyn3      \"b; rs; rs -K 6; b; rsz; rsz -K 6; b; rsz -K 5; b\"", false);
    Comm("alias compress    \"b -l; rw -l; rwz -l; b -l; rwz -l; b -l\"", false);
    Comm("alias compress2   \"b -l; rw -l; rf -l; b -l; rw -l; rwz -l; b -l; rfz -l; rwz -l; b -l\"", false);
    Comm("alias choice      \"fraig_store; resyn; fraig_store; resyn2; fraig_store; fraig_restore\"", false);
    Comm("alias choice2     \"fraig_store; balance; fraig_store; resyn; fraig_store; resyn2; fraig_store; resyn2; fraig_store; fraig_restore\"", false);
    Comm("alias rwsat       \"st; rw -l; b -l; rw -l; rf -l\"", false);
    Comm("alias drwsat2     \"st; drw; b -l; drw; drf; ifraig -C 20; drw; b -l; drw; drf\"", false);
    Comm("alias share       \"st; multi -m; sop; fx; resyn2\"", false);
    Comm("alias addinit     \"read_init; undc; strash; zero\"", false);
    Comm("alias blif2aig    \"undc; strash; zero\"", false);
    Comm("alias v2p         \"&vta_gla; &ps; &gla_derive; &put; w 1.aig; pdr -v\"", false);
    Comm("alias g2p         \"&ps; &gla_derive; &put; w 2.aig; pdr -v\"", false);
    Comm("alias &sw_        \"&put; sweep; st; &get\"", false);
    Comm("alias &fx_        \"&put; sweep; sop; fx; st; &get\"", false);
    Comm("alias &dc3        \"&b; &jf -K 6; &b; &jf -K 4; &b\"", false);
    Comm("alias &dc4        \"&b; &jf -K 7; &fx; &b; &jf -K 5; &fx; &b\"", false);
    Comm("alias src_rw      \"st; rw -l; rwz -l; rwz -l\"", false);
    Comm("alias src_rs      \"st; rs -K 6 -N 2 -l; rs -K 9 -N 2 -l; rs -K 12 -N 2 -l\"", false);
    Comm("alias src_rws     \"st; rw -l; rs -K 6 -N 2 -l; rwz -l; rs -K 9 -N 2 -l; rwz -l; rs -K 12 -N 2 -l\"", false);
    Comm("alias resyn2rs    \"b; rs -K 6; rw; rs -K 6 -N 2; rf; rs -K 8; b; rs -K 8 -N 2; rw; rs -K 10; rwz; rs -K 10 -N 2; b; rs -K 12; rfz; rs -K 12 -N 2; rwz; b\"", false);
    Comm("alias compress2rs \"b -l; rs -K 6 -l; rw -l; rs -K 6 -N 2 -l; rf -l; rs -K 8 -l; b -l; rs -K 8 -N 2 -l; rw -l; rs -K 10 -l; rwz -l; rs -K 10 -N 2 -l; b -l; rs -K 12 -l; rfz -l; rs -K 12 -N 2 -l; rwz -l; b -l\"", false);
    Comm("alias fix_aig     \"logic; undc; strash; zero\"", false);
    Comm("alias fix_blif    \"undc; strash; zero\"", false);
    Comm("alias recadd3     \"st; rec_add3; b; rec_add3; dc2; rec_add3; if -K 8; bidec; st; rec_add3; dc2; rec_add3; if -g -K 6; st; rec_add3\"", false);
}


NetMan::NetMan(): AbcMan(), pNtk(nullptr), isDupl(true) {
}


Abc_Obj_t * Abc_NtkDupObj_KeepName( Abc_Ntk_t * pNtkNew, Abc_Obj_t * pObj, int fCopyName ) {
    Abc_Obj_t * pObjNew;
    // create the new object
    pObjNew = Abc_NtkCreateObj( pNtkNew, (Abc_ObjType_t)pObj->Type );
    // transfer names of the terminal objects
    if ( fCopyName )
    {
        if ( Abc_ObjIsCi(pObj) )
        {
            if ( !Abc_NtkIsNetlist(pNtkNew) )
                Abc_ObjAssignName( pObjNew, Abc_ObjName(Abc_ObjFanout0Ntk(pObj)), NULL );
        }
        else if ( Abc_ObjIsCo(pObj) )
        {
            if ( !Abc_NtkIsNetlist(pNtkNew) )
            {
                if ( Abc_ObjIsPo(pObj) )
                    Abc_ObjAssignName( pObjNew, Abc_ObjName(Abc_ObjFanin0Ntk(pObj)), NULL );
                else
                {
                    assert( Abc_ObjIsLatch(Abc_ObjFanout0(pObj)) );
                    Abc_ObjAssignName( pObjNew, Abc_ObjName(pObj), NULL );
                }
            }
        }
        else if ( Abc_ObjIsBox(pObj) || Abc_ObjIsNet(pObj) || Abc_ObjIsNode(pObj) )
            Abc_ObjAssignName( pObjNew, Abc_ObjName(pObj), NULL );
    }
    // copy functionality/names
    if ( Abc_ObjIsNode(pObj) ) // copy the function if functionality is compatible
    {
        if ( pNtkNew->ntkFunc == pObj->pNtk->ntkFunc ) 
        {
            if ( Abc_NtkIsStrash(pNtkNew) ) 
            {}
            else if ( Abc_NtkHasSop(pNtkNew) || Abc_NtkHasBlifMv(pNtkNew) )
                pObjNew->pData = Abc_SopRegister( (Mem_Flex_t *)pNtkNew->pManFunc, (char *)pObj->pData );
#ifdef ABC_USE_CUDD
            else if ( Abc_NtkHasBdd(pNtkNew) )
                pObjNew->pData = Cudd_bddTransfer((DdManager *)pObj->pNtk->pManFunc, (DdManager *)pNtkNew->pManFunc, (DdNode *)pObj->pData), Cudd_Ref((DdNode *)pObjNew->pData);
#endif
            else if ( Abc_NtkHasAig(pNtkNew) )
                pObjNew->pData = Hop_Transfer((Hop_Man_t *)pObj->pNtk->pManFunc, (Hop_Man_t *)pNtkNew->pManFunc, (Hop_Obj_t *)pObj->pData, Abc_ObjFaninNum(pObj));
            else if ( Abc_NtkHasMapping(pNtkNew) )
                pObjNew->pData = pObj->pData, pNtkNew->nBarBufs2 += !pObj->pData;
            else assert( 0 );
        }
    }
    else if ( Abc_ObjIsNet(pObj) ) // copy the name
    {
    }
    else if ( Abc_ObjIsLatch(pObj) ) // copy the reset value
        pObjNew->pData = pObj->pData;
    pObjNew->fPersist = pObj->fPersist;
    // transfer HAIG
//    pObjNew->pEquiv = pObj->pEquiv;
    // remember the new node in the old node
    pObj->pCopy = pObjNew;
    return pObjNew;
}


Abc_Ntk_t * Abc_NtkDup_KeepName( Abc_Ntk_t * pNtk ) {
    Abc_Ntk_t * pNtkNew; 
    Abc_Obj_t * pObj, * pFanin;
    int i, k;
    if ( pNtk == NULL )
        return NULL;
    // start the network
    pNtkNew = Abc_NtkStartFrom( pNtk, pNtk->ntkType, pNtk->ntkFunc );
    // copy the internal nodes
    if ( Abc_NtkIsStrash(pNtk) )
    {
        // copy the AND gates
        Abc_AigForEachAnd( pNtk, pObj, i )
            pObj->pCopy = Abc_AigAnd( (Abc_Aig_t *)pNtkNew->pManFunc, Abc_ObjChild0Copy(pObj), Abc_ObjChild1Copy(pObj) );
        // relink the choice nodes
        Abc_AigForEachAnd( pNtk, pObj, i )
            if ( pObj->pData )
                pObj->pCopy->pData = ((Abc_Obj_t *)pObj->pData)->pCopy;
        // relink the CO nodes
        Abc_NtkForEachCo( pNtk, pObj, i )
            Abc_ObjAddFanin( pObj->pCopy, Abc_ObjChild0Copy(pObj) );
        // get the number of nodes before and after
        if ( Abc_NtkNodeNum(pNtk) != Abc_NtkNodeNum(pNtkNew) )
            printf( "Warning: Structural hashing during duplication reduced %d nodes (this is a minor bug).\n",
                Abc_NtkNodeNum(pNtk) - Abc_NtkNodeNum(pNtkNew) );
    }
    else
    {
        // duplicate the nets and nodes (CIs/COs/latches already dupped)
        Abc_NtkForEachObj( pNtk, pObj, i )
            if ( pObj->pCopy == NULL )
                Abc_NtkDupObj_KeepName(pNtkNew, pObj, 1);
        // reconnect all objects (no need to transfer attributes on edges)
        Abc_NtkForEachObj( pNtk, pObj, i )
            if ( !Abc_ObjIsBox(pObj) && !Abc_ObjIsBo(pObj) )
                Abc_ObjForEachFanin( pObj, pFanin, k )
                    Abc_ObjAddFanin( pObj->pCopy, pFanin->pCopy );
    }
    // duplicate the EXDC Ntk
    if ( pNtk->pExdc )
        pNtkNew->pExdc = Abc_NtkDup( pNtk->pExdc );
    if ( pNtk->pExcare )
        pNtkNew->pExcare = Abc_NtkDup( (Abc_Ntk_t *)pNtk->pExcare );
    // duplicate timing manager
    if ( pNtk->pManTime )
        Abc_NtkTimeInitialize( pNtkNew, pNtk );
    if ( pNtk->vPhases )
        Abc_NtkTransferPhases( pNtkNew, pNtk );
    if ( pNtk->pWLoadUsed )
        pNtkNew->pWLoadUsed = Abc_UtilStrsav( pNtk->pWLoadUsed );
    // check correctness
    if ( !Abc_NtkCheck( pNtkNew ) )
        fprintf( stdout, "Abc_NtkDup(): Network check has failed.\n" );
    pNtk->pCopy = pNtkNew;
    return pNtkNew;
}


NetMan::NetMan(Abc_Ntk_t * p_ntk, bool is_dupl): AbcMan(), isDupl(is_dupl) {
    if (is_dupl)
        pNtk = Abc_NtkDup_KeepName(p_ntk);
    else
        pNtk = p_ntk;
}


NetMan::NetMan(std::string & fileName): AbcMan() {
    AbcMan::ReadNet(fileName);
    pNtk = AbcMan::GetNet();
}


NetMan::~NetMan() {
    if (isDupl && pNtk != AbcMan::GetNet()) {
        if (pNtk != nullptr) {
            Abc_NtkDelete(pNtk);
            pNtk = nullptr;
        }
    }
}


NetMan::NetMan(const NetMan & net_man): AbcMan(), isDupl(true) {
    pNtk = Abc_NtkDup_KeepName(net_man.pNtk);
}


NetMan::NetMan(NetMan && net_man): AbcMan(), pNtk(net_man.pNtk), isDupl(net_man.isDupl) {
    net_man.isDupl = false;
    net_man.pNtk = nullptr;
}


NetMan & NetMan::operator = (const NetMan & net_man) {
    if (this == &net_man)
        return *this;
    if (isDupl && pNtk != nullptr && pNtk != AbcMan::GetNet() && pNtk != net_man.GetNet())
        Abc_NtkDelete(pNtk);
    pNtk = Abc_NtkDup_KeepName(net_man.GetNet());
    isDupl = true;
    return *this;
}


NetMan & NetMan::operator = (NetMan && net_man) {
    if (this == &net_man)
        return *this;
    if (isDupl && pNtk != nullptr && pNtk != AbcMan::GetNet() && pNtk != net_man.GetNet())
        Abc_NtkDelete(pNtk);
    pNtk = net_man.pNtk;
    isDupl = net_man.isDupl;
    net_man.isDupl = false;
    net_man.pNtk = nullptr;
    return *this;
}


pair <ll, ll> NetMan::GetConstId(bool isVerb) {
    pair <ll, ll> ret(-1, -1);
    auto type = GetNetType();
    Abc_Obj_t * pObj = nullptr;
    ll i = 0;
    Abc_NtkForEachNode(GetNet(), pObj, i) {
        if (type == NET_TYPE::GATE || type == NET_TYPE::SOP) {
            if (Abc_NodeIsConst0(pObj)) {
                if (isVerb)
                    cout << "find const 0: " << pObj << endl;
                if (ret.first == -1)
                    ret.first = GetId(pObj);
            }
            else if (Abc_NodeIsConst1(pObj)) {
                if (isVerb)
                    cout << "find const 1: " << pObj << endl;
                if (ret.second == -1)
                    ret.second = GetId(pObj);
            }
        }
        else if (type == NET_TYPE::AIG) {
            auto pHopObj = static_cast <Hop_Obj_t *> (pObj->pData);
            auto pHopObjR = Hop_Regular(pHopObj);
            if (Hop_ObjIsConst1(pHopObjR)) {
                #ifdef DEBUG
                assert(Hop_ObjFanin0(pHopObjR) == nullptr);
                assert(Hop_ObjFanin1(pHopObjR) == nullptr);
                #endif
                if (!Hop_IsComplement(pHopObj))
                    ret.second = GetId(pObj);
                else 
                    ret.first = GetId(pObj);
            }
        }
        else
            assert(0);
    }
    return ret;
}


pair <ll, ll> NetMan::CreateConst(bool isVerb) {
    auto consts = GetConstId(isVerb);
    pair <ll, ll> ret(consts);
    if (ret.first == -1) {
        auto pObj = Abc_NtkCreateNodeConst0(GetNet());
        RenameAbcObj(pObj, "const0");
        ret.first = GetId(pObj);
        if (isVerb)
            cout << "create const 0: " << pObj << endl;
    }
    if (ret.second == -1) {
        auto pObj = Abc_NtkCreateNodeConst1(GetNet());
        ret.second = GetId(pObj);
        RenameAbcObj(pObj, "const1");
        if (isVerb)
            cout << "create const 1: " << pObj << endl;
    }
    return ret;
}


void NetMan::MergeConst() {
    pair <ll, ll> ret(-1, -1);
    auto type = GetNetType();
    Abc_Obj_t * pObj = nullptr;
    ll i = 0;
    Abc_NtkForEachNode(GetNet(), pObj, i) {
        if (type == NET_TYPE::GATE || type == NET_TYPE::SOP) {
            if (Abc_NodeIsConst0(pObj)) {
                if (ret.first == -1)
                    ret.first = GetId(pObj);
                else {
                    // cout << "merge const 0: " << pObj << endl;
                    Abc_ObjReplace(pObj, GetObj(ret.first));
                }
            }
            else if (Abc_NodeIsConst1(pObj)) {
                if (ret.second == -1)
                    ret.second = GetId(pObj);
                else {
                    // cout << "merge const 1: " << pObj << endl;
                    Abc_ObjReplace(pObj, GetObj(ret.second));
                }
            }
        }
        else
            assert(0);
    }
}


void NetMan::ReArrInTopoOrd() {
    AbcMan::SetMainNetw(pNtk); // abc manage the memory of the old network
    AbcMan::TopoSort();
    pNtk = Abc_NtkDup(AbcMan::GetNet()); // NetMan manage the memory of the duplicated network
}


vector <Abc_Obj_t * > NetMan::TopoSort() const {
    vector <Abc_Obj_t *> nodes;
    nodes.reserve(GetNodeNum());
    SetNetNotTrav();
    for (ll i = 0; i < GetPoNum(); ++i) {
        auto pDriver = GetFanin(GetPo(i), 0);
        if (!GetObjTrav(pDriver))
            TopoSortRec(pDriver, nodes);
    }
    return nodes;
}


void NetMan::TopoSortRec(Abc_Obj_t * pObj, vector <Abc_Obj_t *> & nodes) const {
    if (!IsNode(pObj))
        return;
    if (IsConst(pObj))
        return;
    SetObjTrav(pObj);
    for (ll i = 0; i < GetFaninNum(pObj); ++i) {
        auto pFanin = GetFanin(pObj, i);
        if (!GetObjTrav(pFanin))
            TopoSortRec(pFanin, nodes);
    }
    nodes.emplace_back(pObj);
}


IntVect NetMan::TopoSortWithIds() const {
    IntVect nodes;
    nodes.reserve(GetNodeNum());
    SetNetNotTrav();
    for (ll i = 0; i < GetPoNum(); ++i) {
        auto pDriver = GetFanin(GetPo(i), 0);
        if (!GetObjTrav(pDriver))
            TopoSortRecWithIds(pDriver, nodes);
    }
    return nodes;
}


void NetMan::TopoSortRecWithIds(Abc_Obj_t* pObj, RETURN_VAR IntVect& nodes) const {
    if (!IsNode(pObj))
        return;
    if (IsConst(pObj))
        return;
    SetObjTrav(pObj);
    for (ll i = 0; i < GetFaninNum(pObj); ++i) {
        auto pFanin = GetFanin(pObj, i);
        if (!GetObjTrav(pFanin))
            TopoSortRecWithIds(pFanin, nodes);
    }
    nodes.emplace_back(GetId(pObj));
}


vector <Abc_Obj_t *> NetMan::GetTFI(abc::Abc_Obj_t * pObj) const {
    vector <Abc_Obj_t *> nodes;
    nodes.reserve(GetNodeNum());
    SetNetNotTrav();
    for (ll i = 0; i < GetFaninNum(pObj); ++i) {
        auto pFanin = GetFanin(pObj, i);
        if (!GetObjTrav(pFanin))
            GetTFIRec(pFanin, nodes);
    }
    return nodes;
}


void NetMan::GetTFIRec(abc::Abc_Obj_t * pObj, std::vector <abc::Abc_Obj_t *> & nodes) const {
    if (!IsNode(pObj))
        return;
    SetObjTrav(pObj);
    for (ll i = 0; i < GetFaninNum(pObj); ++i) {
        auto pFanin = GetFanin(pObj, i);
        if (!GetObjTrav(pFanin))
            GetTFIRec(pFanin, nodes);
    }
    nodes.emplace_back(pObj);
}


vector <ll> NetMan::GetTFI(abc::Abc_Obj_t * pObj, const set <ll> & critGraph) const {
    vector <ll> objs;
    objs.reserve(GetNodeNum());
    SetNetNotTrav();
    for (ll i = 0; i < GetFaninNum(pObj); ++i) {
        auto pFanin = GetFanin(pObj, i);
        if (!GetObjTrav(pFanin))
            GetTFIRec(pFanin, objs, critGraph);
    }
    return objs;
}


void NetMan::GetTFIRec(abc::Abc_Obj_t * pObj, std::vector <ll> & objs, const set <ll> & critGraph) const {
    if (critGraph.count(pObj->Id) == 0)
        return;
    // if (!IsNode(pObj))
    //     return;
    SetObjTrav(pObj);
    for (ll i = 0; i < GetFaninNum(pObj); ++i) {
        auto pFanin = GetFanin(pObj, i);
        if (!GetObjTrav(pFanin))
            GetTFIRec(pFanin, objs, critGraph);
    }
    objs.emplace_back(pObj->Id);
}


vector <Abc_Obj_t *> NetMan::GetTFO(abc::Abc_Obj_t * pObj) const {
    vector <Abc_Obj_t *> nodes;
    nodes.reserve(GetNodeNum());
    SetNetNotTrav();
    for (ll i = 0; i < GetFanoutNum(pObj); ++i) {
        auto pFanout = GetFanout(pObj, i);
        if (!GetObjTrav(pFanout))
            GetTFORec(pFanout, nodes);
    }
    reverse(nodes.begin(), nodes.end());
    return nodes;
}


void NetMan::GetTFORec(abc::Abc_Obj_t * pObj, std::vector <abc::Abc_Obj_t *> & nodes) const {
    if (!IsNode(pObj))
        return;
    SetObjTrav(pObj);
    for (ll i = 0; i < GetFanoutNum(pObj); ++i) {
        auto pFanout = GetFanout(pObj, i);
        if (!GetObjTrav(pFanout))
            GetTFORec(pFanout, nodes);
    }
    nodes.emplace_back(pObj);
}


vector <ll> NetMan::GetTFO(abc::Abc_Obj_t * pObj, const set <ll> & critGraph) const {
    vector <ll> objs;
    objs.reserve(GetNodeNum());
    SetNetNotTrav();
    for (ll i = 0; i < GetFanoutNum(pObj); ++i) {
        auto pFanout = GetFanout(pObj, i);
        if (!GetObjTrav(pFanout))
            GetTFORec(pFanout, objs, critGraph);
    }
    reverse(objs.begin(), objs.end());
    return objs;
}


void NetMan::GetTFORec(abc::Abc_Obj_t * pObj, std::vector <ll> & objs, const set <ll> & critGraph) const {
    if (critGraph.count(pObj->Id) == 0)
        return;
    // if (!IsNode(pObj))
    //     return;
    SetObjTrav(pObj);
    for (ll i = 0; i < GetFanoutNum(pObj); ++i) {
        auto pFanout = GetFanout(pObj, i);
        if (!GetObjTrav(pFanout))
            GetTFORec(pFanout, objs, critGraph);
    }
    objs.emplace_back(pObj->Id);
}


std::vector <abc::Abc_Obj_t *> NetMan::GetFanins(abc::Abc_Obj_t * pObj) const {
    vector <Abc_Obj_t *> nodes;
    nodes.reserve(GetFaninNum(pObj));
    for (ll i = 0; i < GetFaninNum(pObj); ++i)
        nodes.emplace_back(GetFanin(pObj, i));
    return nodes;
}


std::vector <abc::Abc_Obj_t *> NetMan::GetFanouts(abc::Abc_Obj_t * pObj) const {
    vector <Abc_Obj_t *> nodes;
    nodes.reserve(GetFanoutNum(pObj));
    for (ll i = 0; i < GetFanoutNum(pObj); ++i)
        nodes.emplace_back(GetFanout(pObj, i));
    return nodes;
}


Abc_Obj_t* NetMan::GetObjByName(const std::string & name, bool considerPO) const {
    auto pName = const_cast<char*>(name.c_str());
    Abc_Obj_t * pObj;
    int Num;
    // try to find the terminal
    if (considerPO) {
        Num = Nm_ManFindIdByName( pNtk->pManName, pName, ABC_OBJ_PO );
        if ( Num >= 0 )
            return Abc_ObjFanin0( Abc_NtkObj( pNtk, Num ) );
    }
    Num = Nm_ManFindIdByName( pNtk->pManName, pName, ABC_OBJ_PI );
    if ( Num >= 0 )
        // return Abc_ObjFanin0( Abc_NtkObj( pNtk, Num ) );
        return Abc_NtkObj( pNtk, Num );
    Num = Nm_ManFindIdByName( pNtk->pManName, pName, ABC_OBJ_NODE );
    if ( Num >= 0 )
        return Abc_NtkObj( pNtk, Num );
    // find the internal node
    if ( pName[0] != 'n' )
    {
        printf( "Name \"%s\" is not found among CO or node names (internal names often look as \"n<num>\").\n", pName );
        return nullptr;
    }
    Num = atoi( pName + 1 );
    if ( Num < 0 || Num >= Abc_NtkObjNumMax(pNtk) )
    {
        printf( "The node \"%s\" with ID %d is not in the current network.\n", pName, Num );
        return nullptr;
    }
    pObj = Abc_NtkObj( pNtk, Num );
    if ( pObj == nullptr )
    {
        printf( "The node \"%s\" with ID %d has been removed from the current network.\n", pName, Num );
        return nullptr;
    }
    if ( !Abc_ObjIsNode(pObj) )
    {
        printf( "Object with ID %d is not a node.\n", Num );
        return nullptr;
    }
    return pObj;
}


void NetMan::Comm(const std::string & cmd, bool isVerb) {
    assert(isDupl == true);
    AbcMan::SetMainNetw(pNtk); // abc manage the memory of the old network
    AbcMan::Comm(cmd, isVerb);
    pNtk = Abc_NtkDup(AbcMan::GetNet()); // NetMan manage the memory of the duplicated network
}


void NetMan::Sweep() {
    assert(isDupl == true);
    AbcMan::SetMainNetw(pNtk); // abc manage the memory of the old network
    AbcMan::Sweep();
    pNtk = Abc_NtkDup(AbcMan::GetNet()); // NetMan manage the memory of the duplicated network
}


void NetMan::SynthWithResyn2Comm() {
    assert(isDupl == true);
    AbcMan::SetMainNetw(pNtk); // abc manage the memory of the old network
    AbcMan::SynthWithResyn2Comm();
    pNtk = Abc_NtkDup(AbcMan::GetNet()); // NetMan manage the memory of the duplicated network
}


void NetMan::SynthWithCompress2rsComm() {
    assert(isDupl == true);
    AbcMan::SetMainNetw(pNtk); // abc manage the memory of the old network
    AbcMan::Comm("st; compress2rs; logic; sop; ps");
    pNtk = Abc_NtkDup(AbcMan::GetNet()); // NetMan manage the memory of the duplicated network
}


void NetMan::ConvToSopWithTwoInps() {
    #ifdef DEBUG
    assert(isDupl == true);
    #endif
    AbcMan::SetMainNetw(pNtk); // abc manage the memory of the old network
    AbcMan::Comm("st; logic; sop;");
    pNtk = Abc_NtkDup(AbcMan::GetNet()); // NetMan manage the memory of the duplicated network
}


void NetMan::SynthAndMap(double maxDelay, bool isVerb) {
    #ifdef DEBUG
    assert(isDupl == true);
    #endif
    AbcMan::SetMainNetw(pNtk); // abc manage the memory of the old network
    AbcMan::SynthAndMap(maxDelay, isVerb);
    pNtk = Abc_NtkDup(AbcMan::GetNet()); // NetMan manage the memory of the duplicated network
}


void NetMan::Compile(double maxDelay) {
    assert(isDupl == true);
    AbcMan::SetMainNetw(pNtk); // abc manage the memory of the old network

    // synthesize and map
    const int OPT_ROUND = 1;
    ostringstream oss("");
    oss << "st; ";
    for (int i = 0; i < OPT_ROUND; ++i)
        oss << "compress2rs; ps;";
    AbcMan::Comm(oss.str());
    oss.str("");
    oss << maxDelay;

    // backup AIG
    Abc_Ntk_t * pTempNtk = Abc_NtkDup(AbcMan::GetNet());

    // try "dch; map -D"
    // AbcMan::Comm("dch; map -D " + oss.str());
    AbcMan::Comm("dch; amap;");
    double area = AbcMan::GetArea();
    double delay = AbcMan::GetDelay();
    // cout << "using \"dch; map -D " << oss.str() << "\", area = " << area << ", delay = " << delay << endl;
    cout << "using \"dch; amap\", area = " << area << ", delay = " << delay << endl;
    auto pRecNtk = Abc_NtkDup(AbcMan::GetNet());
    double recArea = area, recDelay = delay;

    // try "map -D"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    AbcMan::Comm("map -D " + oss.str());
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    cout << "using \"map -D " << oss.str() << "\", area = " << area << ", delay = " << delay << endl;
    // if delay <= maxDelay
    // or if delay of recorded circuit > maxDelay
    if (DoubleLessEqual(delay, maxDelay) || DoubleGreat(recDelay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    // try "dch; map -a"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    AbcMan::Comm("dch; map -a");
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    cout << "using \"dch; map -a\", area = " << area << ", delay = " << delay << endl;
    if (DoubleLessEqual(delay, maxDelay) || DoubleGreat(recDelay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    // try "map -a"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    AbcMan::Comm("map -a");
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    cout << "using \"map -a\", area = " << area << ", delay = " << delay << endl;
    if (DoubleLessEqual(delay, maxDelay) || DoubleGreat(recDelay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    // try "dch; amap"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    // AbcMan::Comm("dch; amap");
    AbcMan::Comm("dch; map -D " + oss.str());
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    // cout << "using \"dch; amap\", area = " << area << ", delay = " << delay << endl;
    cout << "using \"dch; map -D " << oss.str() << "\", area = " << area << ", delay = " << delay << endl;
    if (DoubleLessEqual(delay, maxDelay) || DoubleGreat(recDelay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    // try "amap"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    AbcMan::Comm("amap");
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    cout << "using \"amap\", area = " << area << ", delay = " << delay << endl;
    if (DoubleLessEqual(delay, maxDelay) || DoubleGreat(recDelay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    Abc_NtkDelete(pTempNtk);
    pNtk = pRecNtk; // NetMan manage the memory of the duplicated network
}


void NetMan::CompileNew(double maxDelay) {
    assert(isDupl == true);
    AbcMan::SetMainNetw(pNtk); // abc manage the memory of the old network

    // synthesize and map
    const int OPT_ROUND = 1;
    ostringstream oss("");
    oss << "st; ";
    for (int i = 0; i < OPT_ROUND; ++i)
        oss << "compress2rs; ps;";
    AbcMan::Comm(oss.str());
    oss.str("");
    oss << maxDelay;

    // backup AIG
    Abc_Ntk_t * pTempNtk = Abc_NtkDup(AbcMan::GetNet());

    // try "dch; map -D"
    AbcMan::Comm("dch; map -D " + oss.str());
    double area = AbcMan::GetArea();
    double delay = AbcMan::GetDelay();
    cout << "using \"dch; map -D " << oss.str() << "\", area = " << area << ", delay = " << delay << endl;
    auto pRecNtk = Abc_NtkDup(AbcMan::GetNet());
    double recArea = area, recDelay = delay;

    // try "map -D"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    AbcMan::Comm("map -D " + oss.str());
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    cout << "using \"map -D " << oss.str() << "\", area = " << area << ", delay = " << delay << endl;
    // if delay <= maxDelay
    // or if delay of recorded circuit > maxDelay
    if (DoubleLessEqual(delay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    // try "dch; map -a"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    AbcMan::Comm("dch; map -a");
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    cout << "using \"dch; map -a\", area = " << area << ", delay = " << delay << endl;
    if (DoubleLessEqual(delay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    // try "map -a"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    AbcMan::Comm("map -a");
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    cout << "using \"map -a\", area = " << area << ", delay = " << delay << endl;
    if (DoubleLessEqual(delay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    // try "dch; amap"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    AbcMan::Comm("dch; amap");
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    cout << "using \"dch; amap\", area = " << area << ", delay = " << delay << endl;
    if (DoubleLessEqual(delay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    // try "amap"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    AbcMan::Comm("amap");
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    cout << "using \"amap\", area = " << area << ", delay = " << delay << endl;
    if (DoubleLessEqual(delay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    Abc_NtkDelete(pTempNtk);
    pNtk = pRecNtk; // NetMan manage the memory of the duplicated network
}


void NetMan::CompileExtremeArea(double maxDelay) {
    assert(isDupl == true);
    AbcMan::SetMainNetw(pNtk); // abc manage the memory of the old network

    // synthesize and map
    // this command is from BLASYS
    ostringstream oss("");
    oss << "strash;ifraig;dc2;fraig;rewrite;refactor;resub;rewrite;refactor;resub;rewrite;rewrite -z;rewrite -z;rewrite -z;balance;refactor -z;refactor -N 11;resub -K 10;resub -K 12;resub -K 14;resub -K 16;refactor;balance;";
    AbcMan::Comm(oss.str());

    // backup AIG
    Abc_Ntk_t * pTempNtk = Abc_NtkDup(AbcMan::GetNet());

    // try "dch; map -D"
    oss.str("");
    oss << maxDelay;
    AbcMan::Comm("dch; map -D " + oss.str());
    double area = AbcMan::GetArea();
    double delay = AbcMan::GetDelay();
    cout << "using \"dch; map -D " << oss.str() << "\", area = " << area << ", delay = " << delay << endl;
    auto pRecNtk = Abc_NtkDup(AbcMan::GetNet());
    double recArea = area, recDelay = delay;

    // try "map -D"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    AbcMan::Comm("map -D " + oss.str());
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    cout << "using \"map -D " << oss.str() << "\", area = " << area << ", delay = " << delay << endl;
    if (DoubleLessEqual(delay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    // try "dch; map -a"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    AbcMan::Comm("dch; map -a");
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    cout << "using \"dch; map -a\", area = " << area << ", delay = " << delay << endl;
    if (DoubleLessEqual(delay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    // try "map -a"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    AbcMan::Comm("map -a");
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    cout << "using \"map -a\", area = " << area << ", delay = " << delay << endl;
    if (DoubleLessEqual(delay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    // try "dch; amap"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    AbcMan::Comm("dch; amap");
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    cout << "using \"dch; amap\", area = " << area << ", delay = " << delay << endl;
    if (DoubleLessEqual(delay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    // try "amap"
    AbcMan::SetMainNetw(Abc_NtkDup(pTempNtk));
    AbcMan::Comm("amap");
    area = AbcMan::GetArea();
    delay = AbcMan::GetDelay();
    cout << "using \"amap\", area = " << area << ", delay = " << delay << endl;
    if (DoubleLessEqual(delay, maxDelay)) {
        if (DoubleLess(area, recArea) || 
            (DoubleEqual(area, recArea) && DoubleLess(delay, recDelay))) {
            Abc_NtkDelete(pRecNtk);
            pRecNtk = Abc_NtkDup(AbcMan::GetNet());
            recArea = area;
            recDelay = delay;
        }
    }

    Abc_NtkDelete(pTempNtk);
    pNtk = pRecNtk; // NetMan manage the memory of the duplicated network
}


void NetMan::CompileFast4Area() {
    assert(isDupl == true);
    AbcMan::SetMainNetw(pNtk); // abc manage the memory of the old network

    // synthesize and map
    const int OPT_ROUND = 1;
    ostringstream oss("");
    oss << "st; ";
    for (int i = 0; i < OPT_ROUND; ++i)
        // oss << "resyn2; ps;";
        oss << "compress2rs; ps;";
    oss << " dch; amap;";
    AbcMan::Comm(oss.str());

    auto area = AbcMan::GetArea();
    auto delay = AbcMan::GetDelay();
    cout << "using \"dch; amap\", area = " << area << ", delay = " << delay << endl;

    pNtk = Abc_NtkDup(AbcMan::GetNet()); // NetMan manage the memory of the duplicated network
}


void NetMan::CompileWithYosys(std::string& standCellPath) {
    // try synthesize and map with yosys
    string aigFileName("tmp/aig.blif");
    string mapFileName("tmp/map.v");
    WriteBlif(aigFileName);

    // for genlib
    cout << standCellPath << endl;
    bool useGenlib = (standCellPath.find(".genlib") != string::npos);
    string abcScript("source,abc.rc;st;compress2rs;dch;amap;");
    if (useGenlib) {
        abcScript += "ps;";
        auto yosysScript = "read_blif " + aigFileName + "; synth -flatten; opt; opt_clean -purge; opt; opt_clean -purge; abc -genlib " + standCellPath + " -script +" + abcScript + "; write_verilog -noattr -nohex " + mapFileName;
        auto yosysComm = "yosys -p \"" + yosysScript + "\" | grep area";
        ExecSystComm(yosysComm);
    }
    else { // for liberty
        assert(standCellPath.find(".lib") != string::npos);
        abcScript += "topo;stime;";
        auto yosysScript = "read_blif " + aigFileName + "; synth -flatten; opt; opt_clean -purge; opt; opt_clean -purge; abc -liberty " + standCellPath + " -script +" + abcScript + "; write_verilog -noattr -nohex " + mapFileName;
        auto yosysComm = "yosys -p \"" + yosysScript + "\" | grep Area";
        ExecSystComm(yosysComm);
    }

    // read the circuit back
    assert(isDupl == true);
    AbcMan::SetMainNetw(pNtk); // abc manage the memory of the old network
    AbcMan::Comm("r -m " + mapFileName);
    if (useGenlib)
        AbcMan::Comm("unbuffer; ps;"); 
    else
        AbcMan::Comm("unbuffer; topo; stime;"); 
    pNtk = Abc_NtkDup(AbcMan::GetNet()); // NetMan manage the memory of the duplicated network
}


void NetMan::WriteBlif(const string & fileName) const {
    cout << "write blif to " << fileName << endl;
    FILE * fp = fopen(fileName.c_str(), "w");
    assert(fp != nullptr);
    fprintf(fp, ".model %s\n", GetNet()->pName);
    // dump inputs
    fprintf(fp, ".inputs ");
    for (int i = 0; i < GetPiNum(); ++i)
        fprintf(fp, "%s ", const_cast<char*>(GetPiName(i).c_str()));
    fprintf(fp, "\n");
    // dump outputs
    fprintf(fp, ".outputs ");
    unordered_map<string, int> isPoPrinted;
    for (int i = 0; i < GetPoNum(); ++i) {
        fprintf(fp, "%s ", const_cast<char*>(GetPoName(i).c_str()));
        isPoPrinted[GetPoName(i)] = 0;
    }
    fprintf(fp, "\n");
    // dump nodes
    auto netType = GetNetType();
    assert(netType == NET_TYPE::SOP || netType == NET_TYPE::GATE);
    for (int id = 0; id < GetIdMaxPlus1(); ++id) {
        auto pObj = GetObj(id);
        if (IsNode(pObj)) {
            fprintf(fp, ".names ");
            for (int i = 0; i < GetFaninNum(pObj); ++i)
                fprintf(fp, "%s ", const_cast<char*>(GetName(GetFanin(pObj, i)).c_str()));
            fprintf(fp, "%s\n", const_cast<char*>(GetName(pObj).c_str()));
            if (netType == NET_TYPE::SOP)
                fprintf(fp, "%s\n", static_cast <char *> (pObj->pData));
            else if (netType == NET_TYPE::GATE) {
                fprintf(fp, "# %s\n", Mio_GateReadName((Mio_Gate_t *)pObj->pData));
                fprintf(fp, "%s\n", Mio_GateReadSop((Mio_Gate_t *)pObj->pData));
            }
            else
                assert(0);
            if (isPoPrinted.count(GetName(pObj)))
                isPoPrinted[GetName(pObj)] = 1;
        }
    }
    for (int id = 0; id < GetPoNum(); ++id) {
        auto pObj = GetPo(id);
        if (isPoPrinted[GetName(pObj)] == 0) {
            fprintf(fp, ".names ");
            fprintf(fp, "%s ", const_cast<char*>(GetName(GetFanin(pObj, 0)).c_str()));
            fprintf(fp, "%s\n", const_cast<char*>(GetName(pObj).c_str()));
            assert(Abc_ObjIsComplement(pObj) == 0);
            fprintf(fp, "1 1\n");
        }
    }

    fprintf(fp, ".end\n");
    fclose(fp);
}


static char * Abc_NtkPrintSop( char * pSop ) {
    static char Buffer[1000];
    char * pGet, * pSet;
    pSet = Buffer;
    for ( pGet = pSop; *pGet; pGet++ )
    {        
        if ( *pGet == '\n' )
        {
            *pSet++ = '\\';
            *pSet++ = 'n';
        }
        else
            *pSet++ = *pGet;
    }
    *(pSet-2) = 0;
    return Buffer;
}
static int Abc_NtkCountLogicNodes( Vec_Ptr_t * vNodes ) {
    Abc_Obj_t * pObj;
    int i, Counter = 0;
    Vec_PtrForEachEntry( Abc_Obj_t *, vNodes, pObj, i )
    {
        if ( !Abc_ObjIsNode(pObj) )
            continue;
        if ( Abc_ObjFaninNum(pObj) == 0 && Abc_ObjFanoutNum(pObj) == 0 )
            continue;
        Counter ++;
    }
    return Counter;
}
void Net_WriteDotNtk( Abc_Ntk_t * pNtk, Vec_Ptr_t * vNodes, Vec_Ptr_t * vNodesShow, const char * pFileName, int fGateNames, int fUseReverse, vector<double> * pData, double threshold ) {
    FILE * pFile;
    Abc_Obj_t * pNode, * pFanin;
    char * pSopString;
    int LevelMin, LevelMax, fHasCos, Level, i, k, fHasBdds, fCompl, Prev;
    // int Limit = 500;
    int Limit = 1500;

    assert( Abc_NtkIsStrash(pNtk) || Abc_NtkIsLogic(pNtk) );

    if ( vNodes->nSize < 1 )
    {
        printf( "The set has no nodes. DOT file is not written.\n" );
        return;
    }

    if ( vNodes->nSize > Limit )
    {
        printf( "The set has more than %d nodes. DOT file is not written.\n", Limit );
        return;
    }

    // start the stream
    if ( (pFile = fopen( pFileName, "w" )) == NULL )
    {
        fprintf( stdout, "Cannot open the intermediate file \"%s\".\n", pFileName );
        return;
    }

    // transform logic functions from BDD to SOP
    if ( (fHasBdds = Abc_NtkIsBddLogic(pNtk)) )
    {
        if ( !Abc_NtkBddToSop(pNtk, -1, ABC_INFINITY, 1) )
        {
            printf( "Io_WriteDotNtk(): Converting to SOPs has failed.\n" );
            return;
        }
    }

    // mark the nodes from the set
    Vec_PtrForEachEntry( Abc_Obj_t *, vNodes, pNode, i )
        pNode->fMarkC = 1;
    if ( vNodesShow )
        Vec_PtrForEachEntry( Abc_Obj_t *, vNodesShow, pNode, i )
            pNode->fMarkB = 1;

    // get the levels of nodes
    LevelMax = Abc_NtkLevel( pNtk );
    if ( fUseReverse )
    {
        LevelMin = Abc_NtkLevelReverse( pNtk );
        assert( LevelMax == LevelMin );
        Vec_PtrForEachEntry( Abc_Obj_t *, vNodes, pNode, i )
            if ( Abc_ObjIsNode(pNode) )
                pNode->Level = LevelMax - pNode->Level + 1;
    }

    // find the largest and the smallest levels
    LevelMin = 10000;
    LevelMax = -1;
    fHasCos  = 0;
    Vec_PtrForEachEntry( Abc_Obj_t *, vNodes, pNode, i )
    {
        if ( Abc_ObjIsCo(pNode) )
        {
            fHasCos = 1;
            continue;
        }
        if ( LevelMin > (int)pNode->Level )
            LevelMin = pNode->Level;
        if ( LevelMax < (int)pNode->Level )
            LevelMax = pNode->Level;
    }

    // set the level of the CO nodes
    if ( fHasCos )
    {
        LevelMax++;
        Vec_PtrForEachEntry( Abc_Obj_t *, vNodes, pNode, i )
        {
            if ( Abc_ObjIsCo(pNode) )
                pNode->Level = LevelMax;
        }
    }

    // write the DOT header
    fprintf( pFile, "# %s\n",  "Network structure generated by ABC" );
    fprintf( pFile, "\n" );
    fprintf( pFile, "digraph network {\n" );
    fprintf( pFile, "size = \"7.5,10\";\n" );
//    fprintf( pFile, "size = \"10,8.5\";\n" );
//    fprintf( pFile, "size = \"14,11\";\n" );
//    fprintf( pFile, "page = \"8,11\";\n" );
//  fprintf( pFile, "ranksep = 0.5;\n" );
//  fprintf( pFile, "nodesep = 0.5;\n" );
    fprintf( pFile, "center = true;\n" );
//    fprintf( pFile, "orientation = landscape;\n" );
//  fprintf( pFile, "edge [fontsize = 10];\n" );
//  fprintf( pFile, "edge [dir = none];\n" );
    fprintf( pFile, "edge [dir = back];\n" );
    fprintf( pFile, "\n" );

    // labels on the left of the picture
    fprintf( pFile, "{\n" );
    fprintf( pFile, "  node [shape = plaintext];\n" );
    fprintf( pFile, "  edge [style = invis];\n" );
    fprintf( pFile, "  LevelTitle1 [label=\"\"];\n" );
    fprintf( pFile, "  LevelTitle2 [label=\"\"];\n" );
    // generate node names with labels
    for ( Level = LevelMax; Level >= LevelMin; Level-- )
    {
        // the visible node name
        fprintf( pFile, "  Level%d", Level );
        fprintf( pFile, " [label = " );
        // label name
        fprintf( pFile, "\"" );
        fprintf( pFile, "\"" );
        fprintf( pFile, "];\n" );
    }

    // genetate the sequence of visible/invisible nodes to mark levels
    fprintf( pFile, "  LevelTitle1 ->  LevelTitle2 ->" );
    for ( Level = LevelMax; Level >= LevelMin; Level-- )
    {
        // the visible node name
        fprintf( pFile, "  Level%d",  Level );
        // the connector
        if ( Level != LevelMin )
            fprintf( pFile, " ->" );
        else
            fprintf( pFile, ";" );
    }
    fprintf( pFile, "\n" );
    fprintf( pFile, "}" );
    fprintf( pFile, "\n" );
    fprintf( pFile, "\n" );

    // generate title box on top
    fprintf( pFile, "{\n" );
    fprintf( pFile, "  rank = same;\n" );
    fprintf( pFile, "  LevelTitle1;\n" );
    fprintf( pFile, "  title1 [shape=plaintext,\n" );
    fprintf( pFile, "          fontsize=20,\n" );
    fprintf( pFile, "          fontname = \"Times-Roman\",\n" );
    fprintf( pFile, "          label=\"" );
    fprintf( pFile, "%s", "Network structure visualized by ABC" );
    fprintf( pFile, "\\n" );
    fprintf( pFile, "Benchmark \\\"%s\\\". ", pNtk->pName );
    fprintf( pFile, "Time was %s. ",  Extra_TimeStamp() );
    fprintf( pFile, "\"\n" );
    fprintf( pFile, "         ];\n" );
    fprintf( pFile, "}" );
    fprintf( pFile, "\n" );
    fprintf( pFile, "\n" );

    // generate statistics box
    fprintf( pFile, "{\n" );
    fprintf( pFile, "  rank = same;\n" );
    fprintf( pFile, "  LevelTitle2;\n" );
    fprintf( pFile, "  title2 [shape=plaintext,\n" );
    fprintf( pFile, "          fontsize=18,\n" );
    fprintf( pFile, "          fontname = \"Times-Roman\",\n" );
    fprintf( pFile, "          label=\"" );
    if ( Abc_NtkObjNum(pNtk) == Vec_PtrSize(vNodes) )
        fprintf( pFile, "The network contains %d logic nodes and %d latches.", Abc_NtkNodeNum(pNtk), Abc_NtkLatchNum(pNtk) );
    else
        fprintf( pFile, "The set contains %d logic nodes and spans %d levels.", Abc_NtkCountLogicNodes(vNodes), LevelMax - LevelMin + 1 );
    fprintf( pFile, "\\n" );
    fprintf( pFile, "\"\n" );
    fprintf( pFile, "         ];\n" );
    fprintf( pFile, "}" );
    fprintf( pFile, "\n" );
    fprintf( pFile, "\n" );

    // generate the POs
    if ( fHasCos )
    {
        fprintf( pFile, "{\n" );
        fprintf( pFile, "  rank = same;\n" );
        // the labeling node of this level
        fprintf( pFile, "  Level%d;\n",  LevelMax );
        // generate the PO nodes
        Vec_PtrForEachEntry( Abc_Obj_t *, vNodes, pNode, i )
        {
            if ( !Abc_ObjIsCo(pNode) )
                continue;
            if (pData == nullptr) {
                fprintf( pFile, "  Node%d [label = \"%s%s\"", 
                    pNode->Id, 
                    (Abc_ObjIsBi(pNode)? Abc_ObjName(Abc_ObjFanout0(pNode)):Abc_ObjName(pNode)), 
                    (Abc_ObjIsBi(pNode)? "_in":"") );
            }
            else {
                fprintf( pFile, "  Node%d [label = \"%s%s\\n%lf\"", 
                    pNode->Id, 
                    (Abc_ObjIsBi(pNode)? Abc_ObjName(Abc_ObjFanout0(pNode)):Abc_ObjName(pNode)), 
                    (Abc_ObjIsBi(pNode)? "_in":""),
                    (*pData)[pNode->Id] );
            }
            fprintf( pFile, ", shape = %s", (Abc_ObjIsBi(pNode)? "box":"invtriangle") );
            if (pData == nullptr) {
                if ( pNode->fMarkB )
                    fprintf( pFile, ", style = filled" );
            }
            else {
                if ( pNode->fMarkB || (*pData)[pNode->Id] < threshold || (*pData)[pNode->Id] > 1.0 - threshold )
                    fprintf( pFile, ", style = filled" );
            }
            fprintf( pFile, ", color = coral, fillcolor = coral" );
            fprintf( pFile, "];\n" );
        }
        fprintf( pFile, "}" );
        fprintf( pFile, "\n" );
        fprintf( pFile, "\n" );
    }

    // generate nodes of each rank
    for ( Level = LevelMax - fHasCos; Level >= LevelMin && Level > 0; Level-- )
    {
        fprintf( pFile, "{\n" );
        fprintf( pFile, "  rank = same;\n" );
        // the labeling node of this level
        fprintf( pFile, "  Level%d;\n",  Level );
        Vec_PtrForEachEntry( Abc_Obj_t *, vNodes, pNode, i )
        {
            if ( (int)pNode->Level != Level )
                continue;
            if ( Abc_ObjFaninNum(pNode) == 0 )
                continue;

            if ( Abc_NtkIsStrash(pNtk) )
                pSopString = "";
            else if ( Abc_NtkHasMapping(pNtk) && fGateNames )
                pSopString = Mio_GateReadName((Mio_Gate_t *)pNode->pData);
            else if ( Abc_NtkHasMapping(pNtk) )
                pSopString = Abc_NtkPrintSop(Mio_GateReadSop((Mio_Gate_t *)pNode->pData));
            else
                pSopString = Abc_NtkPrintSop((char *)pNode->pData);
            if (pData == nullptr)
                // fprintf( pFile, "  Node%d [label = \"%d\\n%s\"", pNode->Id, pNode->Id, pSopString );
                fprintf( pFile, "  Node%d [label = \"%s\\n%s\"", pNode->Id, Abc_ObjName(pNode), pSopString );
            else
                fprintf( pFile, "  Node%d [label = \"%s\\n%s\\n%lf\"", pNode->Id, Abc_ObjName(pNode), pSopString, (*pData)[pNode->Id] );

            fprintf( pFile, ", shape = ellipse" );
            if (pData == nullptr) {
                if ( pNode->fMarkB )
                    fprintf( pFile, ", style = filled" );
            }
            else {
                if ( pNode->fMarkB || (*pData)[pNode->Id] < threshold || (*pData)[pNode->Id] > 1.0 - threshold )
                    fprintf( pFile, ", style = filled, color = lightblue" );
            }
            fprintf( pFile, "];\n" );
        }
        fprintf( pFile, "}" );
        fprintf( pFile, "\n" );
        fprintf( pFile, "\n" );
    }

    // generate the PI nodes if any
    if ( LevelMin == 0 )
    {
        fprintf( pFile, "{\n" );
        fprintf( pFile, "  rank = same;\n" );
        // the labeling node of this level
        fprintf( pFile, "  Level%d;\n",  LevelMin );
        // generate the PO nodes
        Vec_PtrForEachEntry( Abc_Obj_t *, vNodes, pNode, i )
        {
            if ( !Abc_ObjIsCi(pNode) )
            {
                // check if the costant node is present
                if ( Abc_ObjFaninNum(pNode) == 0 && Abc_ObjFanoutNum(pNode) > 0 )
                {
                    fprintf( pFile, "  Node%d [label = \"Const%d\"", pNode->Id, Abc_NtkIsStrash(pNode->pNtk) || Abc_NodeIsConst1(pNode) );
                    fprintf( pFile, ", shape = ellipse" );
                    if ( pNode->fMarkB )
                        fprintf( pFile, ", style = filled" );
                    fprintf( pFile, ", color = coral, fillcolor = coral" );
                    fprintf( pFile, "];\n" );
                }
                continue;
            }
            if (pData == nullptr) {
                fprintf( pFile, "  Node%d [label = \"%s\"", 
                    pNode->Id, 
                    (Abc_ObjIsBo(pNode)? Abc_ObjName(Abc_ObjFanin0(pNode)):Abc_ObjName(pNode)) );
            }
            else {
                fprintf( pFile, "  Node%d [label = \"%s\\n%lf\"", 
                    pNode->Id, 
                    (Abc_ObjIsBo(pNode)? Abc_ObjName(Abc_ObjFanin0(pNode)):Abc_ObjName(pNode)),
                    (*pData)[pNode->Id] );
            }
            fprintf( pFile, ", shape = %s", (Abc_ObjIsBo(pNode)? "box":"triangle") );
            if (pData == nullptr) {
                if ( pNode->fMarkB )
                    fprintf( pFile, ", style = filled" );
            }
            else {
                if ( pNode->fMarkB || (*pData)[pNode->Id] < threshold || (*pData)[pNode->Id] > 1.0 - threshold )
                    fprintf( pFile, ", style = filled" );
            }
            fprintf( pFile, ", color = coral, fillcolor = coral" );
            fprintf( pFile, "];\n" );
        }
        fprintf( pFile, "}" );
        fprintf( pFile, "\n" );
        fprintf( pFile, "\n" );
    }

    // generate invisible edges from the square down
    fprintf( pFile, "title1 -> title2 [style = invis];\n" );
    Vec_PtrForEachEntry( Abc_Obj_t *, vNodes, pNode, i )
    {
        if ( (int)pNode->Level != LevelMax )
            continue;
        fprintf( pFile, "title2 -> Node%d [style = invis];\n", pNode->Id );
    }
    // generate invisible edges among the COs
    Prev = -1;
    Vec_PtrForEachEntry( Abc_Obj_t *, vNodes, pNode, i )
    {
        if ( (int)pNode->Level != LevelMax )
            continue;
        if ( !Abc_ObjIsPo(pNode) )
            continue;
        if ( Prev >= 0 )
            fprintf( pFile, "Node%d -> Node%d [style = invis];\n", Prev, pNode->Id );
        Prev = pNode->Id;
    }

    // generate edges
    Vec_PtrForEachEntry( Abc_Obj_t *, vNodes, pNode, i )
    {
        if ( Abc_ObjIsLatch(pNode) )
            continue;
        Abc_ObjForEachFanin( pNode, pFanin, k )
        {
            if ( Abc_ObjIsLatch(pFanin) )
                continue;
            fCompl = 0;
            if ( Abc_NtkIsStrash(pNtk) )
                fCompl = Abc_ObjFaninC(pNode, k);
            // generate the edge from this node to the next
            fprintf( pFile, "Node%d",  pNode->Id );
            fprintf( pFile, " -> " );
            fprintf( pFile, "Node%d",  pFanin->Id );
            fprintf( pFile, " [style = %s", fCompl? "dotted" : "solid" );
//            fprintf( pFile, ", label = \"%c\"", 'a' + k );
            fprintf( pFile, "]" );
            fprintf( pFile, ";\n" );
        }
    }

    fprintf( pFile, "}" );
    fprintf( pFile, "\n" );
    fprintf( pFile, "\n" );
    fclose( pFile );

    // unmark the nodes from the set
    Vec_PtrForEachEntry( Abc_Obj_t *, vNodes, pNode, i )
        pNode->fMarkC = 0;
    if ( vNodesShow )
        Vec_PtrForEachEntry( Abc_Obj_t *, vNodesShow, pNode, i )
            pNode->fMarkB = 0;

    // convert the network back into BDDs if this is how it was
    if ( fHasBdds )
        Abc_NtkSopToBdd(pNtk);
}


void NetMan::WriteDot(const string & fileName, std::vector<double>* pData, double threshold) const {
    Vec_Ptr_t * vNodes;
    vNodes = Abc_NtkCollectObjects( pNtk );
    if (pData != nullptr)
        assert(pData->size() == GetIdMaxPlus1());
    Net_WriteDotNtk( pNtk, vNodes, nullptr, fileName.c_str(), 0, 0, pData, threshold );
    Vec_PtrFree( vNodes );
}


void NetMan::Print(bool showFunct) const {
    if (GetNet()->pName != nullptr)
        cout << GetNet()->pName << endl;
    else
        cout << "network" << endl;
    for (ll i = 0; i < GetIdMaxPlus1(); ++i) {
        if (GetObj(i) == nullptr)
            continue;
        PrintObj(i, showFunct);
    }
}


void NetMan::PrintObjBas(Abc_Obj_t * pObj, string && endWith) const {
    cout << GetName(pObj) << "(" << GetId(pObj) << ")" << endWith;
}


void NetMan::PrintObj(Abc_Obj_t * pObj, bool showFunct) const {
    PrintObjBas(pObj, ":");
    for (ll i = 0; i < GetFaninNum(pObj); ++i)
        PrintObjBas(GetFanin(pObj, i), ",");
    // cout << endl;
    if (showFunct) {
        if (NetMan::GetNetType() == NET_TYPE::SOP) {
            if (IsNode(pObj))
                cout << static_cast <char *> (pObj->pData);
            else
                cout << endl;
        }
        else if (NetMan::GetNetType() == NET_TYPE::GATE) {
            if (IsNode(pObj))
                cout << Mio_GateReadName(static_cast <Mio_Gate_t *> (pObj->pData)) << endl;
            else
                cout << endl;
        }
        else
            assert(0);
    }
}


bool NetMan::IsPIOSame(NetMan & oth_net) const {
    if (this->GetPiNum() != oth_net.GetPiNum())
        return false;
    for (ll i = 0; i < this->GetPiNum(); ++i) {
        if (this->GetPiName(i) != oth_net.GetPiName(i))
            return false;
    }
    if (this->GetPoNum() != oth_net.GetPoNum())
        return false;
    for (ll i = 0; i < this->GetPoNum(); ++i) {
        if (this->GetPoName(i) != oth_net.GetPoName(i))
            return false;
    }
    return true;
}



ll NetMan::GetNodeMffcSize(Abc_Obj_t * pNode) const {
    assert(IsNode(pNode));
    Vec_Ptr_t * vCone = Vec_PtrAlloc( 100 );
    Abc_NodeDeref_rec(pNode);
    Abc_NodeMffcConeSupp(pNode, vCone, nullptr);
    Abc_NodeRef_rec( pNode );
    ll ret = Vec_PtrSize(vCone);
    Vec_PtrFree(vCone);
    return ret;
}


vector <Abc_Obj_t *> NetMan::GetNodeMffc(Abc_Obj_t * pNode) const {
    assert(IsNode(pNode));
    Vec_Ptr_t * vCone = Vec_PtrAlloc( 100 );
    Abc_NodeDeref_rec(pNode);
    Abc_NodeMffcConeSupp(pNode, vCone, nullptr);
    Abc_NodeRef_rec( pNode );
    vector <Abc_Obj_t *> mffc;
    mffc.reserve(Vec_PtrSize(vCone));
    Abc_Obj_t * pObj = nullptr;
    ll i = 0;
    Vec_PtrForEachEntry(Abc_Obj_t *, vCone, pObj, i)
        mffc.emplace_back(pObj);
    Vec_PtrFree(vCone);
    return mffc;
}


int NetMan::NodeDeref_rec(abc::Abc_Obj_t* pRoot, abc::Abc_Obj_t* pNode, unordered_set<abc::Abc_Obj_t*>& divSet) const {
    Abc_Obj_t * pFanin;
    int i, Counter = 1;
    if ( Abc_ObjIsCi(pNode) || divSet.count(pNode) || (pRoot != pNode && IsPoDriver(pNode)) )
        return 0;
    Abc_ObjForEachFanin( pNode, pFanin, i )
    {
        assert( pFanin->vFanouts.nSize > 0 );
        if ( --pFanin->vFanouts.nSize == 0 )
            Counter += NodeDeref_rec(pRoot, pFanin, divSet);
    }
    return Counter;
}


int NetMan::NodeRef_rec(abc::Abc_Obj_t* pRoot, abc::Abc_Obj_t * pNode, unordered_set<abc::Abc_Obj_t*>& divSet) const {
    Abc_Obj_t * pFanin;
    int i, Counter = 1;
    if ( Abc_ObjIsCi(pNode) || divSet.count(pNode) || (pRoot != pNode && IsPoDriver(pNode)) )
        return 0;
    Abc_ObjForEachFanin( pNode, pFanin, i )
    {
        if ( pFanin->vFanouts.nSize++ == 0 )
            Counter += NodeRef_rec(pRoot, pFanin, divSet);
    }
    return Counter;
}


int NetMan::GetSizeGain(int rootId, const IntVect& divIds) const {
    // get divisor set
    AbcObjSet divSet; 
    for (int divId: divIds)
        divSet.emplace(GetObj(divId));

    // get size gain
    auto pRoot = GetObj(rootId);
    assert(IsNode(pRoot));
    int nSizeGain = NodeDeref_rec(pRoot, pRoot, divSet);
    NodeRef_rec(pRoot, pRoot, divSet);
    return nSizeGain;
}


int NetMan::GetSizeGain(const IntVect& targetIds, const IntVect& divIds) const {
    // get divisor set
    AbcObjSet divSet; 
    for (int divId: divIds)
        divSet.emplace(GetObj(divId));
    
    // get size gain
    int nSizeGain = 0;
    IntSet skipNodes;
    for (int targetId: targetIds) {
        auto pTarget = GetObj(targetId);
        assert(IsNode(pTarget));
        if (pTarget->vFanouts.nSize == 0) {
            skipNodes.emplace(targetId);
            continue;
        }
        nSizeGain += NodeDeref_rec(pTarget, pTarget, divSet);
    }

    int nSizeGain2 = 0;
    for (auto iter = targetIds.rbegin(); iter != targetIds.rend(); ++iter) {
        auto pTarget = GetObj(*iter);
        assert(IsNode(pTarget));
        if (skipNodes.count(*iter))
            continue;
        nSizeGain2 += NodeRef_rec(pTarget, pTarget, divSet);
    }
    assert(nSizeGain == nSizeGain2);
    
    return nSizeGain;
}


abc::Abc_Obj_t* NetMan::CreateNode(const AbcObjVect& pFanins, const std::string& sop) {
    auto pNewNode = abc::Abc_NtkCreateNode(GetNet());
    for (const auto pFanin: pFanins)
        Abc_ObjAddFanin(pNewNode, pFanin);
    assert(GetNetType() == NET_TYPE::SOP);
    pNewNode->pData = Abc_SopRegister((Mem_Flex_t *)GetNet()->pManFunc, sop.c_str());
    return pNewNode;
}


int NetMan::CreateNode(const IntVect& faninIds, const std::string & sop) {
    auto pNewNode = abc::Abc_NtkCreateNode(GetNet());
    for (const auto & faninId: faninIds)
        Abc_ObjAddFanin(pNewNode, GetObj(faninId));
    assert(GetNetType() == NET_TYPE::SOP);
    pNewNode->pData = Abc_SopRegister((Mem_Flex_t *)GetNet()->pManFunc, sop.c_str());
    return pNewNode->Id;
}


int NetMan::CreateNodeAIG(const IntVect& faninIds, const std::string & sop) {
    assert(GetNetType() == NET_TYPE::SOP);
    if (sop == "01 1\n10 1\n" || sop == "10 1\n01 1\n" || sop == "00 0\n11 0\n" || sop == "11 0\n00 0\n") { // xor
        auto and0 = CreateNode(faninIds, "01 1\n");
        auto and1 = CreateNode(faninIds, "10 1\n");
        cout << "create xor" << endl;
        return CreateNode({and0, and1}, "00 0\n");
    }
    else if (sop == "01 0\n10 0\n" || sop == "10 0\n01 0\n" || sop == "00 1\n11 1\n" || sop == "11 1\n00 1\n") { // xnor
        auto and0 = CreateNode(faninIds, "00 1\n");
        auto and1 = CreateNode(faninIds, "11 1\n");
        cout << "create xnor" << endl;
        return CreateNode({and0, and1}, "00 0\n");
    }
    else if (faninIds.size() == 3) {
        auto v0 = faninIds[0], v1 = faninIds[1], v2 = faninIds[2];
        auto pSop = const_cast<char*>(sop.c_str());
        bool isComplement = abc::Abc_SopIsComplement(pSop);
        int nCube = abc::Abc_SopGetCubeNum(pSop);
        int nVar = abc::Abc_SopGetVarNum(pSop);
        assert(nVar == 3);
        if (nCube == 1) {
            auto pCube = pSop;
            assert(pCube[0] != '-' && pCube[1] != '-' && pCube[2] != '-');
            auto and0 = CreateNode({v0, v1}, to_string(pCube[0] - '0') + to_string(pCube[1] - '0') + " 1\n");
            return CreateNode({and0, v2}, "1" + to_string(pCube[2] - '0') + (isComplement? " 0\n": " 1\n"));
        }
        else if (nCube == 2) {
            auto pCube0 = pSop;
            vector<ll> cv0;
            for (int i = 0; i < nVar; ++i) {
                if (pCube0[i] != '-')
                    cv0.emplace_back(i);
            }
            auto pCube1 = pSop + nVar + 3;
            vector<ll> cv1;
            for (int i = 0; i < nVar; ++i) {
                if (pCube1[i] != '-')
                    cv1.emplace_back(i);
            }
            // cout << "cv0.size() = " << cv0.size() << ", cv1.size() = " << cv1.size() << endl;
            if (cv0.size() == 1 && cv1.size() == 2) {
                string phaseVar0((pCube0[cv0[0]] == '0')? "1": "0");
                // cout << "phaseVar0 = " << phaseVar0 << endl;
                auto and0 = CreateNode({faninIds[cv1[0]], faninIds[cv1[1]]}, to_string(pCube1[cv1[0]] - '0') + to_string(pCube1[cv1[1]] - '0') + " 1\n");
                // cout << "sop = " << to_string(pCube1[cv1[0]] - '0') + to_string(pCube1[cv1[1]] - '0') + " 1\n" << endl;
                // cout << "sop = " << "0" + phaseVar0 + (isComplement? " 1\n": " 0\n") << endl;
                return CreateNode({and0, faninIds[cv0[0]]}, "0" + phaseVar0 + (isComplement? " 1\n": " 0\n"));
            }
            else if (cv0.size() == 2 && cv1.size() == 1) {
                string phaseVar1((pCube1[cv1[0]] == '0')? "1": "0");
                auto and0 = CreateNode({faninIds[cv0[0]], faninIds[cv0[1]]}, to_string(pCube0[cv0[0]] - '0') + to_string(pCube0[cv0[1]] - '0') + " 1\n");
                // cout << "sop = " << to_string(pCube0[cv0[0]] - '0') + to_string(pCube0[cv0[1]] - '0') + " 1\n" << endl;
                // cout << "sop = " << "0" + phaseVar1 + (isComplement? " 1\n": " 0\n") << endl;
                return CreateNode({and0, faninIds[cv1[0]]}, "0" + phaseVar1 + (isComplement? " 1\n": " 0\n"));
            }
            else
                assert(0);
        }
        else
            assert(0);
    }
    else {
        auto pNewNode = abc::Abc_NtkCreateNode(GetNet());
        for (const auto & faninId: faninIds)
            Abc_ObjAddFanin(pNewNode, GetObj(faninId));
        pNewNode->pData = Abc_SopRegister((Mem_Flex_t *)GetNet()->pManFunc, sop.c_str());
        return pNewNode->Id;
    }
    return 0;
}


// std::vector <ll> NetMan::TempRepl(abc::Abc_Obj_t * pTS, abc::Abc_Obj_t * pSS) {
//     #ifdef DEBUG
//     assert(pTS != pSS);
//     assert(pTS->pNtk == pSS->pNtk);
//     assert(abc::Abc_ObjFanoutNum(pTS));
//     #endif
//     // record fanouts
//     vector <ll> ret = {pTS->Id, pSS->Id};
//     Abc_Obj_t * pFanout = nullptr;
//     ll i = 0;
//     Abc_ObjForEachFanout(pTS, pFanout, i) {
//         ret.emplace_back(pFanout->Id);
//         ll iFanin = Vec_IntFind(&pFanout->vFanins, pTS->Id);
//         assert(iFanin != -1);
//         ret.emplace_back(iFanin);
//     }
//     PrintVect(ret, "\n");
//     // transfer fanouts
//     abc::Abc_ObjTransferFanout(pTS, pSS);
//     return ret;
// }


static inline int Vec_IntFindFrom(Vec_Int_t * p, int Entry, int start) {
    int i = 0;
    for ( i = start; i < p->nSize; i++ )
        if ( p->pArray[i] == Entry )
            return i;
    assert(0);
    return -1;
}

IntVect NetMan::TempRepl(abc::Abc_Obj_t * pTS, abc::Abc_Obj_t * pSS) {
    assert(pTS != pSS);
    assert(pTS->pNtk == pSS->pNtk);
    assert(abc::Abc_ObjFanoutNum(pTS));
    // record fanouts
    IntVect ret = {pTS->Id, pSS->Id};
    Abc_Obj_t * pFanout = nullptr;
    int i = 0;
    set<pair<int, int>> foIFaninPair;
    Abc_ObjForEachFanout(pTS, pFanout, i) {
        ret.emplace_back(pFanout->Id);
        int start = 0;
        int iFanin = Vec_IntFindFrom(&pFanout->vFanins, pTS->Id, start);
        while (foIFaninPair.count(pair(pFanout->Id, iFanin)))
            iFanin = Vec_IntFindFrom(&pFanout->vFanins, pTS->Id, iFanin + 1);
        ret.emplace_back(iFanin);
        foIFaninPair.emplace(pair(pFanout->Id, iFanin));
    }
    // PrintVect(ret, "\n");
    // transfer fanouts
    abc::Abc_ObjTransferFanout(pTS, pSS);
    return ret;
}


void NetMan::Recov(IntVect& replTrace, bool isVerb) {
    assert(replTrace.size() > 2);
    assert(replTrace.size() % 2 == 0);
    auto pTS = GetObj(replTrace[0]), pSS = GetObj(replTrace[1]);
    if (isVerb) cout << "recover " << pTS << " and " << pSS << endl;
    for (int i = 1; i < replTrace.size() / 2; ++i) {
        auto pFanout = GetObj(replTrace[i * 2]);
        auto iFanin = replTrace[i * 2 + 1];
        PatchFanin(pFanout, iFanin, pSS, pTS);
    }
}


static inline int Vec_IntFindRev( Vec_Int_t * p, int Entry ) {
    int i;
    // for ( i = 0; i < p->nSize; i++ )
    for (i = p->nSize - 1; i >= 0; --i)
        if ( p->pArray[i] == Entry )
            return i;
    return -1;
}


static inline int Vec_IntRemoveRev( Vec_Int_t * p, int Entry ) {
    int i;
    // for ( i = 0; i < p->nSize; i++ )
    for (i = p->nSize - 1; i >= 0; --i)
        if ( p->pArray[i] == Entry )
            break;
    if ( i == p->nSize )
        return 0;
    assert( i < p->nSize );
    for ( i++; i < p->nSize; i++ )
        p->pArray[i-1] = p->pArray[i];
    p->nSize--;
    return 1;
}


static inline void Vec_IntPushMem( Mem_Step_t * pMemMan, Vec_Int_t * p, int Entry ) {
    if ( p->nSize == p->nCap )
    {
        int * pArray;
        int i;

        if ( p->nSize == 0 )
            p->nCap = 1;
        if ( pMemMan )
            pArray = (int *)Mem_StepEntryFetch( pMemMan, p->nCap * 8 );
        else
            pArray = ABC_ALLOC( int, p->nCap * 2 );
        if ( p->pArray )
        {
            for ( i = 0; i < p->nSize; i++ )
                pArray[i] = p->pArray[i];
            if ( pMemMan )
                Mem_StepEntryRecycle( pMemMan, (char *)p->pArray, p->nCap * 4 );
            else
                ABC_FREE( p->pArray );
        }
        p->nCap *= 2;
        p->pArray = pArray;
    }
    p->pArray[p->nSize++] = Entry;
}


void NetMan::PatchFanin( Abc_Obj_t * pObj, ll iFanin, Abc_Obj_t * pFaninOld, Abc_Obj_t * pFaninNew ) {
    Abc_Obj_t * pFaninNewR = Abc_ObjRegular(pFaninNew);
    assert( !Abc_ObjIsComplement(pObj) );
    assert( !Abc_ObjIsComplement(pFaninOld) );
    assert( pFaninOld != pFaninNewR );
    assert( pObj->pNtk == pFaninOld->pNtk );
    assert( pObj->pNtk == pFaninNewR->pNtk );
    assert( abc::Abc_ObjFanin(pObj, iFanin) == pFaninOld );

    // remember the attributes of the old fanin
    Vec_IntWriteEntry( &pObj->vFanins, iFanin, pFaninNewR->Id );
    if ( Abc_ObjIsComplement(pFaninNew) )
        Abc_ObjXorFaninC( pObj, iFanin );

    // update the fanout of the fanin
    if ( !Vec_IntRemoveRev( &pFaninOld->vFanouts, pObj->Id ) ) {
        printf( "Node %s is not among", Abc_ObjName(pObj) );
        printf( " the fanouts of its old fanin %s...\n", Abc_ObjName(pFaninOld) );
    }
    Vec_IntPushMem( pObj->pNtk->pMmStep, &pFaninNewR->vFanouts, pObj->Id );
}


void NetMan::Trunc(ll truncBit) {
    cout << "***** truncate " << truncBit << " bits" << endl;
    // truncation
    auto consts = CreateConst();
    #ifdef DEBUG
    assert(truncBit <= GetPoNum());
    #endif
    for (ll poId = 0; poId < truncBit; ++poId) {
        auto pPo = GetPo(poId);
        #ifdef DEBUG
        assert(GetFaninNum(pPo) == 1);
        #endif
        auto pDriv = GetFanin(pPo, 0);
        Abc_ObjPatchFanin(pPo, pDriv, GetObj(consts.first));
    }
    // clean up
    CleanUp();
}


void NetMan::SetBit(ll iBit, bool useConst1) {
    // cout << "***** set bit-" << iBit << " to " << useConst1 << endl;
    auto consts = CreateConst();
    assert(iBit <= GetPoNum());
    auto pPo = GetPo(iBit);
    assert(GetFaninNum(pPo) == 1);
    auto pDriv = GetFanin(pPo, 0);
    Abc_ObjPatchFanin(pPo, pDriv, useConst1? GetObj(consts.second): GetObj(consts.first));
    CleanUp();
}


bool NetMan::CleanUp() {
    bool isUpd = false;
    bool isCont = true;
    while (isCont) {
        isCont = false;
        // delete redundant node
        for (ll nodeId = 0; nodeId < GetIdMaxPlus1(); ++nodeId) {
            if (IsNode(nodeId) && GetFanoutNum(nodeId) == 0) {
                isCont = true;
                auto mffc = GetNodeMffc(GetObj(nodeId));
                for (const auto & pObj: mffc) {
                    // cout << "delete " << pObj << " ";
                    // if (GetNetType() == NET_TYPE::GATE) {
                    //     auto gateName = GetGateName(pObj);
                    //     cout << gateName;
                    //     if (gateName.find("HA1") != -1 || gateName.find("FA1") != -1) {
                    //         if (GetTwinNode(pObj) != nullptr)
                    //             cout << " twin " << GetTwinNode(pObj);
                    //     }
                    // }
                    // cout << endl;
                    DelObj(pObj);
                    isUpd = true;
                }
                break;
            }
        }
    }
    return isUpd;
}


// void NetMan::CleanUpPro() {
//     assert(GetNetType() == NET_TYPE::SOP);
//     // remove redundant nodes
//     Abc_NtkCleanup( pNtk, 1 );

//     // collect inverters
//     vector<Abc_Obj_t *> vInverters;
//     for (ll i = 0; i < GetIdMaxPlus1(); ++i) {
//         if (!IsNode(i)) 
//             continue;
//         // special case: f = !a & !a
//         assert(GetFaninNum(i) <= 2);
//         if (GetFaninNum(i) == 2) {
//             auto pFanin0 = GetFanin(i, 0);
//             auto pFanin1 = GetFanin(i, 1);
//             if (pFanin0 == pFanin1) {
//                 auto sop = (string)(GetSOP(i));
//                 if (sop == "11 0\n" || sop == "00 1\n")
//                 vInverters.emplace_back(GetObj(i));
//             }
//         }
//     }
//     // deal with each inverter
//     for (const auto pInv: vInverters) {
//         cout << "absorb inverter " << pInv << endl;
//         auto pFanin = GetFanin(pInv, 0);
//         ReplaceByComplementedObj(pInv->Id, pFanin->Id);
//     }

//     // remove buffers
//     // collect all buffers
//     vector<Abc_Obj_t *> vBuffers;
//     for (ll i = 0; i < GetIdMaxPlus1(); ++i) {
//         if (!IsNode(i)) 
//             continue;
//         if (IsBuffer(i))
//             vBuffers.emplace_back(GetObj(i));
//         // special case: f = a & a
//         assert(GetFaninNum(i) <= 2);
//         if (GetFaninNum(i) == 2) {
//             auto pFanin0 = GetFanin(i, 0);
//             auto pFanin1 = GetFanin(i, 1);
//             if (pFanin0 == pFanin1) {
//                 auto sop = (string)(GetSOP(i));
//                 if (sop == "11 1\n" || sop == "00 0\n")
//                 vBuffers.emplace_back(GetObj(i));
//             }
//         }
//     }
//     // for each buffer, replace it with its fanin
//     for (const auto pBuf: vBuffers) {
//         // cout << "remove buffer " << pBuf << endl;
//         auto pFanin = GetFanin(pBuf, 0);
//         Abc_ObjReplace(pBuf, pFanin);
//     }

// }


// bool NetMan::ProcHalfAndFullAdd() {
//     // special processing for half/full adder
//     if (GetNetType() != NET_TYPE::GATE)
//         return false;
//     bool isUpd = false;
//     ll idMaxPlus1 = GetIdMaxPlus1();
//     for (ll nodeId = 0; nodeId < idMaxPlus1; ++nodeId) {
//         if (!IsNode(nodeId)) continue;
//         auto pNode = GetObj(nodeId);
//         if (GetGateName(pNode).find("HA1") != -1) {
//             if (GetTwinNode(pNode) == nullptr) {
//                 cout << "cannot find twin for "; PrintObj(pNode, true);
//                 auto sop = string(Mio_GateReadSop((Mio_Gate_t *)pNode->pData));
//                 auto pLib = (Mio_Library_t *)Abc_FrameReadLibGen();
//                 auto pNewNode = Abc_NtkCreateNode(GetNet());
//                 for (ll faninId = 0; faninId < GetFaninNum(pNode); ++faninId)
//                     Abc_ObjAddFanin(pNewNode, GetFanin(pNode, faninId));
//                 Abc_Obj_t * pSub = nullptr;
//                 if (sop == "11 1\n") { // CO=A B
//                     auto pGate = Mio_LibraryReadGateByName(pLib, "CKAN2D1BWP7T30P140HVT", nullptr);
//                     assert(pGate != nullptr);
//                     pNewNode->pData = pGate;
//                     pSub = pNewNode;
//                 }
//                 else if (sop == "10 1\n01 1\n" || sop == "01 1\n10 1\n") { // S=A^B
//                     auto pGate = Mio_LibraryReadGateByName(pLib, "XOR2D0BWP7T30P140HVT", nullptr);
//                     assert(pGate != nullptr);
//                     pNewNode->pData = pGate;
//                     pSub = pNewNode;
//                 }
//                 else {
//                     cout << sop;
//                     assert(0);
//                 }
//                 assert(pSub != nullptr);
//                 cout << "replace " << pNode << " with new node " << pSub << endl;
//                 TransfFanout(pNode, pSub);
//                 isUpd = true;
//             }
//         }
//         else if (GetGateName(pNode).find("FA1") != -1) {
//             if (GetTwinNode(pNode) == nullptr) {
//                 cout << "cannot find twin for "; PrintObj(pNode, true);
//                 auto sop = string(Mio_GateReadSop((Mio_Gate_t *)pNode->pData));
//                 auto pLib = (Mio_Library_t *)Abc_FrameReadLibGen();
//                 auto pNewNode = Abc_NtkCreateNode(GetNet());
//                 for (ll faninId = 0; faninId < GetFaninNum(pNode); ++faninId)
//                     Abc_ObjAddFanin(pNewNode, GetFanin(pNode, faninId));
//                 Abc_Obj_t * pSub = nullptr;
//                 if (sop == "1-1 1\n-11 1\n11- 1\n") { // CO=A B+B CI+A CI
//                     auto pGate = Mio_LibraryReadGateByName(pLib, "MAOI222D0BWP7T30P140HVT", nullptr);
//                     assert(pGate != nullptr);
//                     pNewNode->pData = pGate;
//                     pSub = Abc_NtkCreateNodeInv(GetNet(), pNewNode);
//                 }
//                 else if (sop == "100 1\n010 1\n111 1\n001 1\n") { // S=A^B^CI
//                     auto pGate = Mio_LibraryReadGateByName(pLib, "XOR3D1BWP7T30P140HVT", nullptr);
//                     assert(pGate != nullptr);
//                     pNewNode->pData = pGate;
//                     pSub = pNewNode;
//                 }
//                 else {
//                     cout << sop;
//                     assert(0);
//                 }
//                 assert(pSub != nullptr);
//                 cout << "replace " << pNode << " with new node " << pSub << endl;
//                 TransfFanout(pNode, pSub);
//                 isUpd = true;
//             }
//         }
//     }
//     return isUpd;
// }


// void NetMan::ProcHalfAndFullAddNew() {
//     // special processing for half/full adder
//     if (GetNetType() != NET_TYPE::GATE)
//         return;
//     unordered_set <ll> vis;
//     vector <ll> targNodes;
//     for (ll iNode = 0; iNode < GetIdMaxPlus1(); ++iNode) {
//         if (!IsNode(iNode))
//             continue;
//         if (vis.count(iNode))
//             continue;
//         vis.emplace(iNode);
//         auto pNode = GetObj(iNode);
//         auto pTwin = GetTwinNode(pNode);
//         if (pTwin == nullptr)
//             continue;
//         vis.emplace(pTwin->Id);
//         if (GetGateName(pNode).find("HA1") != -1)
//             targNodes.emplace_back(iNode);
//         else if (GetGateName(pNode).find("FA1") != -1)
//             targNodes.emplace_back(iNode);
//         else
//             assert(0);
//     }

//     for (ll targId: targNodes) {
//         auto pNode = GetObj(targId);
//         auto pTwin = GetTwinNode(pNode);
//         assert(pTwin != nullptr);
//         if (GetGateName(pNode).find("HA1") != -1) {
//             // print
//             // PrintObj(pNode, true); 
//             auto sop0 = string(Mio_GateReadSop((Mio_Gate_t *)pNode->pData));
//             // cout << sop0;
//             // PrintObj(pTwin, true);
//             auto sop1 = string(Mio_GateReadSop((Mio_Gate_t *)pTwin->pData));
//             // cout << sop1;
//             // cout << endl;

//             // pNode sop0 S, pTwin sop1 CO
//             assert(sop0 == "10 1\n01 1\n" && sop1 == "11 1\n");
//             vector <Abc_Obj_t *> fanins;
//             ll nFanin = GetFaninNum(pNode);
//             assert(nFanin == GetFaninNum(pTwin) && nFanin == 2);
//             for (ll iFanin = 0; iFanin < nFanin; ++iFanin) {
//                 assert(GetFanin(pNode, iFanin) == GetFanin(pTwin, iFanin));
//                 fanins.emplace_back(GetFanin(pNode, iFanin));
//             }

//             // create gates
//             auto pNodeCo = CreateGate(std::vector <Abc_Obj_t *> ({fanins[0], fanins[1]}), "CKAN2D1BWP7T30P140HVT");
//             auto pNodeN6 = CreateGate(std::vector <Abc_Obj_t *> ({fanins[0], fanins[1]}), "NR2D0BWP7T30P140HVT");
//             auto pNodeS = CreateGate(std::vector <Abc_Obj_t *> ({pNodeCo, pNodeN6}), "NR2D0BWP7T30P140HVT");

//             // cout << "replace " << pTwin << " with new node " << pNodeCo << endl;
//             TransfFanout(pTwin, pNodeCo);
//             // cout << "replace " << pNode << " with new node " << pNodeS << endl;
//             TransfFanout(pNode, pNodeS);
//         }
//         else if (GetGateName(pNode).find("FA1") != -1) {
//             // print
//             // PrintObj(pNode, true); 
//             auto sop0 = string(Mio_GateReadSop((Mio_Gate_t *)pNode->pData));
//             // cout << sop0;
//             // PrintObj(pTwin, true);
//             auto sop1 = string(Mio_GateReadSop((Mio_Gate_t *)pTwin->pData));
//             // cout << sop1;
//             // cout << endl;

//             // pNode sop0 S, pTwin sop1 CO
//             assert(sop0 == "100 1\n010 1\n111 1\n001 1\n" && sop1 == "1-1 1\n-11 1\n11- 1\n");
//             vector <Abc_Obj_t *> fanins;
//             ll nFanin = GetFaninNum(pNode);
//             assert(nFanin == GetFaninNum(pTwin) && nFanin == 3);
//             for (ll iFanin = 0; iFanin < nFanin; ++iFanin) {
//                 assert(GetFanin(pNode, iFanin) == GetFanin(pTwin, iFanin));
//                 fanins.emplace_back(GetFanin(pNode, iFanin));
//             }

//             // create gates
//             auto pNodeN6 = CreateGate(std::vector <Abc_Obj_t *> ({fanins[1], fanins[2], fanins[1], fanins[2]}), "MOAI22D0BWP7T30P140HVT");
//             auto pNodeS = CreateGate(std::vector <Abc_Obj_t *> ({fanins[0], pNodeN6, fanins[0], pNodeN6}), "MOAI22D0BWP7T30P140HVT");
//             auto pNodeCo = CreateGate(std::vector <Abc_Obj_t *> ({fanins[1], fanins[2], fanins[0], pNodeN6}), "OA22D0BWP7T30P140HVT");

//             // cout << "replace " << pTwin << " with new node " << pNodeCo << endl;
//             TransfFanout(pTwin, pNodeCo);
//             // cout << "replace " << pNode << " with new node " << pNodeS << endl;
//             TransfFanout(pNode, pNodeS);
//         }
//     }
//     CleanUp();
// }

abc::Abc_Obj_t * NetMan::CreateGate(vector <Abc_Obj_t *> && fanins, const std::string & gateName) {
    auto pLib = (Mio_Library_t *)Abc_FrameReadLibGen();
    auto pGate = Mio_LibraryReadGateByName(pLib, const_cast <char *> (gateName.c_str()), nullptr);
    assert(pGate != nullptr);
    auto pNewNode = Abc_NtkCreateNode(GetNet());
    for (const auto & fanin: fanins)
        Abc_ObjAddFanin(pNewNode, fanin);
    pNewNode->pData = pGate;
    return pNewNode;
}


void NetMan::ReplaceByComplementedObj(ll targId, ll subId) {
    assert(GetNetType() == NET_TYPE::SOP);
    auto pTarg = GetObj(targId);
    auto pSub = GetObj(subId);
    assert(pTarg != nullptr && pSub != nullptr);
    // collect fanouts of pTarg
    vector <Abc_Obj_t *> fanouts;
    for (int i = 0; i < GetFanoutNum(pTarg); ++i)
        fanouts.emplace_back(GetFanout(pTarg, i));
    // for each pTarg's fanout (marked as fo), collect the fanin index of pTarg; 
    // meanwhile, make sure that pTarg only appears once in fo's fanin list
    vector <ll> foIFanin;
    for (const auto & fo: fanouts) {
        vector <ll> iFanins;
        for (int i = 0; i < GetFaninNum(fo); ++i) {
            if (GetFanin(fo, i) == pTarg)
                iFanins.emplace_back(i);
        }
        assert(iFanins.size() == 1);
        foIFanin.emplace_back(iFanins[0]);
    }
    // transfer fanouts
    Abc_ObjTransferFanout( pTarg, pSub );
    // fix SOP functions for fanouts
    Abc_Obj_t * pInv = nullptr;
    for (ll i = 0; i < fanouts.size(); ++i) {
        auto fo = fanouts[i];
        auto iFanin = foIFanin[i];
        if (IsObjPo(fo)) {
            if (pInv == nullptr) {
                pInv = Abc_NtkCreateNodeInv( GetNet(), pSub );
                cout << "create inverter for " << pSub << ":" << pInv << endl;
            } 
            // because pTarg's fanouts have been transferred to pSub, so pSub is the only fanin of fo
            assert(GetFaninNum(fo) == 1 && GetFanin(fo, 0) == pSub);
            Abc_ObjPatchFanin( fo, pSub, pInv );
        }
        else {
            auto pSop = static_cast <char *> (fo->pData);
            Abc_SopComplementVar( pSop, iFanin );
        }
    }
    // remove redundant nodes
    Abc_NtkDeleteObj_rec( pTarg, 1 );
}


void NetMan::SetConstantInput(Abc_Obj_t * pNode, Abc_Obj_t * pFanin, int fConst0) {
    int iFanin = Vec_IntFind( &pNode->vFanins, pFanin->Id );
    if (iFanin == -1 ) {
        printf( "Node %s should be among", Abc_ObjName(pFanin) );
        printf( " the fanins of node %s...\n", Abc_ObjName(pNode) );
        assert(0);
        return;
    }

    // construct new sop
    string newSop = "";
    auto pOldSop = static_cast<char*>(pNode->pData);
    bool isOldSopComplemented = Abc_SopIsComplement(pOldSop);
    int nVars = Abc_SopGetVarNum(pOldSop);
    assert(iFanin < nVars);
    char * pCube = nullptr;
    Abc_SopForEachCube(pOldSop, nVars, pCube) {
        // if ((fConst0 && pCube[iFanin] == '0') || (!fConst0 && pCube[iFanin] == '1')) 
        if ((fConst0 && pCube[iFanin] != '1') || (!fConst0 && pCube[iFanin] != '0')) // consider don't care
        {
            // append pCube to newSop except for the iFanin-th bit
            for (int i = 0; i < nVars; ++i) {
                if (i == iFanin)
                    continue;
                newSop += pCube[i];
            }
            newSop += isOldSopComplemented? " 0\n": " 1\n";
        }
    }
    if (newSop == "")
        newSop = isOldSopComplemented? " 1\n": " 0\n";
    // cout << "old fanins: ";
    // for (int i = 0; i < GetFaninNum(pNode); ++i)
    //     cout << GetFanin(pNode, i) << " ";
    // cout << endl;
    // cout << "old sop: " << pOldSop;

    // remove pNode's fanin: pFanin
    if (newSop == " 1\n" || newSop == " 0\n")
        Abc_ObjRemoveFanins(pNode);
    else
        Abc_ObjDeleteFanin( pNode, pFanin );
    // cout << "pFanin: " << pFanin << endl;
    // cout << "pFanin's fanouts: ";
    // for (int i = 0; i < GetFanoutNum(pFanin); ++i)
    //     cout << GetFanout(pFanin, i) << " ";
    // cout << endl;

    // update pNode's sop
    pNode->pData = Abc_SopRegister((Mem_Flex_t *)GetNet()->pManFunc, newSop.c_str());
    // cout << "new fanins: ";
    // for (int i = 0; i < GetFaninNum(pNode); ++i)
    //     cout << GetFanin(pNode, i) << " ";
    // cout << endl;
    // cout << "new sop: " << newSop;
    // cout << endl;
}


void NetMan::PropagateConst(ll startId) {
    assert(GetNetType() == NET_TYPE::SOP);
    // initialize
    auto vNodes = Vec_PtrAlloc( 100 );
    Vec_PtrPush( vNodes, GetObj(startId) );
    // sweep the nodes
    while ( Vec_PtrSize(vNodes) > 0 ) {
        // get any sweepable node
        auto pNode = (Abc_Obj_t *)Vec_PtrPop(vNodes);
        if ( !Abc_ObjIsNode(pNode) )
            continue;
        // get any non-CO fanout of this node
        auto pFanout = Abc_NodeFindNonCoFanout(pNode);
        if ( pFanout == nullptr ) 
            continue;
        cout << "propagate " << pNode << " to its fanout " << pFanout << endl;
        assert( Abc_ObjIsNode(pFanout) );
        // transform the function of the fanout
        if ( Abc_ObjFaninNum(pNode) == 0 ) {
            // cout << "set const input " << pFanout << " " << pNode << " " << Abc_NodeIsConst0(pNode) << endl;
            SetConstantInput( pFanout, pNode, Abc_NodeIsConst0(pNode) );
        }
        else 
        {
            assert( Abc_ObjFaninNum(pNode) == 1 );
            auto pDriver = Abc_ObjFanin0(pNode);
            if ( Abc_NodeIsInv(pNode) )
                Abc_NodeComplementInput( pFanout, pNode );
            Abc_ObjPatchFanin( pFanout, pNode, pDriver );
        }
        // check if the fanout should be added
        if ( Abc_ObjFaninNum(pFanout) < 2 )
            Vec_PtrPush( vNodes, pFanout );
        // check if the node has other fanouts
        if ( Abc_ObjFanoutNum(pNode) > 0 )
            Vec_PtrPush( vNodes, pNode );
        else
            Abc_NtkDeleteObj_rec( pNode, 1 );
    }
    Vec_PtrFree( vNodes );
}


bool NetMan::IsNodeXor(ll nodeId, ll& div0, ll& div1, bool print) {
    div0 = -1; 
    div1 = -1;
    if (GetFaninNum(nodeId) != 2)
        return false;

    ll fanin0 = GetFaninId(nodeId, 0), fanin1 = GetFaninId(nodeId, 1);
    if (fanin0 > fanin1) 
        swap(fanin0, fanin1);

    if (GetFaninNum(fanin0) != 2 || GetFaninNum(fanin1) != 2) 
        return false;

    ll fanin0Fanin0 = GetFaninId(fanin0, 0), fanin0Fanin1 = GetFaninId(fanin0, 1);
    if (fanin0Fanin0 > fanin0Fanin1) 
        swap(fanin0Fanin0, fanin0Fanin1);
    ll fanin1Fanin0 = GetFaninId(fanin1, 0), fanin1Fanin1 = GetFaninId(fanin1, 1);
    if (fanin1Fanin0 > fanin1Fanin1) 
        swap(fanin1Fanin0, fanin1Fanin1);
    if (fanin0Fanin0 == fanin1Fanin0 && fanin0Fanin1 == fanin1Fanin1) {
        div0 = fanin0Fanin0;
        div1 = fanin0Fanin1;
        if (GetSOP(nodeId) == string("00 0\n") || GetSOP(nodeId) == string("00 1\n")) { // OR/NOR
            unordered_set<string> sopSet;
            sopSet.emplace(GetSOP(fanin0));
            sopSet.emplace(GetSOP(fanin1));
            if (sopSet.count(string("01 1\n")) && sopSet.count(string("10 1\n")))
                return true;
            else if (sopSet.count(string("00 1\n")) && sopSet.count(string("11 1\n")))
                return true;
        }
    }

    return false;
}


std::string FixName(const char* name) {
    std::string newName(name);
    for (auto& ch: newName) {
        if (ch == '[' || ch == ']')
            ch = '_';
    }
    return newName;
}

void NetMan::DumpCFile(std::string&& fileName) {
    #ifdef DEBUG
    assert(isDupl == true);
    assert(this->GetNetType() == NET_TYPE::SOP || this->GetNetType() == NET_TYPE::GATE);
    #endif

    std::cout << "write " << fileName << endl;
    FILE * file = fopen(fileName.c_str(), "w");
    abc::Abc_Obj_t * pObj = nullptr;
    int i = 0;

    fprintf(file, "#include <stdbool.h>\n");
    fprintf(file, "void %s(", pNtk->pName);
    Abc_NtkForEachPi(pNtk, pObj, i)
        fprintf(file, "bool %s, ", FixName(Abc_ObjName(pObj)).c_str());
    Abc_NtkForEachPo(pNtk, pObj, i) {
        if (i != Abc_NtkPoNum(pNtk) - 1)
            fprintf(file, "bool& %s, ", FixName(Abc_ObjName(pObj)).c_str());
        else
            fprintf(file, "bool& %s)\n", FixName(Abc_ObjName(pObj)).c_str());
    }
    fprintf(file, "{\n");
    Abc_NtkForEachNode(pNtk, pObj, i) {
        std::ostringstream oss("");
        oss << "bool " << FixName(Abc_ObjName(pObj)).c_str() << " = ";

        if (Abc_NodeIsConst0(pObj)) {
            oss << "0;\n";
        }
        else if (Abc_NodeIsConst1(pObj)) {
            oss << "1;\n";
        }
        else {
            char *pSop = nullptr;
            if (this->GetNetType() == NET_TYPE::SOP)
                pSop = static_cast<char*>(pObj->pData);
            else if (this->GetNetType() == NET_TYPE::GATE)
                pSop = Mio_GateReadSop(static_cast<Mio_Gate_t*>(pObj->pData));
            else
                assert(0);
            int nVars = abc::Abc_SopGetVarNum(pSop);
            if (abc::Abc_SopIsComplement(pSop))
                oss << "!(\n";
            else
                oss << "(\n";
            for (char * pCube = pSop; *pCube; pCube += nVars + 3) {
                if (pCube == pSop)
                    oss << "( ";
                else
                    oss << "|| ( ";
                bool isFirst = true;
                for (int k = 0; pCube[k] != ' '; k++) {
                    abc::Abc_Obj_t * pFanin = Abc_ObjFanin(pObj, k);
                    std::string faninName = std::string(FixName(Abc_ObjName(pFanin)).c_str());
                    if (isFirst) {
                        if (pCube[k] == '0') {
                            isFirst = false;
                            oss << "!" << faninName << " ";
                        }
                        else if (pCube[k] == '1') {
                            isFirst = false;
                            oss << faninName << " ";
                        }
                        else if (pCube[k] == '-')
                            ;
                        else
                            assert(0);
                    }
                    else {
                        if (pCube[k] == '0')
                            oss << "&& !" << faninName << " ";
                        else if (pCube[k] == '1')
                            oss << "&& " << faninName << " ";
                        else if (pCube[k] == '-')
                            ;
                        else
                            assert(0);
                    }
                }
                oss << ")\n";
            }
            oss << ");\n";
        }
        fprintf(file, "%s", oss.str().c_str());
    }
    Abc_NtkForEachPo(pNtk, pObj, i) {
        abc::Abc_Obj_t * pDriver = Abc_ObjFanin0(pObj);
        assert(Abc_ObjIsNode(pDriver));
        fprintf(file, "%s = ", FixName(Abc_ObjName(pObj)).c_str());
        fprintf(file, "%s;\n", FixName(Abc_ObjName(pDriver)).c_str());
    }
    fprintf(file, "}\n");
    fclose(file);
}


void GlobStartAbc() {
    Abc_Start();
    AbcMan abcMan;
    abcMan.LoadAlias();
}


void GlobStopAbc() {
    Abc_Stop();
}


void EvalNetw(NetMan & net, const std::string & outpPath, double err,  ll round) {
    if (!IsPathExist(outpPath))
        CreatePath(outpPath);
    if (net.GetNetType() == NET_TYPE::SOP) {
        // measure and output SOP
        ostringstream oss("");
        oss << outpPath << round << "_" << net.GetNet()->pName << "_error_" << err << "_size_" << net.GetArea() << "_depth_" << net.GetDelay();
        net.WriteNet(oss.str() + ".blif");
    }
    else if (net.GetNetType() == NET_TYPE::GATE) {
        // measure and output original gate-netlist
        auto tempNet = net;
        tempNet.ReArrInTopoOrd();
        ostringstream oss("");
        oss << outpPath << round << "_" << tempNet.GetNet()->pName << "_error_" << err << "_area_" << tempNet.GetArea() << "_delay_" << tempNet.GetDelay();
        tempNet.WriteNet(oss.str() + ".v", true);
        tempNet.ConvToSop();
        tempNet.WriteNet(oss.str() + ".blif", true);
    }
    else
        assert(0);
}