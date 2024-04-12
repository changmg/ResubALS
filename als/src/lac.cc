#include "lac.h"


using namespace std;
using namespace abc;


// bool RacLAC::Check() {
//     auto pSop = const_cast<char*>(sop.c_str());
//     char * pCube = nullptr;
//     int nVars = divs.size();
//     int value = 0;
//     int i = 0;
//     assert(nVars == Abc_SopGetVarNum(pSop));
//     boost::dynamic_bitset<ull> varRedund(nVars);
//     varRedund.set();
//     Abc_SopForEachCube(pSop, nVars, pCube) {
//         Abc_CubeForEachVar(pCube, value, i)
//             varRedund[i] &= (value == '-');
//     }
//     for (i = 0; i < nVars; ++i) {
//         if (varRedund[i]) {
//             // cout << varRedund << endl;
//             // Print();
//             // assert(0);
//             return false;
//         }
//     }
//     return true;
// }


bool ResubLAC::Check() {
    auto pSop = const_cast<char*>(sop.c_str());
    char * pCube = nullptr;
    int nVars = divs.size();
    int value = 0;
    int i = 0;
    assert(nVars == Abc_SopGetVarNum(pSop));
    boost::dynamic_bitset<ull> varRedund(nVars);
    varRedund.set();
    Abc_SopForEachCube(pSop, nVars, pCube) {
        Abc_CubeForEachVar(pCube, value, i)
            varRedund[i] &= (value == '-');
    }
    for (i = 0; i < nVars; ++i) {
        if (varRedund[i]) {
            // cout << varRedund << endl;
            // Print();
            // assert(0);
            return false;
        }
    }
    return true;
}


ll ResubLAC::GetAddedNodeNum() {
    auto pSop = const_cast<char*>(sop.c_str());
    char * pCube = nullptr;
    int nVars = divs.size();
    int value = 0;
    int i = 0;
    assert(nVars == Abc_SopGetVarNum(pSop));
    assert(Abc_SopGetCubeNum(pSop) >= 1);
    ll addedNodeNum = Abc_SopGetCubeNum(pSop) - 1; // number of ORs
    Abc_SopForEachCube(pSop, nVars, pCube) { // number of ANDs
        ll nAnd = 0;
        Abc_CubeForEachVar(pCube, value, i) {
            if (value != '-')
                ++nAnd;
        }
        addedNodeNum += (nAnd? nAnd - 1: 0);
    }
    return addedNodeNum;
}


void LACMan::Gen012ResubLACsPro(NetMan& net, IntVect& nodeIds, unsigned seed, int maxLevelDiff, int nFrame4ResubGen, int maxCandResub) {
    // initialize
    int simulationFrame = nFrame4ResubGen;
    int halfFrame = simulationFrame >> 1;
    int LAC_NUM_LIMIT = maxCandResub;
    assert(net.GetNetType() == NET_TYPE::SOP);
    pLacs.clear();
    net.GetLev();
    Abc_NtkStartReverseLevels(net.GetNet(), 0);

    // simulate
    Simulator smlt(net, seed, simulationFrame);
    if (simulationFrame < 64)
        smlt.InpUnif();
    else
        smlt.InpUnifFast();
    smlt.Sim();

    // collect target nodes
    IntVect targIds;
    for (int nodeId: nodeIds) {
        if (net.IsNode(nodeId) && !net.IsConst(nodeId) && net.GetFaninNum(nodeId) > 1)
            targIds.emplace_back(nodeId);
    }

    // collect divisors
    cout << "collecting divisors" << endl;
    vector<IntVect> divs4Nodes;
    divs4Nodes.resize(net.GetIdMaxPlus1());
    boost::timer::progress_display pd(targIds.size());
    for (int targId: targIds) {
        auto pTarget = net.GetObj(targId);
        GetDivs(pTarget, abc::Abc_ObjRequiredLevel(pTarget) - 1, divs4Nodes[targId]);
        ++pd;
    }
    // generate 0 resubstitution
    for (int targId: targIds) {
        int sizeGain = net.IsTheOnlyPoDriver(targId)? net.GetSizeGain(targId, IntVect{}): net.GetSizeGain(targId, IntVect{}) + 1;
        if (smlt.GetDat(targId)->count() <= halfFrame)
            pLacs.emplace_back(make_shared<ResubLAC>(targId, sizeGain, IntVect{}, string(" 0\n")));
        else
            pLacs.emplace_back(make_shared<ResubLAC>(targId, sizeGain, IntVect{}, string(" 1\n")));
    }
    cout << "generated 0-resubs, #lacs = " << pLacs.size() << endl;
    // generate 1 resubstitution
    for (int targId: targIds) {
        auto& divs = divs4Nodes[targId];
        for (int div: divs) {
            auto diff = (*smlt.GetDat(div) ^ *smlt.GetDat(targId)).count();
            if (diff == 0) {
                int sizeGain = net.GetSizeGain(targId, IntVect{div});
                pLacs.emplace_back(make_shared<ResubLAC>(targId, sizeGain, IntVect{div}, string("1 1\n")));
                if (pLacs.size() > LAC_NUM_LIMIT)
                    break;
            }
            else if (diff == simulationFrame) {
                int sizeGain = net.GetSizeGain(targId, IntVect{div});
                pLacs.emplace_back(make_shared<ResubLAC>(targId, sizeGain, IntVect{div}, string("0 1\n")));
                if (pLacs.size() > LAC_NUM_LIMIT)
                    break;
            }
        }
    }
    cout << "generated 1-resubs, #lacs = " << pLacs.size() << endl;
    // generate 2 resubstitution: try replacing the i-th fanin with another divisor
    if (pLacs.size() <= LAC_NUM_LIMIT) {
    bool _break = false;
    for (int targId: targIds) {
        if (_break)
            break;
        int nFanin = net.GetFaninNum(targId);
        assert(nFanin == 2);
        int fanin0 = net.GetFaninId(targId, 0), fanin1 = net.GetFaninId(targId, 1);
        auto& divs = divs4Nodes[targId];
        for (int i = 0; i < nFanin; ++i) {
            int remainedFanin = (i == 0)? fanin1: fanin0;
            int replacedFanin = net.GetFaninId(targId, i);
            auto faninIds = IntVect{remainedFanin, -1};
            for (int div: divs) {
                if (div == replacedFanin || div == remainedFanin)
                    continue;
                if (_break)
                    break;
                faninIds[1] = div;
                int sizeGain = net.GetSizeGain(targId, faninIds) - 1;
                if (sizeGain >= 1) {
                    for (int comb = 0; comb < 4; ++comb) {
                        int var0 = (comb >> 1) & 1, var1 = comb & 1;
                        auto dat0 = var0? *smlt.GetDat(faninIds[0]): ~(*smlt.GetDat(faninIds[0]));
                        auto dat1 = var1? *smlt.GetDat(faninIds[1]): ~(*smlt.GetDat(faninIds[1]));
                        auto res = dat0 & dat1;
                        auto diff = (res ^ *smlt.GetDat(targId)).count();
                        if (diff == 0) {
                            pLacs.emplace_back(make_shared<ResubLAC>(targId, sizeGain, faninIds, to_string(var0) + to_string(var1) + string(" 1\n")));
                            _break = (pLacs.size() > LAC_NUM_LIMIT);
                        }
                        else if (diff == simulationFrame) {
                            pLacs.emplace_back(make_shared<ResubLAC>(targId, sizeGain, faninIds, to_string(var0) + to_string(var1) + string(" 0\n")));
                            _break = (pLacs.size() > LAC_NUM_LIMIT);
                        }
                    }
                }
            }
        }
    }
    cout << "generated 2-resubs, #lacs = " << pLacs.size() << endl;
    }
    // generate 2 resubstitution: try replacing both fanins with new divisors, using AND-based functions
    if (pLacs.size() <= LAC_NUM_LIMIT) {
    bool _break = false;
    for (int targId: targIds) {
        if (_break)
            break;
        assert(net.GetFaninNum(targId) == 2);
        int fanin0 = net.GetFaninId(targId, 0), fanin1 = net.GetFaninId(targId, 1);
        auto& divs = divs4Nodes[targId];
        IntVect posDiv4And, negDiv4And, posDiv4Or, negDiv4Or;
        posDiv4And.reserve(divs.size());
        negDiv4And.reserve(divs.size());
        posDiv4Or.reserve(divs.size());
        negDiv4Or.reserve(divs.size());
        for (int div: divs) {
            if (div == fanin0 || div == fanin1)
                continue;
            int levelDiff = abc::Abc_ObjReverseLevel(net.GetObj(div)) - abc::Abc_ObjReverseLevel(net.GetObj(targId));
            if (levelDiff > maxLevelDiff)
                continue;
            auto check = (*smlt.GetDat(targId) & *smlt.GetDat(div)) | (~*smlt.GetDat(targId)); // targ = div & xxx, targ = 1 \Rightarrow div = 1
            if (check.all())
                posDiv4And.emplace_back(div);
            check = (*smlt.GetDat(targId) & ~*smlt.GetDat(div)) | (~*smlt.GetDat(targId)); // targ = ~div & xxx, targ = 1 \Rightarrow div = 0
            if (check.all())
                negDiv4And.emplace_back(div);
            check = (~*smlt.GetDat(targId) & ~*smlt.GetDat(div)) | (*smlt.GetDat(targId)); // targ = div | xxx, targ = 0 \Rightarrow div = 0
            if (check.all())
                posDiv4Or.emplace_back(div);
            check = (~*smlt.GetDat(targId) & *smlt.GetDat(div)) | (*smlt.GetDat(targId)); // targ = ~div | xxx, targ = 0 \Rightarrow div = 1
            if (check.all())
                negDiv4Or.emplace_back(div);
        }
        // try targ = div0 & xxx
        for (int div0: posDiv4And) {
            // try targ = div0 & div1
            for (int div1: posDiv4And) {
                int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
                if (sizeGain >= 1) {
                    auto diff = (*smlt.GetDat(div0) & *smlt.GetDat(div1)) ^ *smlt.GetDat(targId);
                    if (diff.none()) {
                        pLacs.emplace_back(make_shared<ResubLAC>(targId, sizeGain, IntVect{div0, div1}, string("11 1\n")));
                        _break = (pLacs.size() > LAC_NUM_LIMIT);
                    }
                }
            }
            // try targ = div0 & ~div1
            for (int div1: negDiv4And) {
                int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
                if (sizeGain >= 1) {
                    auto diff = (*smlt.GetDat(div0) & ~*smlt.GetDat(div1)) ^ *smlt.GetDat(targId);
                    if (diff.none()) {
                        pLacs.emplace_back(make_shared<ResubLAC>(targId, sizeGain, IntVect{div0, div1}, string("10 1\n")));
                        _break = (pLacs.size() > LAC_NUM_LIMIT);
                    }
                }
            }
        }
        // try targ = ~div0 & xxx
        for (int div0: negDiv4And) {
            // try targ = ~div0 & div1
            for (int div1: posDiv4And) {
                int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
                if (sizeGain >= 1) {
                    auto diff = (~*smlt.GetDat(div0) & *smlt.GetDat(div1)) ^ *smlt.GetDat(targId);
                    if (diff.none()) {
                        pLacs.emplace_back(make_shared<ResubLAC>(targId, sizeGain, IntVect{div0, div1}, string("01 1\n")));
                        _break = (pLacs.size() > LAC_NUM_LIMIT);
                    }
                }
            }
            // try targ = ~div0 & ~div1
            for (int div1: negDiv4And) {
                int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
                if (sizeGain >= 1) {
                    auto diff = (~*smlt.GetDat(div0) & ~*smlt.GetDat(div1)) ^ *smlt.GetDat(targId);
                    if (diff.none()) {
                        pLacs.emplace_back(make_shared<ResubLAC>(targId, sizeGain, IntVect{div0, div1}, string("00 1\n")));
                        _break = (pLacs.size() > LAC_NUM_LIMIT);
                    }
                }
            }
        }
        // try targ = div0 | xxx
        for (int div0: posDiv4Or) {
            // try targ = div0 | div1
            for (int div1: posDiv4Or) {
                int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
                if (sizeGain >= 1) {
                    auto diff = (*smlt.GetDat(div0) | *smlt.GetDat(div1)) ^ *smlt.GetDat(targId);
                    if (diff.none()) {
                        pLacs.emplace_back(make_shared<ResubLAC>(targId, sizeGain, IntVect{div0, div1}, string("00 0\n")));
                        _break = (pLacs.size() > LAC_NUM_LIMIT);
                    }
                }
            }
            // try targ = div0 | ~div1
            for (int div1: negDiv4Or) {
                int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
                if (sizeGain >= 1) {
                    auto diff = (*smlt.GetDat(div0) | ~*smlt.GetDat(div1)) ^ *smlt.GetDat(targId);
                    if (diff.none()) {
                        pLacs.emplace_back(make_shared<ResubLAC>(targId, sizeGain, IntVect{div0, div1}, string("01 0\n")));
                        _break = (pLacs.size() > LAC_NUM_LIMIT);
                    }
                }
            }
        }
        // try targ = ~div0 | xxx
        for (int div0: negDiv4Or) {
            // try targ = ~div0 | div1
            for (int div1: posDiv4Or) {
                int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
                if (sizeGain >= 1) {
                    auto diff = (~*smlt.GetDat(div0) | *smlt.GetDat(div1)) ^ *smlt.GetDat(targId);
                    if (diff.none()) {
                        pLacs.emplace_back(make_shared<ResubLAC>(targId, sizeGain, IntVect{div0, div1}, string("10 0\n")));
                        _break = (pLacs.size() > LAC_NUM_LIMIT);
                    }
                }
            }
            // try targ = ~div0 | ~div1
            for (int div1: negDiv4Or) {
                int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
                if (sizeGain >= 1) {
                    auto diff = (~*smlt.GetDat(div0) | ~*smlt.GetDat(div1)) ^ *smlt.GetDat(targId);
                    if (diff.none()) {
                        pLacs.emplace_back(make_shared<ResubLAC>(targId, sizeGain, IntVect{div0, div1}, string("11 0\n")));
                        _break = (pLacs.size() > LAC_NUM_LIMIT);
                    }
                }
            }
        }
    }
    cout << "generated more 2-resubs, #lacs = " << pLacs.size() << endl;
    }

    // clean up
    abc::Abc_NtkStopReverseLevels(net.GetNet());
}


int LACMan::CollectPromisingLACs(Simulator& accSmlt, Simulator& appSmlt, ll errorMargin) {
    // check
    assert(IsPIOSame(accSmlt, appSmlt));
    assert(lacType == LAC_TYPE::RESUB);
    assert(pLacs.size());

    // initialize
    pPromisingLacs.clear();
    pPromisingLacs.resize(appSmlt.GetIdMaxPlus1());
    pPromisingLacPatterns.clear();
    pPromisingLacPatterns.resize(appSmlt.GetIdMaxPlus1());
    int oldTargId = -1;
    unordered_map<int, int> size2BestId;
    int count = 0;

    // assume that the LACs are sorted by targId
    for (auto i = 0u; i < pLacs.size(); ++i) {
        auto pLac = dynamic_pointer_cast<ResubLAC>(pLacs[i]);
        int targId = pLac->GetTargId();
        int sizeGain = pLac->GetSizeGain();
        auto deltaError = pLac->GetErrPro();
        if (oldTargId != -1 && oldTargId != targId) { // need to store promising LACs for the old node and deal with a new node
            pPromisingLacs[oldTargId].reserve(size2BestId.size());
            count += size2BestId.size();
            // if (size2BestId.size())
            //     cout << "node " << oldTargId << ": " << endl;
            for (const auto& [sizeGain, bestLacId]: size2BestId) {
                pPromisingLacs[oldTargId].emplace_back(pLacs[bestLacId]);
            }
            size2BestId.clear();
        }
        if (deltaError <= errorMargin) {
            if (size2BestId.count(sizeGain) == 0)
                size2BestId[sizeGain] = i;
            else {
                auto _id = size2BestId[sizeGain];
                if (pLacs[_id]->GetErrPro() > deltaError)
                    size2BestId[sizeGain] = i;
            }
        }
        oldTargId = targId;
    }

    // store the last node
    pPromisingLacs[oldTargId].reserve(size2BestId.size());
    pPromisingLacPatterns[oldTargId].reserve(size2BestId.size());
    count += size2BestId.size();
    // if (size2BestId.size())
    //     cout << "node " << oldTargId << ": " << endl;
    for (const auto& [sizeGain, bestLacId]: size2BestId) {
        pPromisingLacs[oldTargId].emplace_back(pLacs[bestLacId]);
    }
    cout << "#promising LACs: " << count << endl;
    return count;
}


int LACMan::CollectPromisingLACs(ll nObjs, BigInt& errorMargin) {
    // check
    assert(lacType == LAC_TYPE::RESUB);
    assert(pLacs.size());

    // initialize
    pPromisingLacs.clear();
    pPromisingLacs.resize(nObjs);
    int oldTargId = -1;
    unordered_map<int, int> size2BestId;
    int count = 0;

    // assume that the LACs are sorted by targId
    for (auto i = 0u; i < pLacs.size(); ++i) {
        auto pLac = dynamic_pointer_cast<ResubLAC>(pLacs[i]);
        int targId = pLac->GetTargId();
        int sizeGain = pLac->GetSizeGain();
        auto deltaError = pLac->GetErrPro();
        if (oldTargId != -1 && oldTargId != targId) { // need to store promising LACs for the old node and deal with a new node
            pPromisingLacs[oldTargId].reserve(size2BestId.size());
            count += size2BestId.size();
            // if (size2BestId.size())
            //     cout << "node " << oldTargId << ": " << endl;
            for (const auto& [sizeGain, bestLacId]: size2BestId) {
                pPromisingLacs[oldTargId].emplace_back(pLacs[bestLacId]);
                // dynamic_pointer_cast<ResubLAC>(pLacs[bestLacId])->Print();
            }
            size2BestId.clear();
        }
        if (deltaError <= errorMargin) {
            if (size2BestId.count(sizeGain) == 0)
                size2BestId[sizeGain] = i;
            else {
                auto _id = size2BestId[sizeGain];
                if (pLacs[_id]->GetErrPro() > deltaError)
                    size2BestId[sizeGain] = i;
            }
        }
        oldTargId = targId;
    }

    // store the last node
    pPromisingLacs[oldTargId].reserve(size2BestId.size());
    count += size2BestId.size();
    // if (size2BestId.size())
    //     cout << "node " << oldTargId << ": " << endl;
    for (const auto& [sizeGain, bestLacId]: size2BestId)
        pPromisingLacs[oldTargId].emplace_back(pLacs[bestLacId]);
    cout << "#promising LACs: " << count << endl;
    return count;
}


shared_ptr <LAC> LACMan::GetBestLac() const {
    assert(pLacs.size());
    auto pBestLac = pLacs[0];
    BigInt bestErr = pBestLac->GetErrPro();
    for (ll i = 1; i < pLacs.size(); ++i) {
        auto pLac = pLacs[i];
        BigInt err = pLac->GetErrPro();
        if (bestErr >= err) {
            bestErr = err;
            pBestLac = pLac;
        }
    }
    return pBestLac;
}


shared_ptr <LAC> LACMan::GetBestResubLacConsiderRealSize(NetMan& net) const {
    assert(pLacs.size());
    shared_ptr<LAC> pBestLac = nullptr;
    BigInt bestErr = numeric_limits <BigInt>::max();
    ll bestSizeGain = -1;
    for (ll i = 0; i < pLacs.size(); ++i) {
        auto pLac = pLacs[i];
        auto pResubLac = dynamic_pointer_cast <ResubLAC> (pLac);

        BigInt err = pLac->GetErrPro();
        if (bestErr > err || i == 0) {
            bestErr = err;
            pBestLac = pLac;
            bestSizeGain = pResubLac->GetSizeGain();
        }
        else if (bestErr == err) {
            ll sizeGain = pResubLac->GetSizeGain();
            if (sizeGain > bestSizeGain) {
                pBestLac = pLac;
                bestSizeGain = sizeGain;
            }
        }
    }
    cout << "bestsizeGain = " << bestSizeGain << endl;
    return pBestLac;
}


extern "C" {
    Vec_Ptr_t * Abc_MfsWinMarkTfi(Abc_Obj_t * pNode);
    void Abc_MfsWinSweepLeafTfo_rec(Abc_Obj_t * pObj, int nLevelLimit);
}
void LACMan::GetDivs(Abc_Obj_t * pNode, int nLevDivMax, RETURN_VAR IntVect& divs) {
    const int nWinMax = 300;
    const int nFanoutsMax = 30;
    divs.clear();
    Vec_Ptr_t * vCone, * vDivs;
    Abc_Obj_t * pObj, * pFanout, * pFanin;
    int k, f, m;
    int nDivsPlus = 0, nTrueSupp;

    // mark the TFI with the current trav ID
    Abc_NtkIncrementTravId( pNode->pNtk );
    vCone = Abc_MfsWinMarkTfi( pNode );

    // count the number of PIs
    nTrueSupp = 0;
    Vec_PtrForEachEntry( Abc_Obj_t *, vCone, pObj, k )
        nTrueSupp += Abc_ObjIsCi(pObj);
//    printf( "%d(%d) ", Vec_PtrSize(p->vSupp), m );

    // mark with the current trav ID those nodes that should not be divisors:
    // (1) the node and its TFO
    // (2) the MFFC of the node
    // (3) the node's fanins (these are treated as a special case)
    Abc_NtkIncrementTravId( pNode->pNtk );
    Abc_MfsWinSweepLeafTfo_rec( pNode, nLevDivMax );
//    Abc_MfsWinVisitMffc( pNode );
    Abc_ObjForEachFanin( pNode, pObj, k )
        Abc_NodeSetTravIdCurrent( pObj );

    // at this point the nodes are marked with two trav IDs:
    // nodes to be collected as divisors are marked with previous trav ID
    // nodes to be avoided as divisors are marked with current trav ID

    // start collecting the divisors
    vDivs = Vec_PtrAlloc( nWinMax );
    Vec_PtrForEachEntry( Abc_Obj_t *, vCone, pObj, k )
    {
        if ( !Abc_NodeIsTravIdPrevious(pObj) )
            continue;
        if ( (int)pObj->Level > nLevDivMax )
            continue;
        Vec_PtrPush( vDivs, pObj );
        if ( Vec_PtrSize(vDivs) >= nWinMax )
            break;
    }
    Vec_PtrFree( vCone );

    // explore the fanouts of already collected divisors
    if ( Vec_PtrSize(vDivs) < nWinMax )
    Vec_PtrForEachEntry( Abc_Obj_t *, vDivs, pObj, k )
    {
        // consider fanouts of this node
        Abc_ObjForEachFanout( pObj, pFanout, f )
        {
            // stop if there are too many fanouts
            if ( nFanoutsMax && f > nFanoutsMax )
                break;
            // skip nodes that are already added
            if ( Abc_NodeIsTravIdPrevious(pFanout) )
                continue;
            // skip nodes in the TFO or in the MFFC of node
            if ( Abc_NodeIsTravIdCurrent(pFanout) )
                continue;
            // skip COs
            if ( !Abc_ObjIsNode(pFanout) )
                continue;
            // skip nodes with large level
            if ( (int)pFanout->Level > nLevDivMax )
                continue;
            // skip nodes whose fanins are not divisors  -- here we skip more than we need to skip!!! (revise later)  August 7, 2009
            Abc_ObjForEachFanin( pFanout, pFanin, m )
                if ( !Abc_NodeIsTravIdPrevious(pFanin) )
                    break;
            if ( m < Abc_ObjFaninNum(pFanout) )
                continue;
            // make sure this divisor in not among the nodes
//            Vec_PtrForEachEntry( Abc_Obj_t *, p->vNodes, pFanin, m )
//                assert( pFanout != pFanin );
            // add the node to the divisors
            Vec_PtrPush( vDivs, pFanout );
            // Vec_PtrPushUnique( p->vNodes, pFanout );
            Abc_NodeSetTravIdPrevious( pFanout );
            nDivsPlus++;
            if ( Vec_PtrSize(vDivs) >= nWinMax )
                break;
        }
        if ( Vec_PtrSize(vDivs) >= nWinMax )
            break;
    }

    // sort the divisors by level in the increasing order
    Vec_PtrSort( vDivs, (int (*)(const void *, const void *))Abc_NodeCompareLevelsIncrease );

    // add the fanins of the node
    Abc_ObjForEachFanin( pNode, pFanin, k )
        Vec_PtrPush( vDivs, pFanin );
    
    divs.reserve(Vec_PtrSize(vDivs));
    Vec_PtrForEachEntry(Abc_Obj_t *, vDivs, pObj, k)
        divs.emplace_back(pObj->Id);

    // clean up
    Vec_PtrFree(vDivs);
}


static void GetDivs4NodeRec(Abc_Obj_t* pNode, RETURN_VAR IntVect& divs, BitVect& visited) {
    visited[pNode->Id] = 1;
    Abc_Obj_t* pFanin = nullptr;
    int i = 0;
    Abc_ObjForEachFanin(pNode, pFanin, i) {
        if (!visited[pFanin->Id])
            GetDivs4NodeRec(pFanin, divs, visited);
    }
    divs.emplace_back(pNode->Id);
}


static void GetDivs4Node(Abc_Obj_t* pNode, RETURN_VAR IntVect& divs, BitVect& visited) {
    assert(Abc_ObjIsNode(pNode));
    divs.clear();
    visited.reset();
    Abc_Obj_t* pFanin = nullptr;
    int i = 0;
    Abc_ObjForEachFanin(pNode, pFanin, i) {
        if (!visited[pFanin->Id])
            GetDivs4NodeRec(pFanin, divs, visited);
    }
}


// // multiple-thread lock
// static std::mutex mtx;


// static void GetDivs4SomeNodes(NetMan& net, IntVect& targIds, RETURN_VAR Int2DVect& divs4Nodes, boost::timer::progress_display& pd, int start, int end) {
//     for (int i = start; i < end; ++i) {
//         int targId = targIds[i];
//         auto pTarget = net.GetObj(targId);
//         BitVect visited(net.GetIdMaxPlus1(), 0);
//         GetDivs4Node(pTarget, divs4Nodes[targId], visited);
//         std::unique_lock<std::mutex> lock(mtx);
//         ++pd;
//         lock.unlock();
//     }
// }


// void LACMan::GetDivsParallelly(NetMan& net, IntVect& targIds, RETURN_VAR Int2DVect& divs4Nodes) {
//     cout << "collecting divisors" << endl;

//     // multi-thread parameters
//     int targetNodeNum = targIds.size();
//     int realThread = min(targetNodeNum, nThread);
//     assert(realThread > 0);
//     cout << "real thread number: " << realThread << endl;
//     int chunkSize = targetNodeNum / realThread;
//     int remainder = targetNodeNum % realThread;

//     // multi-thread calculation
//     divs4Nodes.resize(net.GetIdMaxPlus1());
//     boost::timer::progress_display pd(targetNodeNum);
//     vector<thread> threads;
//     int start = 0;
//     for (int i = 0; i < realThread; ++i) {
//         int end = start + chunkSize + (i < remainder? 1: 0);
//         threads.emplace_back(GetDivs4SomeNodes, std::ref(net), std::ref(targIds), std::ref(divs4Nodes), std::ref(pd), start, end);
//         start = end;
//     }
//     for (auto& thread: threads)
//         thread.join();
// }