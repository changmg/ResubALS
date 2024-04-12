#include "error.h"


using namespace std;
using namespace abc;
using namespace boost;


ErrMan::ErrMan(NetMan & netMan0, NetMan & netMan1, unsigned _seed, ll n_frame, DISTR_TYPE distr_type):
    net0(netMan0), net1(netMan1), pSmlt0(nullptr), pSmlt1(nullptr), seed(_seed), nFrame(n_frame), distrType(distr_type) {
    assert(IsPIOSame(net0, net1));
}


void ErrMan::InitForStatErr() {
    if (pSmlt0 != nullptr || pSmlt1 != nullptr) {
        assert(pSmlt0 != nullptr && pSmlt1 != nullptr);
        return;
    }
    pSmlt0 = make_shared <Simulator> (net0, seed, nFrame);
    pSmlt1 = make_shared <Simulator> (net1, seed, nFrame); 
    if (distrType == DISTR_TYPE::UNIF) {
        pSmlt0->InpUnifFast();
        pSmlt1->InpUnifFast();
    }
    else if (distrType == DISTR_TYPE::ENUM) {
        pSmlt0->InpEnum();
        pSmlt1->InpEnum();
    }
    else if (distrType == DISTR_TYPE::MIX) {
        pSmlt0->InpMix();
        pSmlt1->InpMix();
    }
    // else if (distrType == DISTR_TYPE::SELF) {
    //     pSmlt0->InpSelf(selfDefDistr);
    //     pSmlt1->InpSelf(selfDefDistr);
    // }
    else
        assert(0);
    pSmlt0->Sim();
    pSmlt1->Sim();
}


double ErrMan::CalcErrRate() {
    InitForStatErr();
    return pSmlt0->GetErrRate(*pSmlt1);
}


double ErrMan::CalcMeanErrDist(bool isSign) {
    InitForStatErr();
    return pSmlt0->GetMeanErrDist(*pSmlt1, isSign);
}


double ErrMan::CalcMeanSquareErr(bool isSign) {
    InitForStatErr();
    return pSmlt0->GetMeanSquareErr(*pSmlt1, isSign);
}

double ErrMan::CalcMeanHammDist() {
    InitForStatErr();
    return pSmlt0->GetMeanHammDist(*pSmlt1);
}


double ErrMan::CalcSigNoiseRat(bool isSign) {
    InitForStatErr();
    return pSmlt0->GetSigNoiseRat(*pSmlt1, isSign);
}


// double ErrMan::CalcSelfDefErr(bool isSign, const string & selfDefMetr) {
//     InitForStatErr();
//     return pSmlt0->GetSelfDefErr(*pSmlt1, isSign, selfDefMetr);
// }


ull ErrMan::CalcMaxErrDist(bool isSign) {
    assert(net0.GetPiNum() <= 60 && net1.GetPiNum() <= 60);
    assert(isSign == 0);
    auto pNtk0 = Abc_NtkDup(net0.GetNet());
    auto pNtk1 = Abc_NtkDup(net1.GetNet());
    ull res = GET_MEM(pNtk0, pNtk1);
    Abc_NtkDelete(pNtk0);
    Abc_NtkDelete(pNtk1);
    return res;
}


ull ErrMan::GET_MEM(Abc_Ntk_t * pNtk1, Abc_Ntk_t * pNtk2) {
    ull pTemp = 0;
    ll fRemove1, fRemove2;
    assert( Abc_NtkHasOnlyLatchBoxes(pNtk1) );
    assert( Abc_NtkHasOnlyLatchBoxes(pNtk2) );
    fRemove1 = (!Abc_NtkIsStrash(pNtk1) || Abc_NtkGetChoiceNum(pNtk1)) && (pNtk1 = Abc_NtkStrash(pNtk1, 0, 0, 0));
    fRemove2 = (!Abc_NtkIsStrash(pNtk2) || Abc_NtkGetChoiceNum(pNtk2)) && (pNtk2 = Abc_NtkStrash(pNtk2, 0, 0, 0));
    if ( pNtk1 && pNtk2 )
        pTemp = NtkMiterComp( pNtk1, pNtk2 );
    if ( fRemove1 )  Abc_NtkDelete( pNtk1 );
    if ( fRemove2 )  Abc_NtkDelete( pNtk2 );
    return pTemp;
}


static void Abc_NtkMiterPrepare( Abc_Ntk_t * pNtk1, Abc_Ntk_t * pNtk2, Abc_Ntk_t * pNtkMiter) {
    Abc_Obj_t * pObj, * pObjNew; ll i;
    Abc_AigConst1(pNtk1)->pCopy = Abc_AigConst1(pNtkMiter);
    Abc_AigConst1(pNtk2)->pCopy = Abc_AigConst1(pNtkMiter);

    // create new PIs and remember them in the old PIs
    Abc_NtkForEachPi(pNtk1, pObj, i)
    {
        pObjNew = Abc_NtkCreatePi(pNtkMiter);
        // remember this PI in the old PIs
        pObj->pCopy = pObjNew;
        pObj = Abc_NtkPi(pNtk2, i);  
        pObj->pCopy = pObjNew;
            // add name
        Abc_ObjAssignName( pObjNew, Abc_ObjName(pObj), NULL );
    }

        // pObjNew = Abc_NtkCreatePo(pNtkMiter);
        // Abc_ObjAssignName(pObjNew, "miter", NULL);

    Abc_NtkForEachLatch( pNtk1, pObj, i )
    {
        pObjNew = Abc_NtkDupBox( pNtkMiter, pObj, 0 );
        // add names
        Abc_ObjAssignName( pObjNew, Abc_ObjName(pObj), "_1" );
        Abc_ObjAssignName( Abc_ObjFanin0(pObjNew),  Abc_ObjName(Abc_ObjFanin0(pObj)), "_1" );
        Abc_ObjAssignName( Abc_ObjFanout0(pObjNew), Abc_ObjName(Abc_ObjFanout0(pObj)), "_1" );
    }
    Abc_NtkForEachLatch( pNtk2, pObj, i )
    {
        pObjNew = Abc_NtkDupBox( pNtkMiter, pObj, 0 );
        // add name
        Abc_ObjAssignName( pObjNew, Abc_ObjName(pObj), "_2" );
        Abc_ObjAssignName( Abc_ObjFanin0(pObjNew),  Abc_ObjName(Abc_ObjFanin0(pObj)), "_2" );
        Abc_ObjAssignName( Abc_ObjFanout0(pObjNew), Abc_ObjName(Abc_ObjFanout0(pObj)), "_2" );
    }
}


static void Abc_NtkMiterAddOne( Abc_Ntk_t * pNtk, Abc_Ntk_t * pNtkMiter ) {
    Abc_Obj_t * pNode;
    ll i;
    assert( Abc_NtkIsDfsOrdered(pNtk) );
    Abc_AigForEachAnd( pNtk, pNode, i )
        pNode->pCopy = Abc_AigAnd( (Abc_Aig_t *)pNtkMiter->pManFunc, Abc_ObjChild0Copy(pNode), Abc_ObjChild1Copy(pNode) );
}


ull ErrMan::NtkMiterComp(Abc_Ntk_t * pNtk1, Abc_Ntk_t * pNtk2) {
    const ll fImplic = 0, fComb = 0, nPartSize = 0, fMulti = 0;
    char Buffer[1000];
    Abc_Ntk_t * pNtkMiter;

    assert( Abc_NtkIsStrash(pNtk1) );
    assert( Abc_NtkIsStrash(pNtk2) );

    // start the new network
    pNtkMiter = Abc_NtkAlloc(ABC_NTK_STRASH, ABC_FUNC_AIG, 1);
    sprintf( Buffer, "%s_%s_miter", pNtk1->pName, pNtk2->pName );
    pNtkMiter->pName = Extra_UtilStrsav(Buffer);

    // perform strashing
    Abc_NtkMiterPrepare( pNtk1, pNtk2, pNtkMiter );
    Abc_NtkMiterAddOne( pNtk1, pNtkMiter ); 
    Abc_NtkMiterAddOne( pNtk2, pNtkMiter );
    ull x = NtkMiterFinalize( pNtk1, pNtk2, pNtkMiter, fComb, nPartSize, fImplic, fMulti );
    
    Abc_AigCleanup((Abc_Aig_t *)pNtkMiter->pManFunc);

    // make sure that everything is okay
    Abc_NtkDelete( pNtkMiter );
    return x;
}


ull ErrMan::NtkMiterFinalize( Abc_Ntk_t * pNtk1, Abc_Ntk_t * pNtk2, Abc_Ntk_t * pNtkMiter, ll fComb, ll nPartSize, ll fImplic, ll fMulti ) {
    Vec_Ptr_t * vPairs;
    Abc_Obj_t * pNode;
    ll i;
    assert( nPartSize == 0 || fMulti == 0 );
    // collect the PO pairs from both networks
    vPairs = Vec_PtrAlloc( 100 );
    // collect the PO nodes for the miter
    Abc_NtkForEachPo( pNtk1, pNode, i )
    {
        Vec_PtrPush( vPairs, Abc_ObjChild0Copy(pNode) );
        pNode = Abc_NtkPo( pNtk2, i );
        Vec_PtrPush( vPairs, Abc_ObjChild0Copy(pNode) );
    }
    Abc_NtkForEachLatch( pNtk1, pNode, i )
        Abc_ObjAddFanin( Abc_ObjFanin0(pNode)->pCopy, Abc_ObjChild0Copy(Abc_ObjFanin0(pNode)) );
    Abc_NtkForEachLatch( pNtk2, pNode, i )
        Abc_ObjAddFanin( Abc_ObjFanin0(pNode)->pCopy, Abc_ObjChild0Copy(Abc_ObjFanin0(pNode)) );
    
    // add the miter
    // Abc_Obj_t * NewPO = Abc_NtkPo(pNtkMiter, 0);
    // Abc_ObjAddFanin( NewPO, pMiter );
    
    Abc_Obj_t ** X = new Abc_Obj_t *[vPairs->nSize / 2];
    Abc_Obj_t ** Y = new Abc_Obj_t *[vPairs->nSize / 2];
    for(i = 0; i < vPairs->nSize; i += 2)
    {
        X[i / 2] = (Abc_Obj_t *) vPairs->pArray[i];
        Y[i / 2] = (Abc_Obj_t *) vPairs->pArray[i + 1];
    }
            // for(i=0; i<vPairs->nSize; i += 2)
            // {
            //     Abc_Obj_t * NewPo = Abc_NtkCreatePo(pNtkMiter);
            //     Abc_ObjAddFanin(NewPo, (Abc_Obj_t *)vPairs->pArray[i]);
            // }
            // for(i=0; i<vPairs->nSize; i += 2)
            // {
            //     Abc_Obj_t * NewPo = Abc_NtkCreatePo(pNtkMiter);
            //     Abc_ObjAddFanin(NewPo, (Abc_Obj_t *)vPairs->pArray[i+1]);
            // }
            // Ckt_WriteBlif(pNtkMiter, "first.blif"); 
    Abc_Obj_t ** R = X_subtract_Y_abs(pNtkMiter, X, Y, vPairs->nSize / 2);     
    delete[] X;
    delete[] Y;
    ull res = GETMEM(pNtkMiter, R, vPairs->nSize / 2);
    delete[] R;
            // cout << res << endl;

    Vec_PtrFree( vPairs );
    return res;
}


ull ErrMan::GETMEM(Abc_Ntk_t * pNtk, Abc_Obj_t *R[], ll n) {
    Abc_Obj_t * ConstNode[2];
    Abc_Obj_t ** mem = new Abc_Obj_t *[n];
    ConstNode[1] = Abc_AigConst1(pNtk);
    ConstNode[0] = Abc_ObjNot(ConstNode[1]);
    bool CurrentState = true, PreviousState = true;
    ll round = 0;
    std::vector<ll> ConstMEMInput(n, 0);
    ConstMEMInput[n - 1] = 1;
    Abc_Obj_t * tempR;
    while(true)
    {
        // for(ll k = 0; k < n; k++) cout << ConstMEMInput[n-1-k];
        // cout << round << endl;
        Abc_Obj_t * result = Abc_NtkCreatePo(pNtk);
        for(ll k = 0; k < n; k++) mem[k] = ConstNode[ConstMEMInput[k]];
        // compXY = (R > mem)
        tempR = X_lt_Y(pNtk, mem, R, n);
        Abc_ObjAddFanin(result, tempR);
        // Ckt_WriteBlif(pNtk, "merge_SAT.blif");
        PreviousState = CurrentState;
        CurrentState = SATSolver(pNtk);
        if(round == n) break;
        if(CurrentState)
        {
            if(round < n - 1) ConstMEMInput[n - 1 - ++round] = 1;
            else break;
        }
        else 
        {
            if (round < n-1) {ConstMEMInput[n - 1 - round++] = 0; ConstMEMInput[n - 1 - round] = 1;} 
            else ConstMEMInput[n - 1 - round++] = 0;
        }
        Abc_NtkDeleteObj(result);
    }
    delete[] mem; 
    ull res = 0;
    for(ll k = 0; k < n; k++) {
        res <<= 1;
        res += ConstMEMInput[n-1-k];
        // cout << res << endl;
    }
    if(CurrentState) res++; 
    return res;
}


Abc_Obj_t ** ErrMan::X_subtract_Y_abs(Abc_Ntk_t * pNtk, Abc_Obj_t * X[], Abc_Obj_t * Y[], ll n) {
    if(n <= 0) return nullptr;
    Abc_Obj_t ** R = new Abc_Obj_t *[n];
    Abc_Obj_t ** Cout = new Abc_Obj_t *[n];
    R[0] = Abc_AigXor((Abc_Aig_t *) pNtk->pManFunc, X[0], Y[0]);
    Cout[0] = Abc_AigAnd((Abc_Aig_t *) pNtk->pManFunc, Abc_ObjNot(X[0]), Y[0]);
    for(ll i=1; i<n; i++)
    {
        R[i] = Abc_AigXor((Abc_Aig_t *) pNtk->pManFunc, 
            Abc_AigXor((Abc_Aig_t *) pNtk->pManFunc, Cout[i-1], X[i]), 
            Y[i]
        );
        Abc_Obj_t * temp1 = Abc_AigAnd((Abc_Aig_t *) pNtk->pManFunc, Abc_ObjNot(X[i]), Y[i]);
        Abc_Obj_t * temp2 = Abc_AigOr((Abc_Aig_t *) pNtk->pManFunc, Abc_ObjNot(X[i]), Y[i]);
        Cout[i] = Abc_AigOr((Abc_Aig_t *) pNtk->pManFunc,
            Abc_AigAnd((Abc_Aig_t *) pNtk->pManFunc, temp1, Abc_ObjNot(Cout[i-1])),
            Abc_AigAnd((Abc_Aig_t *) pNtk->pManFunc, temp2, Cout[i-1])
        );
    }

    // 2's complement of R ( = R' + 1 )
    Abc_Obj_t ** R_2Complement = new Abc_Obj_t *[n];
    Abc_Obj_t ** Cout_2 = new Abc_Obj_t * [n];
    R_2Complement[0] = R[0];
    Cout_2[0] = Abc_ObjNot(R[0]);
    for(ll k=1; k<n; k++)
    {
        R_2Complement[k] = Abc_AigXor((Abc_Aig_t *) pNtk->pManFunc, Cout_2[k-1], Abc_ObjNot(R[k]));
        Cout_2[k] = Abc_AigAnd((Abc_Aig_t *) pNtk->pManFunc, Cout_2[k-1], Abc_ObjNot(R[k]));
    }
    Abc_Obj_t ** res = new Abc_Obj_t *[n];
    for(ll k=0; k<n; k++)
        res[k] = Abc_AigMux((Abc_Aig_t *) pNtk->pManFunc, Cout[n-1], R_2Complement[k], R[k]);
    delete[] R_2Complement;
    delete[] Cout;
    delete[] Cout_2;
    delete[] R;
    return res;
}


Abc_Obj_t * ErrMan::X_lt_Y(Abc_Ntk_t * pNtk, Abc_Obj_t * X[], Abc_Obj_t * Y[], ll n) {
    /* implementation */
    if(n <= 0) return nullptr;
    std::vector<Abc_Obj_t *> PreBitsComp(n, nullptr);
    PreBitsComp[0] = Abc_AigAnd((Abc_Aig_t *) pNtk->pManFunc, Abc_ObjNot(X[0]), Y[0]);
    for(ll i=1; i<n; i++)
    {
        PreBitsComp[i] = 
            Abc_AigOr((Abc_Aig_t *) pNtk->pManFunc, 
                Abc_AigAnd((Abc_Aig_t *) pNtk->pManFunc, Abc_ObjNot(X[i]), Y[i]),
                Abc_AigAnd((Abc_Aig_t *) pNtk->pManFunc, PreBitsComp[i-1],  
                    Abc_ObjNot(Abc_AigXor((Abc_Aig_t *) pNtk->pManFunc, X[i], Y[i]))
                )
            );
    }
    return PreBitsComp[n - 1];
}


bool ErrMan::SATSolver(Abc_Ntk_t * pNtk) {
    ll RetValue = -1;
    ll fVerbose = 0;
    ll nConfLimit = 0;
    ll nInsLimit = 0;
    assert(pNtk != nullptr);
    assert( Abc_NtkIsStrash(pNtk) );
    RetValue = Abc_NtkMiterSat( pNtk, (ABC_INT64_T)nConfLimit, (ABC_INT64_T)nInsLimit, fVerbose, NULL, NULL );
    if (pNtk->pModel != nullptr)
        ABC_FREE(pNtk->pModel);
    if (RetValue == -1)
        assert(0);
    else if (RetValue == 0)
        return true;
    else
        return false;
    return 0;
}


// double CalcErrPro(NetMan& net0, NetMan& net1, bool isSign, unsigned seed, ll nFrame, METR_TYPE metrType, DISTR_TYPE distrType) {
//     ErrManPro errMan(net0, net1, isSign, seed, nFrame, metrType, distrType);
//     errMan.InitMit();
//     auto err = errMan.CalcErr();
//     return err;
// }


double CalcErr(NetMan & netMan0, NetMan & netMan1, bool isSign, unsigned seed, ll nFrame, METR_TYPE metrType, DISTR_TYPE distrType) {
    ErrMan errMan(netMan0, netMan1, seed, nFrame, distrType);
    if (metrType == METR_TYPE::ER)
        return errMan.CalcErrRate();
    else if (metrType == METR_TYPE::MED)
        return errMan.CalcMeanErrDist(isSign);
    else if (metrType == METR_TYPE::MSE)
        return errMan.CalcMeanSquareErr(isSign);
    else if (metrType == METR_TYPE::MHD)
        return errMan.CalcMeanHammDist();
    else if (metrType == METR_TYPE::SNR)
        return errMan.CalcSigNoiseRat(isSign);
    else if (metrType == METR_TYPE::MAXED)
        return errMan.CalcMaxErrDist(isSign);
    // else if (metrType == METR_TYPE::SELF)
    //     return errMan.CalcSelfDefErr(isSign, selfDefMetr);
    else {
        assert(0);
        return 0;
    }
}


// double GetMSEFromSNR(NetMan & net, bool isSign, unsigned seed, ll nFrame, DISTR_TYPE distrType, double snr) {
//     Simulator smlt(net, seed, nFrame);
//     if (distrType == DISTR_TYPE::ENUM)
//         smlt.InpEnum();
//     else if (distrType == DISTR_TYPE::UNIF)
//         smlt.InpUnifFast();
//     else if (distrType == DISTR_TYPE::MIX)
//         smlt.InpMix();
//     else
//         assert(0);
//     smlt.Sim();
//     BigInt sumAcc2 = 0;
//     for (ll i = 0; i < nFrame; ++i) {
//         BigInt accOut = smlt.GetOutpPro(i, isSign);
//         sumAcc2 += accOut * accOut;
//     }
//     return static_cast <double> (BigFlt(sumAcc2) / BigFlt(nFrame) / BigFlt(pow(BigFlt(10),  BigFlt(snr) / 10)));
// }


// multiple-thread lock
std::mutex mtx;


static void CalcSomeLACLocalERs(Simulator& appSmlt, LACMan& lacMan, boost::dynamic_bitset<ull>& IsErroneousPattern, timer::progress_display& pd, ll startIndex, ll endIndex) {
    for (ll lacId = startIndex; lacId < endIndex; ++lacId) {
        auto pLac = lacMan.GetLac(lacId);
        ll targId = pLac->GetTargId();

        // calculate $\partial n / \partial LAC$
        auto & specLac = *dynamic_pointer_cast <ResubLAC>(pLac);
        auto divIds = specLac.GetDivIds(); 
        auto sop = specLac.GetSop();
        BitVect newValue(appSmlt.GetFrameNumb(), 0);
        GetNewValue(appSmlt, divIds, sop, newValue);
        auto isChanged = (*appSmlt.GetDat(targId)) ^ newValue;
        auto estimatedDiff = isChanged & (~IsErroneousPattern);
        BigInt er = estimatedDiff.count();
        pLac->SetErrPro(er);

        // update progress
        std::unique_lock<std::mutex> lock(mtx);
        ++pd;
        lock.unlock();
    }
}


void VECBEEMan::ComputeLocalErrorRate(Simulator& accSmlt, Simulator& appSmlt, LACMan & lacMan) {
    // check
    assert(IsPIOSame(accSmlt, appSmlt));
    assert(lacType == LAC_TYPE::RESUB);
    assert(metrType == METR_TYPE::ER);
    assert((nFrame & 63) == 0);

    // get erroneous patterns
    dynamic_bitset<ull> IsErroneousPattern(nFrame, 0);
    for (ll i = 0; i < accSmlt.GetPoNum(); ++i) {
        ll poIdAcc = accSmlt.GetPoId(i);
        ll poIdApp = appSmlt.GetPoId(i);
        IsErroneousPattern |= (*accSmlt.GetDat(poIdAcc) ^ *appSmlt.GetDat(poIdApp));
    }

    // compute local ER
    cout << "calculating LACs' local ERs" << endl;
    int lacNum = lacMan.GetLacNum();
    int realThread = min(nThread, lacNum);
    assert(realThread > 0);
    cout << "using " << realThread << " threads" << endl;
    ll chunkSize = lacNum / realThread;
    ll remainder = lacNum % realThread;
    timer::progress_display pd(lacMan.GetLacNum());
    vector<thread> threads;
    ll start = 0;
    for (ll i = 0; i < realThread; ++i) {
        ll end = start + chunkSize + (i < remainder? 1: 0);
        threads.emplace_back(CalcSomeLACLocalERs, std::ref(appSmlt), std::ref(lacMan), std::ref(IsErroneousPattern), std::ref(pd), start, end);
        start = end;
    }
    for (auto& thread: threads)
        thread.join();
}


static void EstSomeLACErrorBounds(Simulator& miterSmlt, LACMan& lacMan, const IntVect& appId2MiterId, vector<vector<BitVect>>& bdPo2Nodes, timer::progress_display& pd, ll startIndex, ll endIndex) {
    int nPo = miterSmlt.GetPoNum();
    for (ll lacId = startIndex; lacId < endIndex; ++lacId) {
        auto pLac = lacMan.GetResubLac(lacId);
        // cout << "deal with LAC " << lacId << endl;
        // pLac->Print();
        int targIdInApp = pLac->GetTargId();
        int targId = appId2MiterId[targIdInApp];
        assert(targId != -1);
        // calculate $\partial n / \partial LAC$
        // auto & specLac = *dynamic_pointer_cast <ResubLAC>(pLac);
        auto divIdsInApp = pLac->GetDivIds(); 
        IntVect divIds; divIds.clear();
        for (auto divIdInApp: divIdsInApp) {
            assert(appId2MiterId[divIdInApp] != -1);
            divIds.emplace_back(appId2MiterId[divIdInApp]);
        }
        auto sop = pLac->GetSop();
        BitVect newValue(miterSmlt.GetFrameNumb(), 0);
        GetNewValue(miterSmlt, divIds, sop, newValue);
        auto isChanged = (*miterSmlt.GetDat(targId)) ^ newValue;
        // cout << miterSmlt.GetName(targId) << endl;
        // cout << "current value: " << *miterSmlt.GetDat(targId) << endl;
        // cout << "new value:" << newValue << endl;
        BigInt deltaErrorBound(0);
        for (int k = nPo - 1; k >= 0; --k) {
            deltaErrorBound <<= 1;
            // auto affectPoK = isChanged & bdPo2Nodes[k][targId];
            // the k-th miter PO is currently 0, and the LAC influences the k-th miter PO
            auto temp = (isChanged & bdPo2Nodes[k][targId]) & (~*miterSmlt.GetDat(miterSmlt.GetPoId(k)));
            deltaErrorBound += temp.count();
        }
        pLac->SetErrPro(deltaErrorBound);

        // update progress
        std::unique_lock<std::mutex> lock(mtx);
        ++pd;
        lock.unlock();
    }
}


void VECBEEMan::EstimateErrorBoundForEachLAC(Simulator& miterSmlt, LACMan& lacMan, const BigInt& upperBound, bool enableFastErrEst, const IntVect& miterId2AppId, const IntVect& appId2MiterId) {
    assert(lacType == LAC_TYPE::RESUB);
    // assert(metrType == METR_TYPE::MED);
    
    auto topoNodes = miterSmlt.TopoSort();
    NetMan& miterNet = miterSmlt;
    // miterNet.WriteBlif("./tmp/miter.blif");
    // miterNet.WriteDot("./tmp/miter.dot");
    // PrintVect(topoNodes, "\n");
    if (enableFastErrEst)
        FindAppDisjCut(miterNet);
    else
        FindDisjCut(miterNet, topoNodes);
    CalcBoolDiffCut2Node(miterSmlt, topoNodes);
    CalcBoolDiffPo2NodeParallelly(miterSmlt, topoNodes);

    cout << "estimating LACs' MED bounds" << endl;
    int lacNum = lacMan.GetLacNum();
    int realThread = min(nThread, lacNum);
    assert(realThread > 0);
    cout << "using " << realThread << " threads" << endl;
    int chunkSize = lacNum / realThread;
    int remainder = lacNum % realThread;
    timer::progress_display pd(lacMan.GetLacNum());
    vector<thread> threads;
    ll start = 0;
    for (ll i = 0; i < realThread; ++i) {
        ll end = start + chunkSize + (i < remainder? 1: 0);
        threads.emplace_back(EstSomeLACErrorBounds, std::ref(miterSmlt), std::ref(lacMan), std::ref(appId2MiterId), std::ref(bdPo2Nodes), std::ref(pd), start, end);
        start = end;
    }
    for (auto& thread: threads)
        thread.join();

    // for (ll lacId = 0; lacId < lacMan.GetLacNum(); ++lacId) {
    //     auto pLac = lacMan.GetResubLac(lacId);
    //     pLac->Print();
    // }
}


// void VECBEEMan::EstimateErrorBoundForEachLAC(Simulator& accSmlt, Simulator& appSmlt, LACMan& lacMan, const BigInt& upperBound, bool enableFastErrEst) {
//     assert(IsPIOSame(accSmlt, appSmlt));
//     assert(lacType == LAC_TYPE::RESUB);
//     assert((nFrame & 63) == 0);
//     if (metrType == METR_TYPE::ER)
//         EstimateERBound(accSmlt, appSmlt, lacMan);
//     else if (metrType == METR_TYPE::MED) {
//         auto topoNodes = appSmlt.TopoSort();
//         NetMan& appNet = appSmlt;
//         if (enableFastErrEst)
//             // FindAppDisjCut(appNet);
//             FindAppDisjCutNew(appNet, topoNodes);
//         else
//             FindDisjCut(appNet, topoNodes);
//         CalcBoolDiffCut2NodeParallelly(appSmlt, topoNodes);
//         CalcBoolDiffPo2NodeParallelly(appSmlt, topoNodes);
//         EstimateMEDBound(accSmlt, appSmlt, lacMan);
//     }
//     else
//         assert(0);
// }


// void VECBEEMan::ComputeBooleanDifferenceOfPos2Nodes(Simulator& appSmlt, AbcObjVect& topoNodes, bool enableFastErrEst) {
//     assert((nFrame & 63) == 0);
//     NetMan& appNet = appSmlt;
//     if (enableFastErrEst)
//         FindAppDisjCut(appNet);
//     else
//         FindDisjCut(appNet, topoNodes);
//     CalcBoolDiffCut2Node(appSmlt, topoNodes); // perf: multiple threads
//     CalcBoolDiffPo2NodeParallelly(appSmlt, topoNodes);
// }


static void EstSomeLACERBounds(Simulator& appSmlt, LACMan& lacMan, boost::dynamic_bitset<ull>& isCorrectPattern, timer::progress_display& pd, ll startIndex, ll endIndex) {
    for (ll lacId = startIndex; lacId < endIndex; ++lacId) {
        auto pLac = lacMan.GetLac(lacId);
        ll targId = pLac->GetTargId();

        // calculate $\partial n / \partial LAC$
        auto & specLac = *dynamic_pointer_cast <ResubLAC>(pLac);
        auto divIds = specLac.GetDivIds(); 
        auto sop = specLac.GetSop();
        BitVect newValue(appSmlt.GetFrameNumb(), 0);
        GetNewValue(appSmlt, divIds, sop, newValue);
        auto isChanged = (*appSmlt.GetDat(targId)) ^ newValue;
        auto estimatedDiff = isChanged & isCorrectPattern;
        BigInt er = estimatedDiff.count();
        pLac->SetErrPro(er);

        // update progress
        std::unique_lock<std::mutex> lock(mtx);
        ++pd;
        lock.unlock();
    }
}


void VECBEEMan::EstimateERBound(Simulator& accSmlt, Simulator& appSmlt, LACMan& lacMan) {
    // get correct patterns
    dynamic_bitset<ull> isCorrectPattern(nFrame);
    isCorrectPattern.set();
    for (ll i = 0; i < accSmlt.GetPoNum(); ++i) {
        ll poIdAcc = accSmlt.GetPoId(i);
        ll poIdApp = appSmlt.GetPoId(i);
        isCorrectPattern &= (~(*accSmlt.GetDat(poIdAcc) ^ *appSmlt.GetDat(poIdApp)));
    }

    // estimate ER bound
    cout << "estimating LACs' ER bounds" << endl;
    int lacNum = lacMan.GetLacNum();
    int realThread = min(nThread, lacNum);
    assert(realThread > 0);
    cout << "using " << realThread << " threads" << endl;
    ll chunkSize = lacNum / realThread;
    ll remainder = lacNum % realThread;
    timer::progress_display pd(lacMan.GetLacNum());
    vector<thread> threads;
    ll start = 0;
    for (ll i = 0; i < realThread; ++i) {
        ll end = start + chunkSize + (i < remainder? 1: 0);
        threads.emplace_back(EstSomeLACERBounds, std::ref(appSmlt), std::ref(lacMan), std::ref(isCorrectPattern), std::ref(pd), start, end);
        start = end;
    }
    for (auto& thread: threads)
        thread.join();
}


static void EstSomeLACMEDBounds(Simulator& appSmlt, LACMan& lacMan, vector<vector<BitVect>>& bdPo2Nodes, timer::progress_display& pd, ll startIndex, ll endIndex) {
    int nPo = appSmlt.GetPoNum();
    for (ll lacId = startIndex; lacId < endIndex; ++lacId) {
        auto pLac = lacMan.GetLac(lacId);
        ll targId = pLac->GetTargId();

        // calculate $\partial n / \partial LAC$
        auto & specLac = *dynamic_pointer_cast <ResubLAC>(pLac);
        auto divIds = specLac.GetDivIds(); 
        auto sop = specLac.GetSop();
        dynamic_bitset <ull> newValue(appSmlt.GetFrameNumb(), 0);
        GetNewValue(appSmlt, divIds, sop, newValue);
        auto isChanged = (*appSmlt.GetDat(targId)) ^ newValue;
        BigInt deltaMEDBound(0);
        for (int k = nPo - 1; k >= 0; --k) {
            deltaMEDBound <<= 1;
            auto affectPoK = isChanged & bdPo2Nodes[k][targId];
            deltaMEDBound += affectPoK.count();
        }
        pLac->SetErrPro(deltaMEDBound);

        // update progress
        std::unique_lock<std::mutex> lock(mtx);
        ++pd;
        lock.unlock();
    }
}


void VECBEEMan::EstimateMEDBound(Simulator& accSmlt, Simulator& appSmlt, LACMan& lacMan) {
    cout << "estimating LACs' MED bounds" << endl;
    int lacNum = lacMan.GetLacNum();
    int realThread = min(nThread, lacNum);
    assert(realThread > 0);
    cout << "using " << realThread << " threads" << endl;
    ll chunkSize = lacNum / realThread;
    ll remainder = lacNum % realThread;
    timer::progress_display pd(lacMan.GetLacNum());
    vector<thread> threads;
    ll start = 0;
    for (ll i = 0; i < realThread; ++i) {
        ll end = start + chunkSize + (i < remainder? 1: 0);
        threads.emplace_back(EstSomeLACMEDBounds, std::ref(appSmlt), std::ref(lacMan), std::ref(bdPo2Nodes), std::ref(pd), start, end);
        start = end;
    }
    for (auto& thread: threads)
        thread.join();
}


void VECBEEMan::BatchErrEstPro(NetMan& accNet, NetMan& appNet, LACMan& lacMan, const BigInt& uppBound, bool enableFastErrEst, BigInt& backErrInt) {
    assert(IsPIOSame(accNet, appNet));
    assert(lacType == LAC_TYPE::RESUB);
    auto topoNodes = appNet.TopoSort();
    if (enableFastErrEst)
        FindAppDisjCut(appNet);
    else
        FindDisjCut(appNet, topoNodes);
    Simulator accSmlt(accNet, seed, nFrame);
    Simulator appSmlt(appNet, seed, nFrame);
    if (distrType == DISTR_TYPE::UNIF) {
        accSmlt.InpUnifFast();
        appSmlt.InpUnifFast();
    }
    else if (distrType == DISTR_TYPE::ENUM) {
        accSmlt.InpEnum();
        appSmlt.InpEnum();
    }
    else if (distrType == DISTR_TYPE::MIX) {
        accSmlt.InpMix();
        appSmlt.InpMix();
    }
    else
        assert(0);
    accSmlt.Sim();
    appSmlt.Sim();
    CalcBoolDiffCut2NodeParallelly(appSmlt, topoNodes);
    CalcBoolDiffPo2NodeParallelly(appSmlt, topoNodes);
    if (metrType == METR_TYPE::MED || metrType == METR_TYPE::MSE)
        CalcLACNonERErrsNew(accSmlt, appSmlt, lacMan, uppBound, backErrInt);
    else if (metrType == METR_TYPE::ER || metrType == METR_TYPE::MHD)
        CalcLACERErrs(accSmlt, appSmlt, lacMan);
}


void VECBEEMan::FindDisjCut(NetMan& net, AbcObjVect& topoNodes) {
    cout << "finding disjoint cuts" << endl;
    assert(disjCuts.empty());
    assert(cutNtks.empty());
    assert(topoNodes.size());
    assert(topoNodes[0]->pNtk == net.GetNet());

    // init
    cutNtks.resize(net.GetIdMaxPlus1());
    disjCuts.resize(net.GetIdMaxPlus1());
    poMarks.resize(net.GetIdMaxPlus1(), dynamic_bitset <ull>(net.GetPoNum(), 0));
    for (ll i = 0; i < net.GetIdMaxPlus1(); ++i) {
        if (net.GetObj(i) == nullptr)
            continue;
        poMarks[i].reset();
    }

    // update topo ids
    topoIds.resize(net.GetIdMaxPlus1());
    for (ll i = 0; i < topoNodes.size(); ++i)
        topoIds[topoNodes[i]->Id] = i;
    ll topoId = -1;
    for (ll i = 0; i < net.GetPiNum(); ++i)
        topoIds[net.GetPiId(i)] = topoId--;
    topoId = topoNodes.size();
    for (ll i = 0; i < net.GetPoNum(); ++i)
        topoIds[net.GetPoId(i)] = topoId++;

    // determine the POs that each node will affect
    for (ll i = 0; i < net.GetPoNum(); ++i)
        poMarks[net.GetPoId(i)].set(i);
    for (auto it = topoNodes.rbegin(); it != topoNodes.rend(); ++it) {
        auto pObj = *it;
        if (pObj == nullptr)
            continue;
        ll i = net.GetId(pObj);
        for (ll j = 0; j < net.GetFanoutNum(pObj); ++j)
            poMarks[i] |= poMarks[net.GetFanoutId(pObj, j)];
    }

    // collect disjoint cuts and the corresponding cut networks
    timer::progress_display pd(net.GetIdMaxPlus1());
    for (ll i = 0; i < net.GetIdMaxPlus1(); ++i) {
        auto pObj = net.GetObj(i);
        if (!net.IsNode(pObj)) {
            ++pd;
            continue;
        }
        // cout << "finding " << pObj << endl;
        Abc_NtkIncrementTravId(net.GetNet());
        FindDisjCutOfNode(pObj, disjCuts[i]);
        for (const auto & node: topoNodes) {
            if (Abc_NodeIsTravIdCurrent(node))
                cutNtks[i].emplace_back(node);
        }
        for (ll j = 0; j < net.GetPoNum(); ++j) {
            auto pPo = net.GetPo(j);
            if (Abc_NodeIsTravIdCurrent(pPo))
                cutNtks[i].emplace_back(pPo);
        }
        ++pd;
    }
}


void VECBEEMan::FindAppDisjCut(NetMan & net) {
    cout << "finding approximate disjoint cuts" << endl;
    assert(disjCuts.empty());
    assert(cutNtks.empty());

    // init
    cutNtks.resize(net.GetIdMaxPlus1());
    disjCuts.resize(net.GetIdMaxPlus1());

    // collect disjoint cuts and the corresponding cut networks
    for (int iNode = 0; iNode < net.GetIdMaxPlus1(); ++iNode) {
        if (!net.IsNode(iNode))
            continue;
        if (net.IsConst(iNode))
            continue;
        for (int iFanout = 0; iFanout < net.GetFanoutNum(iNode); ++iFanout) {
            auto pFanout = net.GetFanout(iNode, iFanout);
            disjCuts[iNode].emplace_back(pFanout);
            cutNtks[iNode].emplace_back(pFanout);
        }
    }
}


// void VECBEEMan::FindAppDisjCutNew(NetMan & net, AbcObjVect& topoNodes) {
//     const int TARGET_LEVEL = 16;
//     cout << "finding approximate disjoint cuts" << endl;
//     assert(disjCuts.empty());
//     assert(cutNtks.empty());

//     // init
//     cutNtks.resize(net.GetIdMaxPlus1());
//     disjCuts.resize(net.GetIdMaxPlus1());

//     // collect disjoint cuts and the corresponding cut networks
//     for (int targetNodeId = 0; targetNodeId < net.GetIdMaxPlus1(); ++targetNodeId) {
//         if (!net.IsNode(targetNodeId))
//             continue;
//         if (net.IsConst(targetNodeId))
//             continue;
//         // expand the frontier for TARGET_LEVEL times
//         BitVect visited(net.GetIdMaxPlus1(), 0);
//         IntVect frontier = {targetNodeId};
//         for (int level = 0; level < TARGET_LEVEL; ++level) {
//             // for each node in the frontier, expand its fanouts
//             IntVect newFrontier;
//             for (auto frontierNode: frontier) {
//                 if (net.IsObjPo(frontierNode))
//                     newFrontier.emplace_back(frontierNode);
//                 for (int iFanout = 0; iFanout < net.GetFanoutNum(frontierNode); ++iFanout) {
//                     auto pFanout = net.GetFanout(frontierNode, iFanout);
//                     if (!visited[pFanout->Id]) {
//                         visited[pFanout->Id] = 1;
//                         newFrontier.emplace_back(pFanout->Id);
//                     }
//                 }
//             }
//             frontier = newFrontier;
//         }
//         // the nodes in the frontier form the disjCut
//         for (auto frontierNode: frontier)
//             disjCuts[targetNodeId].emplace_back(net.GetObj(frontierNode));
//         // collect the cutNtk
//         for (auto pNode: topoNodes) {
//             if (visited[pNode->Id])
//                 cutNtks[targetNodeId].emplace_back(pNode);
//         }
//         for (ll j = 0; j < net.GetPoNum(); ++j) {
//             auto pPo = net.GetPo(j);
//             if (visited[pPo->Id])
//                 cutNtks[targetNodeId].emplace_back(pPo);
//         }
//     }
// }


void VECBEEMan::FindDisjCutOfNode(Abc_Obj_t* pObj, AbcObjList& disjCut) {
    disjCut.clear();
    ExpandCut(pObj, disjCut);
    Abc_Obj_t * pObjExpd = nullptr;
    while ((pObjExpd = ExpandWhich(disjCut)) != nullptr) {
        ExpandCut(pObjExpd, disjCut);
    }
}


void VECBEEMan::ExpandCut(Abc_Obj_t* pObj, AbcObjList& disjCut) {
    abc::Abc_Obj_t * pFanout = nullptr;
    ll i = 0;
    Abc_ObjForEachFanout(pObj, pFanout, i) {
        if (!abc::Abc_NodeIsTravIdCurrent(pFanout)) {
            if (abc::Abc_ObjFanoutNum(pFanout) || abc::Abc_ObjIsPo(pFanout)) {
                abc::Abc_NodeSetTravIdCurrent(pFanout);
                disjCut.emplace_back(pFanout);
            }
        }
    } 
}


Abc_Obj_t* VECBEEMan::ExpandWhich(AbcObjList& disjCut) {
    for (auto ppAbcObj1 = disjCut.begin(); ppAbcObj1 != disjCut.end(); ++ppAbcObj1) {
        auto ppAbcObj2 = ppAbcObj1;
        for (++ppAbcObj2; ppAbcObj2 != disjCut.end(); ++ppAbcObj2) {
            assert(poMarks[(*ppAbcObj1)->Id].size() == poMarks[(*ppAbcObj2)->Id].size());
            assert((*ppAbcObj1)->Id != (*ppAbcObj2)->Id);
            assert(topoIds[(*ppAbcObj1)->Id] != topoIds[(*ppAbcObj2)->Id]);
            auto isJoint = poMarks[(*ppAbcObj1)->Id] & poMarks[(*ppAbcObj2)->Id];
            if (isJoint.any()) {
                abc::Abc_Obj_t * pRet = nullptr;
                if (topoIds[(*ppAbcObj1)->Id] < topoIds[(*ppAbcObj2)->Id]) {
                    pRet = *ppAbcObj1;
                    disjCut.erase(ppAbcObj1);
                }
                else {
                    pRet = *ppAbcObj2;
                    disjCut.erase(ppAbcObj2);
                }
                return pRet;
            }
        }
    }
    return nullptr;
}


void VECBEEMan::CalcBoolDiffCut2Node(Simulator& appSmlt, AbcObjVect& topoNodes) {
    cout << "calculating boolean difference of cuts with regard to nodes" << endl;
    assert(topoNodes.size());
    assert(topoNodes[0]->pNtk == appSmlt.GetNet());
    timer::progress_display pd(topoNodes.size());
    bdCut2Nodes.resize(appSmlt.GetIdMaxPlus1());
    for (const auto & pObj: topoNodes) {
        ll i = appSmlt.GetId(pObj);
        if (!appSmlt.IsNode(pObj) || appSmlt.IsConst(pObj)) {
            ++pd;
            continue;
        }
        appSmlt.CalcLocBoolDiff(pObj, disjCuts[i], cutNtks[i], bdCut2Nodes[i]);
        ++pd;
    }
}


static void UpdSopForBoolAndPartDiff(Simulator& appSmlt, Abc_Obj_t* pObj, char* pSop, unordered_map<int, BitVect>& tempDat) {
    int nVars = abc::Abc_SopGetVarNum(pSop);
    int nFrame = appSmlt.GetFrameNumb();
    BitVect product(nFrame, 0);
    for (char * pCube = pSop; *pCube; pCube += nVars + 3) {
        bool isFirst = true;
        for (ll i = 0; pCube[i] != ' '; i++) {
            Abc_Obj_t * pFanin = Abc_ObjFanin(pObj, i);
            BitVect &datFi = tempDat.count(pFanin->Id)? tempDat[pFanin->Id]: *appSmlt.GetDat(pFanin->Id);
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
        assert(!isFirst);
        if (pCube == pSop)
            tempDat[pObj->Id] = product;
        else
            tempDat[pObj->Id] |= product;
    }
    // complement
    if (abc::Abc_SopIsComplement(pSop))
        tempDat[pObj->Id].flip();
}


static void CalcSomeLocBoolDiff(AbcObjVect& targetNodes, Simulator& appSmlt, vector<AbcObjList>& disjCuts, vector<AbcObjVect>& cutNtks, vector<vector<BitVect>>& bdCut2Nodes, timer::progress_display& pd, int start, int end) {
    assert(appSmlt.GetNetType() == NET_TYPE::SOP);
    for (int index = start; index < end; ++index) {
        auto pTarget = targetNodes[index];
        auto& cutNtk = cutNtks[pTarget->Id];
        auto& bdCut2Node = bdCut2Nodes[pTarget->Id];
        auto& disjCut = disjCuts[pTarget->Id];
        assert(pTarget->pNtk == appSmlt.GetNet());
        unordered_map<int, BitVect> tempDat;
        tempDat.clear();
        // flip the node
        tempDat[pTarget->Id] = ~(*appSmlt.GetDat(pTarget));
        // simulate
        for (auto& pInner: cutNtk) {
            assert(!Abc_ObjIsPi(pInner));
            assert(!Abc_NodeIsConst(pInner));
            if (Abc_ObjIsPo(pInner)) {
                assert(!abc::Abc_ObjIsComplement(pInner));
                Abc_Obj_t* pDriver = abc::Abc_ObjFanin0(pInner);
                if (tempDat.count(pDriver->Id))
                    tempDat[pInner->Id] = tempDat[pDriver->Id];
                else
                    tempDat[pInner->Id] = *appSmlt.GetDat(pDriver);
            }
            else
                UpdSopForBoolAndPartDiff(appSmlt, pInner, static_cast<char *>(pInner->pData), tempDat);
        }
        // get boolean difference from the node to its disjoint cuts
        bdCut2Node.resize(disjCut.size());
        int i = 0;
        for (auto pCut: disjCut) {
            bdCut2Node[i] = *appSmlt.GetDat(pCut) ^ tempDat[pCut->Id];
            ++i;
        }
        std::unique_lock<std::mutex> lock(mtx);
        ++pd;
        lock.unlock();
    }
}


void VECBEEMan::CalcBoolDiffCut2NodeParallelly(Simulator& appSmlt, AbcObjVect& topoNodes) {
    cout << "calculating boolean difference of cuts with regard to nodes" << endl;
    assert(topoNodes.size());
    assert(topoNodes[0]->pNtk == appSmlt.GetNet());

    // collect target nodes
    AbcObjVect targetNodes;
    targetNodes.reserve(topoNodes.size());
    for (const auto & pObj: topoNodes) {
        if (appSmlt.IsNode(pObj) && !appSmlt.IsConst(pObj))
            targetNodes.emplace_back(pObj);
    }

    // compute multiple thread parameters
    int targetNodeNum = targetNodes.size();
    int realThread = min(targetNodeNum, nThread);
    assert(realThread > 0);
    cout << "real thread number: " << realThread << endl;
    int chunkSize = targetNodeNum / realThread;
    int remainder = targetNodeNum % realThread;

    // multi-thread calculation
    bdCut2Nodes.resize(appSmlt.GetIdMaxPlus1());
    timer::progress_display pd(targetNodeNum);
    vector<thread> threads;
    int start = 0;
    for (int i = 0; i < realThread; ++i) {
        int end = start + chunkSize + (i < remainder? 1: 0);
        threads.emplace_back(CalcSomeLocBoolDiff, std::ref(targetNodes), std::ref(appSmlt), std::ref(disjCuts), std::ref(cutNtks), std::ref(bdCut2Nodes), std::ref(pd), start, end);
        start = end;
    }
    for (auto& thread: threads)
        thread.join();
}


void VECBEEMan::CalcBoolDiffPo2Node(Simulator& appSmlt, AbcObjVect& topoNodes) {
    cout << "calculating boolean difference of POs with regard to nodes" << endl;
    assert(topoNodes.size());
    assert(topoNodes[0]->pNtk == appSmlt.GetNet());
    ll nPo = appSmlt.GetPoNum();
    bdPo2Nodes.resize(nPo);
    // timer::progress_display pd(nPo);
    for (ll o = 0; o < nPo; ++o) {
        // init boolean difference
        auto & bdPo2Node = bdPo2Nodes[o];
        bdPo2Node.resize(appSmlt.GetIdMaxPlus1(), dynamic_bitset <ull> (nFrame, 0));
        // for each PO, update boolean difference
        for (ll n = 0; n < appSmlt.GetPoNum(); ++n) {
            auto pNodeN = appSmlt.GetPo(n);
            auto nId = appSmlt.GetId(pNodeN);
            if (n == o)
                bdPo2Node[nId].set(); 
            else
                bdPo2Node[nId].reset(); 
        }
        // for each node, update boolean difference
        for (auto it = topoNodes.rbegin(); it != topoNodes.rend(); ++it) {
            auto pNodeN = *it;
            if (!appSmlt.IsNode(pNodeN))
                continue;
            ll n = appSmlt.GetId(pNodeN);
            bdPo2Node[n].reset();
            ll i = 0;
            for (auto pCut: disjCuts[n]) {
                bdPo2Node[n] |= bdPo2Node[pCut->Id] & bdCut2Nodes[n][i];
                ++i;
            } 
        }
        // ++pd;
    }
}


static void CalcSomeBoolDiffPo2Node(Simulator & appSmlt, std::vector < std::vector < boost::dynamic_bitset <ull> > > & bdPo2Nodes, vector <Abc_Obj_t *> & topoNodes, std::vector < std::list <abc::Abc_Obj_t *> > & disjCuts, std::vector < std::vector < boost::dynamic_bitset <ull> > > & bdCut2Nodes, timer::progress_display & pd, ll start, ll end) {
    ll nFrame = appSmlt.GetFrameNumb();
    for (ll o = start; o < end; ++o) {
        auto & bdPo2Node = bdPo2Nodes[o];
        // init boolean difference
        bdPo2Node.resize(appSmlt.GetIdMaxPlus1(), dynamic_bitset <ull> (nFrame, 0));
        // for each PO, update boolean difference
        for (ll n = 0; n < appSmlt.GetPoNum(); ++n) {
            auto pNodeN = appSmlt.GetPo(n);
            auto nId = appSmlt.GetId(pNodeN);
            if (n == o)
                bdPo2Node[nId].set(); 
            else
                bdPo2Node[nId].reset(); 
        }
        // for each node, update boolean difference
        for (auto it = topoNodes.rbegin(); it != topoNodes.rend(); ++it) {
            auto pNodeN = *it;
            if (!appSmlt.IsNode(pNodeN))
                continue;
            ll n = appSmlt.GetId(pNodeN);
            bdPo2Node[n].reset();
            ll i = 0;
            assert(disjCuts[n].size() == bdCut2Nodes[n].size());
            for (auto pCut: disjCuts[n]) {
                bdPo2Node[n] |= bdPo2Node[pCut->Id] & bdCut2Nodes[n][i];
                ++i;
            } 
        }
        std::unique_lock<std::mutex> lock(mtx);
        ++pd;
        lock.unlock();
    }
}


void VECBEEMan::CalcBoolDiffPo2NodeParallelly(Simulator& appSmlt, AbcObjVect& topoNodes) {
    cout << "calculating boolean difference of POs with regard to nodes" << endl;
    assert(topoNodes.size());
    assert(topoNodes[0]->pNtk == appSmlt.GetNet());
    int nPo = appSmlt.GetPoNum();
    bdPo2Nodes.resize(nPo);

    int realThread = min(nPo, nThread);
    assert(realThread > 0);
    cout << "real thread number: " << realThread << endl;
    ll chunkSize = nPo / realThread;
    ll remainder = nPo % realThread;

    timer::progress_display pd(nPo);
    vector<thread> threads;
    ll start = 0;
    for (ll i = 0; i < realThread; ++i) {
        ll end = start + chunkSize + (i < remainder? 1: 0);
        threads.emplace_back(CalcSomeBoolDiffPo2Node, std::ref(appSmlt), std::ref(bdPo2Nodes), std::ref(topoNodes), std::ref(disjCuts), std::ref(bdCut2Nodes), std::ref(pd), start, end);
        start = end;
    }
    for (auto& thread: threads)
        thread.join();
}


static void GetNewValueForBlock(Simulator & smlt, const std::vector <ll> & faninIds, const std::string & sop, ull & value, ll iBlock) {
    if (sop == " 0\n") {
        value = 0;
        return;
    }
    if (sop == " 1\n") {
        value = numeric_limits <ull>::max();
        return;
    }
    char * pSop = const_cast <char *> (sop.c_str());
    ll nVars = Abc_SopGetVarNum(pSop);
    assert(nVars == faninIds.size());

    ull product = 0;
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
                        product = ~GetBlockFromDynBitset(*smlt.GetDat(faninId), iBlock);
                    }
                    else
                        product &= ~GetBlockFromDynBitset(*smlt.GetDat(faninId), iBlock);
                    break;
                case '1':
                    if (isFirst) {
                        isFirst = false;
                        product = GetBlockFromDynBitset(*smlt.GetDat(faninId), iBlock);
                    }
                    else
                        product &= GetBlockFromDynBitset(*smlt.GetDat(faninId), iBlock);
                    break;
                default:
                    assert(0);
            }
        }
        if (isFirst) {
            isFirst = false;
            // product.set();
            product = numeric_limits <ull>::max();
        }
        assert(!isFirst);
        if (pCube == pSop)
            value = product;
        else
            value |= product;
    }

    // complement
    if (abc::Abc_SopIsComplement(pSop))
        // value.flip();
        value = ~value;
}


static BigInt GetValue(vector<BitVect> & dat, int iPatt, bool isSign, int msb) {
    int lsb = 0;
    int shift = msb - lsb;
    assert(msb < 200);
    BigInt ret(0);
    for (ll k = msb; k >= lsb; --k) {
        ret <<= 1;
        if (dat[k][iPatt])
            ++ret;
    }
    if (isSign && ret >= (BigInt(1) << shift))
        ret = -((BigInt(1) << (shift + 1)) - ret);
    return ret;
}


static void CalcSomeLACNonERErrsNew(Simulator& appSmlt, LACMan& lacMan, vector<vector<BitVect> >& bdPo2Nodes, vector<BigInt>& YAcc, BigInt& runMin, timer::progress_display& pd, int startIndex, int endIndex, METR_TYPE metrType, bool& earlyStop, BigInt& backErrInt) {
    int nPo = appSmlt.GetPoNum();
    int nFrame = appSmlt.GetFrameNumb();

    for (int lacId = startIndex; lacId < endIndex; ++lacId) {
        // get LAC
        auto pLac = lacMan.GetResubLac(lacId);
        int targId = pLac->GetTargId();

        // get local boolean difference
        BitVect newValue(nFrame, 0);
        GetNewValue(appSmlt, pLac->GetDivIds(), pLac->GetSop(), newValue);
        auto isChanged = (*appSmlt.GetDat(targId)) ^ newValue;

        // calculate error
        BigInt ser = 0;
        if (isChanged.none())
            ser = backErrInt;
        else {
            // get temp outputs
            vector<BitVect> tempOutps(nPo);
            for (int j = 0; j < nPo; ++j) {
                auto poId = appSmlt.GetPoId(j);
                tempOutps[j] = *appSmlt.GetDat(poId) ^ (isChanged & bdPo2Nodes[j][targId]); 
            }
            // get error
            if (metrType == METR_TYPE::MED) {
                for (int iPatt = 0; iPatt < nFrame; ++iPatt) {
                    auto YNew = GetValue(tempOutps, iPatt, false, nPo - 1);
                    ser += abs(YNew - YAcc[iPatt]);
                    if (ser > runMin)
                        break;
                }
            }
            else if (metrType == METR_TYPE::MSE) {
                for (int iPatt = 0; iPatt < nFrame; ++iPatt) {
                    auto YNew = GetValue(tempOutps, iPatt, false, nPo - 1);
                    ser += (YNew - YAcc[iPatt]) * (YNew - YAcc[iPatt]);
                    if (ser > runMin)
                        break;
                }
            }
            else
                assert(0);
        }
        
        // // early stop
        // if (earlyStop)
        //     return;
        // if (ser <= backErrInt) {
        //     std::unique_lock<std::mutex> lock(mtx);
        //     earlyStop = true;
        //     lock.unlock();
        //     cout << "early stop" << endl;
        //     pLac->SetErrPro(ser);
        //     return;
        // }

        // running min
        std::unique_lock<std::mutex> lock(mtx);
        runMin = min(runMin, ser);
        ++pd;
        lock.unlock();

        // set error
        pLac->SetErrPro(ser);
    }
}


void VECBEEMan::CalcLACNonERErrsNew(Simulator& accSmlt, Simulator& appSmlt, LACMan& lacMan, const BigInt& uppBound, BigInt& backErrInt) {
    cout << "calculating LAC errors" << endl;

    assert(IsPIOSame(accSmlt, appSmlt));
    assert((nFrame & 63) == 0);
    assert(uppBound >= 0);
    assert(lacType == LAC_TYPE::RESUB);
    assert(!isSign);
    assert(appSmlt.GetPoNum() < 200);

    vector<BigInt> YAcc(nFrame, 0);
    for (ll iPatt = 0; iPatt < nFrame; ++iPatt)
        YAcc[iPatt] = accSmlt.GetOutpPro(iPatt, isSign);

    int lacNum = lacMan.GetLacNum();
    int realThread = min(nThread, lacNum);
    assert(realThread > 0);
    cout << "real thread number: " << realThread << endl;
    int chunkSize = lacNum / realThread;
    int remainder = lacNum % realThread;

    timer::progress_display pd(lacMan.GetLacNum());
    BigInt runMin = uppBound + 1;
    vector<thread> threads;
    bool earlyStop = false;
    int start = 0;
    for (int i = 0; i < realThread; ++i) {
        int end = start + chunkSize + (i < remainder? 1: 0);
        threads.emplace_back(CalcSomeLACNonERErrsNew, std::ref(appSmlt), std::ref(lacMan), std::ref(bdPo2Nodes), std::ref(YAcc), std::ref(runMin), std::ref(pd), start, end, metrType, std::ref(earlyStop), std::ref(backErrInt));
        start = end;
    }
    for (auto& thread: threads)
        thread.join();
}


static void CalcSomeLACERErrs(Simulator & appSmlt, LACMan & lacMan, std::vector < std::vector < boost::dynamic_bitset <ull> > > & bdPo2Nodes, timer::progress_display & pd, int startIndex, int endIndex, METR_TYPE metrType, std::vector <boost::dynamic_bitset<ull>> & IsCurrentPODifferent) {
    int nPo = appSmlt.GetPoNum();
    int nFrame = appSmlt.GetFrameNumb();

    for (int lacId = startIndex; lacId < endIndex; ++lacId) {
        auto pLac = lacMan.GetLac(lacId);
        int targId = pLac->GetTargId();

        // calculate $\partial n / \partial LAC$
        boost::dynamic_bitset <ull> isChanged(nFrame, 0);
        auto & specLac = *dynamic_pointer_cast <ResubLAC>(pLac);
        auto divIds = specLac.GetDivIds(); 
        auto sop = specLac.GetSop();
        dynamic_bitset <ull> newValue(nFrame, 0);
        GetNewValue(appSmlt, divIds, sop, newValue);
        isChanged = (*appSmlt.GetDat(targId)) ^ newValue;

        // calculate error
        BigInt err = 0;
        BitVect diff(nFrame, 0);
        if (metrType == METR_TYPE::ER) {
            for (int j = 0; j < nPo; ++j)
                diff |= (IsCurrentPODifferent[j] ^ (isChanged & bdPo2Nodes[j][targId]));
            err = diff.count();
        }
        else if (metrType == METR_TYPE::MHD) {
            for (int j = 0; j < nPo; ++j) {
                diff = (IsCurrentPODifferent[j] ^ (isChanged & bdPo2Nodes[j][targId]));
                err += diff.count();
            }
        }
        pLac->SetErrPro(err);

        // update progress
        std::unique_lock<std::mutex> lock(mtx);
        ++pd;
        lock.unlock();
    }
}


void VECBEEMan::CalcLACERErrs(Simulator & accSmlt, Simulator & appSmlt, LACMan & lacMan) {
    assert(IsPIOSame(accSmlt, appSmlt));
    assert((nFrame & 63) == 0);
    ll nPo = appSmlt.GetPoNum();
    assert(nPo < 200);
    assert(metrType == METR_TYPE::ER || metrType == METR_TYPE::MHD);

    cout << "calculating LAC errors" << endl;

    // compute whether the POs in the current circuit is correct or not
    vector <dynamic_bitset<ull>> IsCurrentPODifferent(nPo);
    for (ll j = 0; j < nPo; ++j) {
        ll poIdAcc = accSmlt.GetPoId(j);
        ll poIdApp = appSmlt.GetPoId(j);
        IsCurrentPODifferent[j] = *accSmlt.GetDat(poIdAcc) ^ *appSmlt.GetDat(poIdApp); 
    }

    // compute in parallel
    int lacNum = lacMan.GetLacNum();
    int realThread = min(nThread, lacNum);
    assert(realThread > 0);
    cout << "real thread number: " << realThread << endl;
    ll chunkSize = lacNum / realThread;
    ll remainder = lacNum % realThread;
    timer::progress_display pd(lacMan.GetLacNum());
    vector<thread> threads;
    ll start = 0;
    for (ll i = 0; i < realThread; ++i) {
        ll end = start + chunkSize + (i < remainder? 1: 0);
        threads.emplace_back(CalcSomeLACERErrs, std::ref(appSmlt), std::ref(lacMan), std::ref(bdPo2Nodes), std::ref(pd), start, end, metrType, std::ref(IsCurrentPODifferent));
        start = end;
    }
    for (auto& thread: threads)
        thread.join();
}

