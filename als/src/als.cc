#include "als.h"


using namespace std;
using namespace abc;
using namespace boost;


void ALSOpt::Print() {
    cout << endl << "******************** options ********************" << endl;
    cout << (isSign? "use signed output": "use unsigned output") << endl;
    cout << "enable fast error estimation = " << (enableFastErrEst? "YES": "NO") << endl;
    cout << "source seed = " << sourceSeed << endl;
    cout << "lacType = " << lacType << endl;
    cout << "distrType = " << distrType << endl;
    cout << "metrType = " << metrType << endl;
    cout << "#simulation patterns = " << nFrame << endl;
    cout << "#simulation patterns for candidate AppResub generation = " << nFrame4ResubGen << endl;
    cout << "max #candidate AppResubs = " << maxCandResub << endl;
    cout << "#threads = " << nThread << endl;
    // cout << "max level difference = " << maxLevelDiff << endl;
    cout << metrType << " upper bound = " << errUppBound << endl;
    if (metrType == METR_TYPE::MED) {
        AbcMan abcMan;
        int nPo = Abc_NtkPoNum(abcMan.GetNet());
        cout << "NMED upper bound = " << BigFlt(errUppBound) / BigFlt((BigInt(1) << nPo) - 1) << endl;
    }
    cout << "standard cell path: " << standCellPath << endl;
    cout << "output path: " << outpPath << endl;
    cout << "*************************************************" << endl << endl;
}


ALSMan::ALSMan(ALSOpt & opt): isSign(opt.isSign), enableFastErrEst(opt.enableFastErrEst), sourceSeed(opt.sourceSeed), seed(0), lacType(opt.lacType), distrType(opt.distrType), metrType(opt.metrType), nFrame(opt.nFrame), nFrame4ResubGen(opt.nFrame4ResubGen), maxCandResub(opt.maxCandResub), nThread(opt.nThread), maxLevelDiff(INT_MAX), round(0), errUppBound(opt.errUppBound), maxDelay(DBL_MAX), accNet(NetMan(opt.pNtk, true)), standCellPath(opt.standCellPath), outpPath(opt.outpPath) {
    if (accNet.GetNetType() == NET_TYPE::GATE) {
        cout << "convert gate netlist into AIG" << endl;
        accNet.Comm("st; compress2rs; logic; sop; ps;");
    }
    randGen.seed(sourceSeed);
}


void ALSMan::RunMultipleSelection() {
    assert(lacType == LAC_TYPE::RESUB);

    // initialize
    cout << endl << "******************** init status ********************" << endl;
    assert(maxDelay == DBL_MAX);
    maxDelay = Eval(accNet, 0.0, true, true);
    cout << "max delay = " << maxDelay << endl;
    accNet.Comm("st; logic; sop; ps;");
    auto currNet = accNet;
    cout << "*****************************************************" << endl << endl;

    // multiple selection
    auto start = chrono::system_clock::now();
    double err = 0.0;
    while (true) {
        cout << "------------------- round " << ++round << "-------------------" << endl;
        err = ApplyMultipleLACs(currNet);
        auto duration = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start);
        cout << "actual runtime = " << double(duration.count()) * chrono::microseconds::period::num / chrono::microseconds::period::den << " sec" << endl;
        cout << "------------------------------------------------" << endl << endl;
        if (DoubleGreat(err, errUppBound))
            break;
    }

    // single selection
    while (true) {
        cout << "------------------- round " << ++round << "-------------------" << endl;
        err = ApplyTheBestLAC(currNet);
        auto duration = chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start);
        cout << "actual runtime = " << double(duration.count()) * chrono::microseconds::period::num / chrono::microseconds::period::den << " sec" << endl;
        cout << "------------------------------------------------" << endl << endl;
        if (DoubleGreat(err, errUppBound))
            break;
    }
}


double ALSMan::ApplyTheBestLAC(NetMan & net) {
    // backup network
    auto backNet = net;

    // use new seed
    seed = NewSeed();
    auto backErr = CalcErr(accNet, net, isSign, seed, nFrame, metrType, distrType);
    cout << "base " << metrType << " = " << backErr << endl;
    if (DoubleGreat(backErr, errUppBound)) {
        Eval(backNet, backErr, true);
        cout << "WARNING: exceed error bound due to unstable " << metrType << " measurement" << endl;
        return DBL_MAX;
    }
    
    // create constant nodes
    net.CreateConst(true);

    // get target nodes
    auto nodes = net.TopoSortWithIds();

    // generate LACs
    LACMan lacMan(lacType, nThread);
    lacMan.Gen012ResubLACsPro(net, nodes, seed, maxLevelDiff, nFrame4ResubGen, maxCandResub);
    cout << "#lacs = " << lacMan.GetLacNum() << endl;

    // error estimation
    const BigInt uppBound = (BigInt)(BigFlt(nFrame) * BigFlt(errUppBound)) + 8192;
    BigInt backErrInt = (BigInt)(BigFlt(nFrame) * BigFlt(backErr)) - 1;
    {// this block is for memory management
    VECBEEMan vecbeeMan(isSign, seed, nFrame, metrType, lacType, distrType, nThread);
    vecbeeMan.BatchErrEstPro(accNet, net, lacMan, uppBound, enableFastErrEst, backErrInt);
    }

    // apply best LAC
    assert(lacType == LAC_TYPE::RESUB);
    auto pBestLac = lacMan.GetBestResubLacConsiderRealSize(net);
    ApplyLacPro(net, pBestLac, backErr);

    // get error
    double err = CalcErr(accNet, net, isSign, seed, nFrame, metrType, distrType);
    cout << "current " << metrType << " = " << err << endl;
    if (metrType == METR_TYPE::MED)
        cout << "NMED = " << (BigFlt)(err) / BigFlt((BigInt(1) << net.GetPoNum()) - 1) << endl;

    // check error estimation
    if (DoubleLess(err, errUppBound)) {
        auto estErr = double(BigFlt(pBestLac->GetErrPro()) / BigFlt(nFrame));
        if (!enableFastErrEst && !DoubleEqual(estErr, err, 1e-4)) {
            cout << "================== WARNING: wrong error estimation ================" << endl;
            cout << setprecision(10) << estErr << "\t" << err << endl;
            assert(0);
        }
    }

    // simplify without errors
    const ll SYNTH_ROUND = 10;
    if (round % SYNTH_ROUND == 0 || backErr > 0.5 * errUppBound)
        ExactSimpl(net);
    else {
        net.CleanUp();
        net.MergeConst();
        net.PrintStat();
    }

    // measure, synthesis & mapping, output
    const int EVAL_ROUND = 10;
    if (DoubleLessEqual(err, errUppBound)) {
        if (round % EVAL_ROUND == 0 || err > 0.5 * errUppBound)
            Eval(net, err);
    }
    else {
        net = backNet;
        Eval(backNet, backErr, true);
    }

    // return error
    return err;
}


double ALSMan::ApplyMultipleLACs(NetMan & net) {
    // check, change seed, & backup network
    assert(lacType == LAC_TYPE::RESUB);
    assert(distrType == DISTR_TYPE::UNIF);
    assert(isSign == false);
    seed = NewSeed();
    auto backNet = net;
    
    // create constant nodes & get target nodes in topological order
    net.CreateConst(true);
    auto targetNodes = net.TopoSortWithIds();

    // generate LACs
    LACMan lacMan(lacType, nThread);
    lacMan.Gen012ResubLACsPro(net, targetNodes, seed, maxLevelDiff, nFrame4ResubGen, maxCandResub);
    cout << "#generated lacs = " << lacMan.GetLacNum() << endl;

    // build error distance miter
    IntVect miterId2AppId, appId2MiterId;
    NetManPtr pMiterNet = nullptr;
    if (metrType == METR_TYPE::ER)
        pMiterNet = BuildErrorRateMiter(accNet, net, miterId2AppId, appId2MiterId);
    else if (metrType == METR_TYPE::MED)
        pMiterNet = BuildErrorDistanceMiter(accNet, net, miterId2AppId, appId2MiterId);
    else
        pMiterNet = BuildMiterWithYosys(accNet, net, miterId2AppId, appId2MiterId);

    // simulate
    Simulator miterSmlt(*pMiterNet, seed, nFrame);
    miterSmlt.InpUnifFast();
    miterSmlt.Sim();
    auto backErr = miterSmlt.GetError();
    cout << "base " << metrType << " = " << backErr << endl;
    if (DoubleGreat(backErr, errUppBound)) {
        Eval(backNet, backErr, true);
        cout << "WARNING: exceed error bound due to unstable " << metrType << " measurement" << endl;
        return DBL_MAX;
    }

    // error estimation
    BigInt errUppBoundInt = (BigInt)(BigFlt(nFrame) * BigFlt(errUppBound)) + 1;
    {// this block is for memory management
    VECBEEMan vecbeeMan(isSign, seed, nFrame, metrType, lacType, distrType, nThread);
    // vecbeeMan.EstimateErrorBoundForEachLAC(accSmlt, appSmlt, lacMan, errUppBoundInt, enableFastErrEst);
    vecbeeMan.EstimateErrorBoundForEachLAC(miterSmlt, lacMan, errUppBoundInt, enableFastErrEst, miterId2AppId, appId2MiterId);
    }

    // collect promising LACs for each node
    BigInt errorMargin = (BigInt)(BigFlt(errUppBound - backErr) * BigFlt(nFrame));
    cout << "error margin = " << errorMargin << endl;
    assert(errorMargin >= 0);
    int nPromisingLacs = lacMan.CollectPromisingLACs(net.GetIdMaxPlus1(), errorMargin);
    if (nPromisingLacs == 0) {
        cout << "no promising LACs" << endl;
        return DBL_MAX;
    }

    // find promising set of LACs, $\Pi$
    LACPtrVect goodLACs;
    FindGoodLACsSmallMemory(net, targetNodes, lacMan, errorMargin, goodLACs);

    // apply best LAC & traditional logic synthesis
    // ApplyLacsConsideringErrors(net, goodLACs, accSmlt);
    ApplyLacs(net, goodLACs);

    // get error
    // auto err = ComputeError(accSmlt, net);
    double err = CalcErr(accNet, net, isSign, seed, nFrame, metrType, distrType);
    cout << "current " << metrType << " = " << err << endl;
    if (metrType == METR_TYPE::MED)
        cout << "NMED = " << (BigFlt)(err) / BigFlt((BigInt(1) << net.GetPoNum()) - 1) << endl;

    // traditional logic synthesis
    ExactSimpl(net);

    // measure, synthesis & mapping, output
    if (DoubleLessEqual(err, errUppBound))
        Eval(net, err);
    else {
        net = backNet;
        Eval(backNet, backErr, true);
    }

    // return error
    return err;
}


void ALSMan::FindGoodLACs(NetMan& net, const IntVect& targetNodes, LACMan& lacMan, BigInt& errorMargin, RETURN_VAR LACPtrVect& goodLACs) {
    assert(net.GetIdMaxPlus1() == lacMan.GetPromisingLacSize());

    // nodeIds[i] is the node id of the ith item in the network
    // item id starts from 1
    IntVect nodeIds;
    nodeIds.reserve(targetNodes.size());
    nodeIds.emplace_back(-1);
    for (auto id: targetNodes) {
        if (lacMan.GetPromisingLacNum(id))
            nodeIds.emplace_back(id);
    }

    // formulate as a knapsack problem
    // M is the number of items
    // W is the capacity of the knapsack
    // for each item i, there are choice[i] choices
    // weight[i][j] is the weight of the jth choice of item i
    // value[i][j] is the value of the jth choice of item i
    // item id is from 1 to M
    const int M = nodeIds.size() - 1;
    const BigInt W = errorMargin;
    IntVect choice(M + 1, 0);
    BigInt2DVect weight(M + 1);
    LL2DVect value(M + 1);
    int count = 0;
    for (int i = 1; i <= M; ++i) {
        int id = nodeIds[i];
        choice[i] = lacMan.GetPromisingLacNum(id);
        weight[i].reserve(choice[i]);
        value[i].reserve(choice[i]);
        for (int j = 0; j < choice[i]; ++j) {
            auto pLac = lacMan.GetPromisingResubLac(id, j);
            auto error = pLac->GetErrPro();
            // assert(error <= LL_WARN);
            weight[i].emplace_back(error);
            value[i].emplace_back(pLac->GetSizeGain());
        }
        count += choice[i];
    }
    // cout << "total number of choices = " << count << endl;

    // print the knapsack problem
    // cout << "M = " << M << ", W = " << W << endl;
    // for (int i = 1; i <= M; ++i) {
    //     cout << "item " << i << " (node " << nodeIds[i] << ")" << " has " << choice[i] << " choices: " << endl;
    //     for (int j = 0; j < choice[i]; ++j)
    //         cout << "choice " << j << ": weight = " << weight[i][j] << ", value = " << value[i][j] << endl;
    // }
    
    // try another way to find good LACs: solve the dual problem
    // get an upper bound of the maximum value
    int V = 0;
    for (int i = 1; i <= M; ++i) {
        ll maxVal = 0;
        for (int j = 0; j < choice[i]; ++j)
            maxVal = max(maxVal, value[i][j]);
        V += maxVal;
    }
    cout << "V = " << V << endl;
    // dp[i][v] is the minimum weight of the first i items with value v
    BigInt2DVect dpDual(M + 1, BigIntVect(V + 1, MAX_BIGINT));
    // chooseDual[i][v] is the choice of the ith item with value v, -1 means not choosing
    Int2DVect chooseDual(M + 1, IntVect(V + 1, -1));
    // dynamic programming
    for (int i = 1; i <= M; ++i) {
        for (int v = 0; v <= V; ++v) {
            dpDual[i][v] = dpDual[i - 1][v];
            chooseDual[i][v] = -1;
            for (int j = 0; j < choice[i]; ++j) {
                if (v > value[i][j]) {
                    auto source = dpDual[i - 1][v - value[i][j]];
                    auto newDp = (source == MAX_BIGINT)? source: source + weight[i][j];
                    if (dpDual[i][v] > newDp) {
                        // assert(newDp <= LL_WARN);
                        dpDual[i][v] = newDp;
                        chooseDual[i][v] = j;
                    }
                }
                else if (v == value[i][j]) {
                    if (dpDual[i][v] > weight[i][j]) {
                        dpDual[i][v] = weight[i][j];
                        chooseDual[i][v] = j;
                    }
                }
            }
        }
    }
    // // print dpDual
    // for (int v = 0; v <= V; ++v)
    //     cout << "dpDual[M][" << v << "] = " << dpDual[M][v] << endl; 
    // get the maximum value under the constraint weight <= W
    ll maxValDual = 0;
    for (int v = 1; v <= V; ++v) {
        if (dpDual[M][v] <= W)
            maxValDual = v;
    }
    cout << "maximum value = " << maxValDual << ", weight = " << dpDual[M][maxValDual] << ", estimated delta error = " << (BigFlt)(dpDual[M][maxValDual]) / (BigFlt)(nFrame) << endl;
    assert(dpDual[M][maxValDual] != MAX_BIGINT);
    // get the selected items and their choices
    ll v = maxValDual;
    for (int i = M; i >= 1; --i) {
        if (chooseDual[i][v] != -1) {
            // cout << "item " << i << " (node " << nodeIds[i] << ")" << " is chosen, choice = " << chooseDual[i][v] << ", weight = " << weight[i][chooseDual[i][v]] << ", value = " << value[i][chooseDual[i][v]] << endl;
            goodLACs.emplace_back(lacMan.GetPromisingLac(nodeIds[i], chooseDual[i][v]));
            v -= value[i][chooseDual[i][v]];
        }
    }
    // print good LACs
    // cout << "good LACs: " << endl;
    // for (auto pLac: goodLACs)
    //     dynamic_pointer_cast<ResubLAC>(pLac)->Print();
}


void ALSMan::FindGoodLACsSmallMemory(NetMan& net, const IntVect& targetNodes, LACMan& lacMan, BigInt& errorMargin, RETURN_VAR LACPtrVect& goodLACs) {
    assert(net.GetIdMaxPlus1() == lacMan.GetPromisingLacSize());

    // nodeIds[i] is the node id of the ith item in the network
    // item id starts from 1
    IntVect nodeIds;
    nodeIds.reserve(targetNodes.size());
    nodeIds.emplace_back(-1);
    for (auto id: targetNodes) {
        if (lacMan.GetPromisingLacNum(id))
            nodeIds.emplace_back(id);
    }

    // formulate as a knapsack problem
    // M is the number of items
    // W is the capacity of the knapsack
    // for each item i, there are choice[i] choices
    // weight[i][j] is the weight of the jth choice of item i
    // value[i][j] is the value of the jth choice of item i
    // item id is from 1 to M
    const int M = nodeIds.size() - 1;
    const BigInt W = errorMargin;
    IntVect choice(M + 1, 0);
    BigInt2DVect weight(M + 1);
    LL2DVect value(M + 1);
    int count = 0;
    for (int i = 1; i <= M; ++i) {
        int id = nodeIds[i];
        choice[i] = lacMan.GetPromisingLacNum(id);
        weight[i].reserve(choice[i]);
        value[i].reserve(choice[i]);
        for (int j = 0; j < choice[i]; ++j) {
            auto pLac = lacMan.GetPromisingResubLac(id, j);
            auto error = pLac->GetErrPro();
            // assert(error <= LL_WARN);
            weight[i].emplace_back(error);
            value[i].emplace_back(pLac->GetSizeGain());
        }
        count += choice[i];
    }
    // cout << "total number of choices = " << count << endl;

    // print the knapsack problem
    // cout << "M = " << M << ", W = " << W << endl;
    // for (int i = 1; i <= M; ++i) {
    //     cout << "item " << i << " (node " << nodeIds[i] << ")" << " has " << choice[i] << " choices: " << endl;
    //     for (int j = 0; j < choice[i]; ++j)
    //         cout << "choice " << j << ": weight = " << weight[i][j] << ", value = " << value[i][j] << endl;
    // }
    
    // try another way to find good LACs: solve the dual problem
    // get an upper bound of the maximum value
    int V = 0;
    for (int i = 1; i <= M; ++i) {
        ll maxVal = 0;
        for (int j = 0; j < choice[i]; ++j)
            maxVal = max(maxVal, value[i][j]);
        V += maxVal;
    }
    cout << "V = " << V << endl;
    // dpDual[i][v] is the minimum weight of the first i items with value v
    BigIntVect dpDual(V + 1, MAX_BIGINT);
    // chooseDual[i][v] is the choice of the ith item with value v, -1 means not choosing
    Int2DVect chooseDual(M + 1, IntVect(V + 1, -1));
    // dynamic programming
    cout << "solve the dual problem..." << endl;
    timer::progress_display pd(M);
    for (int i = 1; i <= M; ++i) {
        for (int v = V; v >= 0; --v) {
            // dpDual[v] = dpDual[v];
            chooseDual[i][v] = -1;
            for (int j = 0; j < choice[i]; ++j) {
                if (v > value[i][j]) {
                    const auto& source = dpDual[v - value[i][j]];
                    auto newDp = (source == MAX_BIGINT)? source: source + weight[i][j];
                    if (dpDual[v] > newDp) {
                        // assert(newDp <= LL_WARN);
                        dpDual[v] = newDp;
                        chooseDual[i][v] = j;
                    }
                }
                else if (v == value[i][j]) {
                    if (dpDual[v] > weight[i][j]) {
                        dpDual[v] = weight[i][j];
                        chooseDual[i][v] = j;
                    }
                }
            }
        }
        ++pd;
    }
    // // print dpDual
    // for (int v = 0; v <= V; ++v)
    //     cout << "dpDual[M][" << v << "] = " << dpDual[M][v] << endl; 
    // get the maximum value under the constraint weight <= W
    ll maxValDual = 0;
    for (int v = 1; v <= V; ++v) {
        if (dpDual[v] <= W)
            maxValDual = v;
    }
    cout << "maximum value = " << maxValDual << ", weight = " << dpDual[maxValDual] << ", estimated delta error = " << (BigFlt)(dpDual[maxValDual]) / (BigFlt)(nFrame) << endl;
    assert(dpDual[maxValDual] != MAX_BIGINT);
    // get the selected items and their choices
    ll v = maxValDual;
    for (int i = M; i >= 1; --i) {
        if (chooseDual[i][v] != -1) {
            // cout << "item " << i << " (node " << nodeIds[i] << ")" << " is chosen, choice = " << chooseDual[i][v] << ", weight = " << weight[i][chooseDual[i][v]] << ", value = " << value[i][chooseDual[i][v]] << endl;
            goodLACs.emplace_back(lacMan.GetPromisingLac(nodeIds[i], chooseDual[i][v]));
            v -= value[i][chooseDual[i][v]];
        }
    }
    // print good LACs
    // cout << "good LACs: " << endl;
    // for (auto pLac: goodLACs)
    //     dynamic_pointer_cast<ResubLAC>(pLac)->Print();
}


unsigned ALSMan::NewSeed() {
    boost::uniform_int <> unDistr(numeric_limits <int>::min(), numeric_limits <int>::max());
    unsigned _seed = static_cast <unsigned> (unDistr(randGen));
    cout << "new seed = " << _seed << endl;
    return _seed;
}


void ALSMan::ApplyLacPro(NetMan & net, std::shared_ptr <LAC> pLac, double backErr) {
    if (lacType == LAC_TYPE::RESUB) {
        assert(net.GetNetType() == NET_TYPE::SOP);
        net.GetLev();
        auto pSpecLac = dynamic_pointer_cast <ResubLAC> (pLac);
        auto targId = pSpecLac->GetTargId();
        auto faninIds = pSpecLac->GetDivIds();
        auto sop = pSpecLac->GetSop();
        cout << "replace " << net.GetObj(targId);
        cout << "(l=" << net.GetObjLev(targId) << ") with old fanins (";
        for (ll i = 0; i < net.GetFaninNum(targId); ++i)
            cout << net.GetFaninId(targId, i) << "(l=" << net.GetObjLev(net.GetFanin(targId, i)) << "),";
        cout << ")";
        cout << " by ";
        cout << "(";
        for (const auto & faninId: faninIds)
            cout << faninId << "(l=" << net.GetObjLev(faninId) << "),";
        cout << ")";
        cout << " with estimated error " << double(BigFlt(pSpecLac->GetErrPro()) / BigFlt(nFrame)) << endl;
        cout << " using function:" << endl;
        cout << sop;

        auto consts = net.CreateConst();
        if (sop == " 0\n") {
            net.Replace(targId, consts.first);
            net.PropagateConst(consts.first);
        }
        else if (sop == " 1\n") {
            net.Replace(targId, consts.second);
            net.PropagateConst(consts.second);
        }
        else if (sop == "1 1\n") {
            assert(faninIds.size() == 1);
            net.Replace(targId, faninIds[0]);
        }
        else if (sop == "0 1\n") {
            assert(faninIds.size() == 1);
            net.ReplaceByComplementedObj(targId, faninIds[0]);
        }
        else {
            auto newNodeId = net.CreateNodeAIG(faninIds, sop);
            net.Replace(targId, newNodeId);
        }
    }
    else
        assert(0);
}


void ALSMan::ApplyLacs(NetMan& net, LACPtrVect& lacs) {
    assert(lacs.size());
    assert(lacType == LAC_TYPE::RESUB);
    assert(net.GetNetType() == NET_TYPE::SOP);

    // transfer fanouts
    cout << "apply " << lacs.size() << " LACs" << endl;
    for (const auto & pLac: lacs) {
        auto pResubLac = dynamic_pointer_cast<ResubLAC> (pLac);
        int targId = pResubLac->GetTargId();
        auto divIds = pResubLac->GetDivIds();
        auto sop = pResubLac->GetSop();
        // pResubLac->Print();
        int newNodeId = net.CreateNode(divIds, sop);
        net.TransfFanout(targId, newNodeId);
    }

    // remove redundant nodes
    net.CleanUp();
    net.Check();
}


// void ALSMan::ApplyLacsConsideringErrors(NetMan& net, LACPtrVect& lacs, Simulator& accSmlt) {
//     assert(lacs.size());
//     assert(lacType == LAC_TYPE::RESUB);
//     assert(net.GetNetType() == NET_TYPE::SOP);

//     // sort LACs by error (ascending)
//     sort(lacs.begin(), lacs.end(), [](const LACPtr& pLac1, const LACPtr& pLac2) {
//         auto pResubLac1 = dynamic_pointer_cast<ResubLAC> (pLac1);
//         auto pResubLac2 = dynamic_pointer_cast<ResubLAC> (pLac2);
//         return pResubLac1->GetErrPro() < pResubLac2->GetErrPro();
//     });
//     cout << "lacs in pool: " << lacs.size() << endl;

//     // try applying LACs
//     int newLacNum = lacs.size();
//     while (true) {
//         cout << "try applying " << newLacNum << " LACs, ";
//         lacs.resize(newLacNum);
//         assert(newLacNum >= 1);
//         auto replInfo = TempApplyLacs(net, lacs, lacType, false);
//         auto err = ComputeError(accSmlt, net); 
//         cout << "causing " << metrType << " = " << err << endl;
//         // if there is <=1 LAC, break
//         if (newLacNum <= 1)
//             break;
//         // here, #LACs >= 2
//         bool stable = CheckError(net, err, 666);
//         // if stable and error is within the bound, break
//         if (stable && err <= errUppBound)
//             break;
//         // if unstable or error is beyond the bound, roll back
//         Recov(net, replInfo, false);
//         newLacNum >>= 1;
//     }

//     // remove redundant nodes
//     net.CleanUp();
//     net.Check();
// }


void ALSMan::ExactSimpl(NetMan & net) {
    cout << endl << "******************** simplify ********************" << endl;
    cout << "before simplification: ";
    net.PrintStat();

    cout << "after simplification: ";
    if (net.GetNetType() == NET_TYPE::GATE) {
        net.SynthAndMap(maxDelay, false);
    }
    else {
        net.Comm("st; compress2rs; logic; sop; ps;");
        // net.SynthWithResyn2Comm();
    }
    net.MergeConst();
    cout << "**********************************************************" << endl << endl;
}


double ALSMan::Eval(NetMan& net, double err, bool useYosys, bool IsInitialCircuit) {
    assert(net.GetNetType() == NET_TYPE::SOP);

    // measure and output SOP
    ostringstream oss("");
    oss << outpPath << round << "_" << net.GetNet()->pName << "_" << metrType << "_" << err << "_size_" << net.GetArea() << "_depth_" << net.GetDelay();
    net.WriteBlif(oss.str() + ".blif");

    // synthesize and technology mapping with ABC
    static double recArea = numeric_limits<double>::max();
    static double recDelay = numeric_limits<double>::max();
    double tempNetArea = DBL_MAX, tempNetDelay = DBL_MAX;
    {
        auto tempNet = net;
        if (IsInitialCircuit)
            tempNet.Comm("st; dch; amap;");
            // tempNet.Comm("st; map;"); // for EPFL only, because the initial circuit is delay optimized
        else
            tempNet.Compile(maxDelay);
            // tempNet.CompileNew(maxDelay); // for EPFL only, because the initial circuit is delay optimized
        // update best
        tempNetArea = tempNet.GetArea();
        tempNetDelay = tempNet.GetDelay();
        if (DoubleLess(tempNetArea, recArea) || 
            (DoubleEqual(tempNetArea, recArea) && DoubleLess(tempNetDelay, recDelay)) ) {
            recArea = tempNetArea;
            recDelay = tempNetDelay;
        }
    }

    // synthesize and technology mapping with Yosys
    // if (!IsInitialCircuit && useYosys) { // (optional,) for obtaining better results
    //     auto tempNet= net;
    //     tempNet.CompileWithYosys(standCellPath);
    //     tempNetArea = tempNet.GetArea();
    //     tempNetDelay = tempNet.GetDelay();
    //     if (DoubleLess(tempNetDelay, maxDelay) || DoubleGreat(recDelay, maxDelay)) { 
    //         if (DoubleLess(tempNetArea, recArea) || 
    //         (DoubleEqual(tempNetArea, recArea) && DoubleLess(tempNetDelay, recDelay)) ) {
    //             recArea = tempNetArea;
    //             recDelay = tempNetDelay;
    //         }
    //     }
    // }

    // output and return
    cout << "current best: area = " << recArea << ", delay = " << recDelay << endl;
    return recDelay;
}


static void AppendNet(Abc_Ntk_t* pResNtk, Abc_Ntk_t* pNtkAcc, Abc_Ntk_t* pNtkApp, Abc_Ntk_t* pNtkMit, int netMark, RETURN_VAR IntVect& miterId2AppId, RETURN_VAR IntVect& appId2MiterId) {
    // init
    // netMark = 0, accNet; netMark = 1, appNet; netMark = 2, mitNet;
    Abc_Ntk_t* pNtkDeal = nullptr;
    if (netMark == 0)
        pNtkDeal = pNtkAcc;
    else if (netMark == 1)
        pNtkDeal = pNtkApp;
    else if (netMark == 2)
        pNtkDeal = pNtkMit;
    else
        assert(0);
    assert(!Abc_NtkIsStrash(pNtkDeal));
    Abc_Obj_t* pObj = nullptr;
    Abc_Obj_t* pFanin = nullptr;
    int i = 0, k = 0;
    Abc_NtkCleanCopy(pNtkDeal);

    // deal with PIs
    if (netMark == 0) {
        Abc_NtkForEachPi(pNtkDeal, pObj, i)
            abc::Abc_NtkDupObj(pResNtk, pObj, 1);
    }
    else if (netMark == 1) {
        Abc_NtkForEachPi(pNtkDeal, pObj, i) {
            auto pPi = abc::Abc_NtkPi(pResNtk, i);
            pObj->pCopy = pPi;
            miterId2AppId[pPi->Id] = pObj->Id;
            appId2MiterId[pObj->Id] = pPi->Id;
        }
    }
    else if (netMark == 2) {
        Abc_NtkForEachPo(pNtkAcc, pObj, i)
            abc::Abc_NtkPi(pNtkDeal, i)->pCopy = abc::Abc_ObjChild0Copy(pObj);
        int nWidth = abc::Abc_NtkPoNum(pNtkAcc);
        Abc_NtkForEachPo(pNtkApp, pObj, i)
            abc::Abc_NtkPi(pNtkDeal, i + nWidth)->pCopy = abc::Abc_ObjChild0Copy(pObj);
    }
    else
        assert(0);
    // duplicate nodes
    Abc_NtkForEachNode(pNtkDeal, pObj, i) {
        if (pObj->pCopy == nullptr) {
            auto pNewNode = abc::Abc_NtkDupObj(pResNtk, pObj, 0);
            if (netMark == 1) {
                miterId2AppId[pNewNode->Id] = pObj->Id;
                appId2MiterId[pObj->Id] = pNewNode->Id;
            }
            if (netMark == 0)
                RenameAbcObj(pObj->pCopy, string(Abc_ObjName(pObj)) + "_acc");
            else if (netMark == 1)
                RenameAbcObj(pObj->pCopy, string(Abc_ObjName(pObj)) + "_app");
            else if (netMark == 2)
                RenameAbcObj(pObj->pCopy, string(Abc_ObjName(pObj)) + "_dev");
            else
                assert(0);
        }
    }
    // reconnect all nodes
    Abc_NtkForEachNode(pNtkDeal, pObj, i) {
        Abc_ObjForEachFanin(pObj, pFanin, k)
            Abc_ObjAddFanin(pObj->pCopy, pFanin->pCopy);
    }
    // deal with POs
    if (netMark == 2) {
        Abc_NtkForEachPo(pNtkDeal, pObj, i)
            abc::Abc_NtkDupObj(pResNtk, pObj, 1);
        Abc_NtkForEachPo(pNtkDeal, pObj, i)
            Abc_ObjAddFanin(pObj->pCopy, abc::Abc_ObjChild0Copy(pObj));
    }
}


NetManPtr ALSMan::BuildErrorRateMiter(NetMan& accNet, NetMan& appNet, RETURN_VAR IntVect& miterId2AppId, RETURN_VAR IntVect& appId2MiterId) {
    // check
    assert(IsPIOSame(accNet, appNet));
    assert(accNet.GetNetType() == NET_TYPE::SOP && appNet.GetNetType() == NET_TYPE::SOP);

    // start empty network
    auto pResNet = make_shared<NetMan>();
    pResNet->StartSopNet();

    // start XOR-OR unit
    int nBit = accNet.GetPoNum();
    auto pDevNet = BuildXorOrCircuit(nBit);
    // pXorOrNet->Sweep();

    // init miterId2AppId and appId2MiterId
    miterId2AppId.resize(accNet.GetIdMaxPlus1() + appNet.GetIdMaxPlus1() + pDevNet->GetIdMaxPlus1(), -1);
    appId2MiterId.resize(appNet.GetIdMaxPlus1(), -1);

    // copy accNet
    AppendNet(pResNet->GetNet(), accNet.GetNet(), appNet.GetNet(), pDevNet->GetNet(), 0, miterId2AppId, appId2MiterId);
    AppendNet(pResNet->GetNet(), accNet.GetNet(), appNet.GetNet(), pDevNet->GetNet(), 1, miterId2AppId, appId2MiterId);
    AppendNet(pResNet->GetNet(), accNet.GetNet(), appNet.GetNet(), pDevNet->GetNet(), 2, miterId2AppId, appId2MiterId);

    // return
    return pResNet;
}


NetManPtr ALSMan::BuildXorOrCircuit(int nBits) {
    // check
    assert(nBits > 0);

    // start empty network
    auto pResNet = make_shared<NetMan>();
    pResNet->StartSopNet();
    auto pAbcNet = pResNet->GetNet();

    // create PIs
    // cout << "create PIs" << endl;
    AbcObjVect piA(nBits, nullptr), piB(nBits, nullptr);
    for (int i = 0; i < nBits; ++i) {
        // cout << "create PIA " << i << endl;
        piA[i] = abc::Abc_NtkCreatePi(pAbcNet);
    }
    for (int i = 0; i < nBits; ++i) {
        // cout << "create PIB " << i << endl;
        piB[i] = abc::Abc_NtkCreatePi(pAbcNet);
    }

    // create XOR gates
    // cout << "create XOR gates" << endl;
    AbcObjVect xorGates(nBits, nullptr);
    for (int i = 0; i < nBits; ++i)
        xorGates[i] = pResNet->CreateNode(AbcObjVect{piA[i], piB[i]}, "01 1\n10 1\n");
    
    // create OR gate
    // cout << "create OR gate" << endl;
    string orFunc = "";
    for (int i = 0; i < nBits; ++i)
        orFunc += "0";
    orFunc += " 0\n";
    auto pOr = pResNet->CreateNode(xorGates, orFunc);

    // create POs
    // for (int i = 0; i < nBits; ++i) {
    //     auto pPo = abc::Abc_NtkCreatePo(pAbcNet);
    //     abc::Abc_ObjAddFanin(pPo, diff[i]);
    // }
    auto pPo = abc::Abc_NtkCreatePo(pAbcNet);
    abc::Abc_ObjAddFanin(pPo, pOr);

    // return
    return pResNet;
}


NetManPtr ALSMan::BuildErrorDistanceMiter(NetMan& accNet, NetMan& appNet, RETURN_VAR IntVect& miterId2AppId, RETURN_VAR IntVect& appId2MiterId) {
    // check
    assert(IsPIOSame(accNet, appNet));
    assert(accNet.GetNetType() == NET_TYPE::SOP && appNet.GetNetType() == NET_TYPE::SOP);

    // start empty network
    auto pResNet = make_shared<NetMan>();
    pResNet->StartSopNet();

    // start absolute difference unit 
    int nBit = accNet.GetPoNum();
    auto pAbsDiffNet = BuildAbsoluteDifferenceCircuit(nBit);
    pAbsDiffNet->Sweep();

    // init miterId2AppId and appId2MiterId
    miterId2AppId.resize(accNet.GetIdMaxPlus1() + appNet.GetIdMaxPlus1() + pAbsDiffNet->GetIdMaxPlus1(), -1);
    appId2MiterId.resize(appNet.GetIdMaxPlus1(), -1);

    // copy accNet
    AppendNet(pResNet->GetNet(), accNet.GetNet(), appNet.GetNet(), pAbsDiffNet->GetNet(), 0, miterId2AppId, appId2MiterId);
    AppendNet(pResNet->GetNet(), accNet.GetNet(), appNet.GetNet(), pAbsDiffNet->GetNet(), 1, miterId2AppId, appId2MiterId);
    AppendNet(pResNet->GetNet(), accNet.GetNet(), appNet.GetNet(), pAbsDiffNet->GetNet(), 2, miterId2AppId, appId2MiterId);

    // return
    return pResNet;
}


NetManPtr ALSMan::BuildAbsoluteDifferenceCircuit(int nBits) {
    // check
    assert(nBits > 0);

    // start empty network
    auto pResNet = make_shared<NetMan>();
    pResNet->StartSopNet();
    auto pAbcNet = pResNet->GetNet();

    // create constant nodes
    // cout << "create constant nodes" << endl;
    auto pConst0 = abc::Abc_NtkCreateNodeConst0(pAbcNet);
    auto pConst1 = abc::Abc_NtkCreateNodeConst1(pAbcNet);

    // create PIs
    // cout << "create PIs" << endl;
    AbcObjVect piA(nBits, nullptr), piB(nBits, nullptr);
    for (int i = 0; i < nBits; ++i) {
        // cout << "create PIA " << i << endl;
        piA[i] = abc::Abc_NtkCreatePi(pAbcNet);
    }
    for (int i = 0; i < nBits; ++i) {
        // cout << "create PIB " << i << endl;
        piB[i] = abc::Abc_NtkCreatePi(pAbcNet);
    }

    // extend A and B by a sign bit
    // cout << "extend A and B by a sign bit" << endl;
    AbcObjVect extA(nBits + 1, nullptr), extB(nBits + 1, nullptr);
    for (int i = 0; i < nBits; ++i)
        extA[i] = piA[i];
    extA[nBits] = pConst0;
    for (int i = 0; i < nBits; ++i)
        extB[i] = piB[i];
    extB[nBits] = pConst0;

    // get ~B
    // cout << "get ~B" << endl;
    AbcObjVect negB(nBits + 1, nullptr);
    for (int i = 0; i <= nBits; ++i)
        negB[i] = abc::Abc_NtkCreateNodeInv(pAbcNet, extB[i]);

    // sum = A + ~B + 1
    // ripple carry adder
    AbcObjVect sum(nBits + 1, nullptr), carry(nBits + 1, nullptr);
    // set C[0] = 1
    carry[0] = pConst1;
    for (int i = 0; i <= nBits; ++i) {
        // sum[i] = A[i] ^ B[i] ^ C[i]
        auto aXorB = pResNet->CreateNode(AbcObjVect{extA[i], negB[i]}, "01 1\n10 1\n");
        sum[i] = pResNet->CreateNode(AbcObjVect{aXorB, carry[i]}, "01 1\n10 1\n");
        // discard the carry out of the last bit
        if (i == nBits)
            break;
        // C[i + 1] = (A[i] & B[i]) | (C[i] & (A[i] ^ B[i]))
        auto aAndB = pResNet->CreateNode(AbcObjVect{extA[i], negB[i]}, "11 1\n");
        auto aXorBAndC = pResNet->CreateNode(AbcObjVect{aXorB, carry[i]}, "11 1\n");
        carry[i + 1] = pResNet->CreateNode(AbcObjVect{aAndB, aXorBAndC}, "00 0\n");
    }

    // get comp[nBits - 1: 0] = (~sum[nBits - 1: 0] + 1)
    AbcObjVect comp(nBits, nullptr), carry2(nBits, nullptr);
    carry2[0] = pConst1;
    for (int i = 0; i < nBits; ++i) {
        auto pNegSum = abc::Abc_NtkCreateNodeInv(pAbcNet, sum[i]);
        comp[i] = pResNet->CreateNode(AbcObjVect{pNegSum, carry2[i]}, "01 1\n10 1\n");
        if (i == nBits - 1)
            break;
        carry2[i + 1] = pResNet->CreateNode(AbcObjVect{pNegSum, carry2[i]}, "11 1\n");
    }

    // diff[nBits - 1: 0] = sum[nBits]? comp[nBits - 1: 0]: sum[nBits - 1: 0]
    AbcObjVect diff(nBits, nullptr);
    for (int i = 0; i < nBits; ++i)
        diff[i] = abc::Abc_NtkCreateNodeMux(pAbcNet, sum[nBits], comp[i], sum[i]);

    // create POs
    for (int i = 0; i < nBits; ++i) {
        auto pPo = abc::Abc_NtkCreatePo(pAbcNet);
        abc::Abc_ObjAddFanin(pPo, diff[i]);
    }

    // return
    return pResNet;
}


NetManPtr ALSMan::BuildMiterWithYosys(NetMan& accNet, NetMan& appNet, RETURN_VAR IntVect& miterId2AppId, RETURN_VAR IntVect& appId2MiterId) {
    // check
    assert(IsPIOSame(accNet, appNet));
    assert(accNet.GetNetType() == NET_TYPE::SOP && appNet.GetNetType() == NET_TYPE::SOP);

    // start empty network
    auto pResNet = make_shared<NetMan>();
    pResNet->StartSopNet();

    // start absolute difference unit 
    int nBit = accNet.GetPoNum();
    auto pDevNet = BuildDeviationCircuit(nBit);

    // init miterId2AppId and appId2MiterId
    miterId2AppId.resize(accNet.GetIdMaxPlus1() + appNet.GetIdMaxPlus1() + pDevNet->GetIdMaxPlus1(), -1);
    appId2MiterId.resize(appNet.GetIdMaxPlus1(), -1);

    // copy accNet
    AppendNet(pResNet->GetNet(), accNet.GetNet(), appNet.GetNet(), pDevNet->GetNet(), 0, miterId2AppId, appId2MiterId);
    AppendNet(pResNet->GetNet(), accNet.GetNet(), appNet.GetNet(), pDevNet->GetNet(), 1, miterId2AppId, appId2MiterId);
    AppendNet(pResNet->GetNet(), accNet.GetNet(), appNet.GetNet(), pDevNet->GetNet(), 2, miterId2AppId, appId2MiterId);

    // return
    return pResNet;
}


static void CreateBehLevDeviation(METR_TYPE metrType, bool isSign, int nBits, const string& fileName) {
    FILE* f = fopen(fileName.c_str(), "w");
    fprintf(f, "module deviation(a, b, f);\n");
    fprintf(f, "parameter width = %d;\n", nBits);
    if (isSign) {
        fprintf(f, "input signed [width - 1: 0] a;\n");
        fprintf(f, "input signed [width - 1: 0] b;\n");
    }
    else {
        fprintf(f, "input [width - 1: 0] a;\n");
        fprintf(f, "input [width - 1: 0] b;\n");
    }
    if (metrType == METR_TYPE::MSE)
        fprintf(f, "output [width * 2 - 1: 0] f;\n");
    else if (metrType == METR_TYPE::MHD) {
        int poWidth = (int)(log2(nBits)) + 1;
        fprintf(f, "output [%d: 0] f;\n", poWidth - 1);
    }
    else
        assert(0);
    
    if (metrType == METR_TYPE::MSE) {
        fprintf(f, "wire [width - 1: 0] diff;\n");
        fprintf(f, "assign diff = (a > b)? (a - b): (b - a);\n");
        fprintf(f, "assign f = diff * diff;\n");
    }
    else if (metrType == METR_TYPE::MHD) {
        fprintf(f, "wire [width - 1: 0] diff;\n");
        fprintf(f, "assign diff = a ^ b;\n");
        fprintf(f, "assign f = 1'b0");
        for (int i = 0; i < nBits; ++i)
            fprintf(f, " + diff[%d]", i);
        fprintf(f, ";\n");
    }
    else
        assert(0);
    fprintf(f, "endmodule\n");
    fclose(f);
}


NetManPtr ALSMan::BuildDeviationCircuit(int nBits) {
    // if the miter has been built or loaded, return the miter
    static std::shared_ptr<NetMan> pResNet = nullptr;
    if (pResNet != nullptr)
        return pResNet;

    // check & create folder
    assert(nBits > 0);
    const string folder = "input/miter/";
    CreatePath(folder);

    // get the name of the miter file
    ostringstream fileNameBase;
    if (metrType == METR_TYPE::MSE)
        fileNameBase << folder << (isSign? "signed_": "unsigned_") << "mse_width_" << nBits;
    else if (metrType == METR_TYPE::MHD)
        fileNameBase << folder << "mhd_width_" << nBits;
    else
        assert(0);
    auto behName = fileNameBase.str() + "_beh.v";
    auto finalName = fileNameBase.str() + "_sop.blif";

    // if miter file exists, load the miter file
    if (IsPathExist(finalName)) {
        AbcMan abcMan;
        abcMan.ReadNet(finalName);
        pResNet = make_shared<NetMan>(abcMan.GetNet(), true);
    }
    else { // if miter file doesn't exist, synthesize a new miter file
        CreateBehLevDeviation(metrType, isSign, nBits, behName);
        ostringstream comm;
        comm << "yosys -q -p \"read_verilog " << behName << "; synth; abc -script +source,abc.rc;st;compress2rs;ps; write_blif " << finalName << "\"";
        ExecSystComm(comm.str());
        AbcMan abcMan;
        abcMan.ReadNet(finalName);
        pResNet = make_shared<NetMan>(abcMan.GetNet(), true); 
    }

    // return
    return pResNet;
}



double ALSMan::ComputeError(Simulator& accSmlt, NetMan& net) {
    Simulator appSmlt(net, seed, nFrame);
    appSmlt.InpUnifFast();
    appSmlt.Sim();
    double err = -1;
    if (metrType == METR_TYPE::ER)
        err = accSmlt.GetErrRate(appSmlt);
    else if (metrType == METR_TYPE::MED)
        err = accSmlt.GetMeanErrDist(appSmlt, isSign);
    else
        assert(0);
    return err;
}


double ALSMan::ComputeError(Simulator& accSmlt, Simulator& appSmlt) {
    double err = -1;
    if (metrType == METR_TYPE::ER)
        err = accSmlt.GetErrRate(appSmlt);
    else if (metrType == METR_TYPE::MED)
        err = accSmlt.GetMeanErrDist(appSmlt, isSign);
    else
        assert(0);
    return err;
}


// bool ALSMan::CheckError(NetMan& net, double err, int seedChange) {
//     if (net.GetPoNum() <= 16)
//         return true;
//     seed += seedChange;
//     auto newErr = CalcErr(accNet, net, isSign, seed, nFrame, metrType, distrType);
//     // if (DoubleEqual(err, 0) || DoubleEqual(newErr, 0)) {
//     //     auto dev = fabs(err - newErr);
//     //     if (dev > 0.2) {
//     //         cout << "old " << metrType << " = " << err << ", new " << metrType << " = " << newErr << ", absolute deviation = " << dev << endl;
//     //         cout << "WARNING: unstable " << metrType << " measurement" << endl;
//     //         seed -= seedChange;
//     //         return false;
//     //     }
//     // }
//     // else {
//     //     auto dev = fabs(err - newErr) / err;
//     //     if (dev > 0.2) {
//     //         cout << "old " << metrType << " = " << err << ", new " << metrType << " = " << newErr << ", relative deviation = " << dev << endl;
//     //         cout << "WARNING: unstable " << metrType << " measurement" << endl;
//     //         seed -= seedChange;
//     //         return false;
//     //     }
//     // }
//     auto diff = fabs(err - newErr);
//     if (diff > 1e-3 && diff > errVariationTolerance) {
//         cout << "WARNING: unstable " << metrType << " measurement: old = " << err << ", new = " << newErr << ", |diff| = " << diff << ", tolerance = " << errVariationTolerance << endl;
//         seed -= seedChange;
//         return false;
//     }
//     seed -= seedChange;
//     return true;
// }


Int2DVect TempApplyLacs(NetMan& net, LACPtrVect& lacs, LAC_TYPE lacType, bool isVerb) {
    assert(lacs.size());
    Int2DVect replTraces;
    for (const auto & pLac: lacs) {
        abc::Abc_Obj_t* pTS = nullptr, *pSS = nullptr;
        pTS = net.GetObj(pLac->GetTargId());
        if (lacType == LAC_TYPE::RESUB) {
            auto& specLac = *dynamic_pointer_cast<ResubLAC> (pLac);
            auto newNodeId = net.CreateNode(specLac.GetDivIds(), specLac.GetSop());
            pSS = net.GetObj(newNodeId);
            if (isVerb) {cout << pTS << " is rewritten, "; net.PrintObj(pSS, true);}
        }
        else
            assert(0);
        auto replTrace = net.TempRepl(pTS, pSS);
        replTraces.emplace_back(replTrace);
    }
    return replTraces;
}


void Recov(NetMan& net, Int2DVect& replTraces, bool isVerb) {
    for (auto it = replTraces.rbegin(); it != replTraces.rend(); ++it)
        net.Recov(*it, isVerb);
}


