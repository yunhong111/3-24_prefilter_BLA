/**
trie_cuckoo_noise_07_21.cpp
create by: Yunhong
create time: 07/21/2015
*/


#include "Trie_cuckoo_noise_07_21.h"
// 1 bit flag, 2 bits actions, f bits fingerprint for filters

const float load_factor = 0.91;
const int per_key_prefilter = 40;

const int actionSize = 4;

// time
double runTime = 0.0;

// packet generation rate: pps
const float rateParameter0 = 297298.0;

float get_total_memory()
{
	// KB
	float storage = 360;
	
	// convert to bits
	storage *= 8.0*1024;

	cout<<"total_memory: "<<storage<<endl;
	return storage;
}

float get_total_prefilter_memory()
{
	size_t slot_num = 2000; // 5006 11640+
	float storage = slot_num*per_key_prefilter/load_factor;

	cout<<"+++++++++++++++++++++++++++++"<<(storage/(8*1024))<<endl;

	cout<<"prefilter+black memory: "<<storage<<endl;
	return storage;
}

// prefilter is a cuckoo table with size of
// table: key+action+1 bit flag = (32+5) + 2 + 1 = 40
// load factor = 
float get_prefilter_size(size_t other_key_size)
{
	float prefilter_size = other_key_size;

	// get prefilter size(bits) = slots*per_slot_size
	float storage = prefilter_size*per_key_prefilter/load_factor;

	cout<<"get_prefilter_size: "<<storage<<endl;
	return storage;
}

float get_black_memory(size_t other_key_size)
{
	float prefilter_memory = get_prefilter_size(other_key_size);

	float total_pre_memory = get_total_prefilter_memory();

	float black_memory = total_pre_memory - prefilter_memory;

	if(black_memory <= 0)
	{
		cout<<"get_black_memory::error::total prefilter memory is too small!";
		exit(0);
	}

	cout<<"black_memory: "<<black_memory<<endl;
	return black_memory;
}

//
float get_act_filter_size()
{
	float total_memory = get_total_memory();

	float total_pre_memory = get_total_prefilter_memory();

	float act_filter_memory = total_memory - total_pre_memory;

	if(act_filter_memory <= 0)
	{
		cout<<"get_act_filter_size::error::total memory is too small!";
		exit(0);
	}

	cout<<"act_filter_memory: "<<act_filter_memory<<endl;
	return act_filter_memory;
}

void init_act_filter(vector<string> &keys, int& finger)
{
	float storage = get_act_filter_size(); // bits
	
    //-----------------------------------------------------------
    // Parameters for cuckoo filter
    float loadFactor = load_factor;
    int slotNo = 4;
    size_t keySize = keys.size();
    size_t bucketNo = size_t(keySize/(loadFactor*slotNo))+1;
    float fingerprint0 = storage*loadFactor/(keys.size())-3;
    int fingerprint = storage*loadFactor/(keys.size())-3;
    cout<<(storage/(8*1024))<<" ++++++++++++++++++++++++++++ "<<keys.size()<<endl;
    finger = fingerprint;
    long maxNumKicks = 1500;
    cout<<"* Fingerprint length: "<<fingerprint0<<" "<<fingerprint<<endl;
    //mL0 =vector<vector<string> > (bucketNo, vector<string>(slotNo, " "));

    // --------------------------------------------------------
    //init cuckoo filter
    cuckooFilter.ClearTable();
    cuckooFilter.cuckooFilterInit(bucketNo,fingerprint,slotNo,
                                       maxNumKicks); 
}

void add_act_filter(vector<string> &keys,vector<int> &keyActions)
{
	int finger;
	init_act_filter(keys, finger);
                
    // define variables
    cout<<"* Cuckoo filter adding ..."<<endl;
    bool flagAdd;
    int countFail = 0;
    string str;
    size_t keySize = keys.size();

    for(size_t i = 0; i < keySize; i++)
    {
        str = keys[i];
        //add keys to cuckoo filter
        flagAdd = cuckooFilter.AddKey(str,(keyActions[i]));
        if(flagAdd == 0)
        {
            countFail++;
            cout << "* 664 Cuckoo filter flag_add fail"<<endl;
            cout<<"* Order: "<<i<<"  ";
            break;
        }

    }
    cout<<"* Line 670 Count fail num: "<<countFail<<endl;
}

void init_prefilter1(CuckooTable& cuckooPretable, size_t other_key_size)
{
	float prefilter_memory= get_prefilter_size(other_key_size);

	// prefilter table slot number
	size_t prefilter_size = prefilter_memory*load_factor/per_key_prefilter;
	cout<<"prefilter slots: "<<prefilter_size<<endl;
	
	// init prefilter table
	float a = load_factor;
	int f = 12;
	int bc = 4;
	long MaxNumKicks = 1000;
    long m = prefilter_size/(a*bc)+1;
    
    cuckooPretable.ClearTable();
    cuckooPretable.CuckooTableInit(m,f,bc,
                                        MaxNumKicks);
}

// other_key: ipv4/prefix
void add_prefilter1(CuckooTable& cuckooPretable, strings& other_keys, ints& other_keyactions)
{
	init_prefilter1(cuckooPretable, other_keys.size());

	for(size_t i = 0; i < other_keys.size(); i++)
	{
		bool addFlag = cuckooPretable.AddKeyPrefix(other_keys[i],32, other_keyactions[i]);
		if(!addFlag)
		{
			cout<<"2 add fail..."<<endl;
		}
	}
}

void init_black_table(size_t other_key_size)
{
	float black_memory = get_black_memory(other_key_size);

	// black table slot number
	size_t black_size = black_memory*load_factor/per_key_prefilter;
	cout<<"black table slots: "<<black_size<<endl;
	
	// init black table
	float a = load_factor;
	int f = 12;
	int bc = 4;
	long MaxNumKicks = 1000;
    long m = black_size/(a*bc)+1;
    
    cuckooBlackKeyTable.ClearTable();
    cuckooBlackKeyTable.CuckooTableInit(m,f,bc,
                                        MaxNumKicks);
}

void add_aggr_table(strings& aggregateKeys)
{
	// ----------------------------------------
    // Add aggregate keys to cuckooTable
    int fingerprint = 12;
    long maxNumKicks = 1000;
    int slotNo = 4;
    size_t aggregateKeySize = aggregateKeys.size();
    long bucketNo = long(aggregateKeySize/(load_factor*slotNo));
    cuckooAggrKeyTable.CuckooTableInit(bucketNo,fingerprint,slotNo,maxNumKicks);

    for(size_t i = 0; i < aggregateKeySize; i++)
    {
        size_t found = aggregateKeys[i].find('/');
        string str = aggregateKeys[i].substr(0,found);
        string prefixstr = aggregateKeys[i].substr(found+1,
                                                   aggregateKeys[i].size()-found);

        cuckooAggrKeyTable.AddKey(str,str2num(prefixstr));
    }
}
void load_keys(char* infile_name, vector<string> &keys, vector<int> &keyActions)
{
	string key_ipv4;
	int key_action;
	ifstream infile_key(infile_name);
	//ofstream outfile("../test/key_actions_0_1");
	cout<<"load_keys+++++++++++++++"<<endl;
	while(infile_key >> key_ipv4>>key_action)
	{
		keys.push_back(key_ipv4);
		keyActions.push_back(key_action);

		//cout<<key_ipv4<<" "<<key_action<<endl;
	}
	infile_key.close();
}

double nextTime(double rateParameter)
{
    double nTime = -log(1.0f - (double) rand() / double(RAND_MAX + 1.0)) / rateParameter;

    if(nTime == INFINITY)
    {
        double nTime = -log(1.0f - (double) rand() / (RAND_MAX + 1.0)) / rateParameter;
        cout<<"* nextTime: INFINITY!"<<endl;
        exit(0);
    }
    return nTime;
}

int main(int argc, char * argv[])
{
    cout << "\n";
    cout << "PRIME_OPENMP\n";
    cout << "  C++/OpenMP version\n";

    cout << "\n";
    cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
    cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";

    // ---------------------------
    // Global variables
    CUCKOO_BLACK_SIZE = strtof(argv[6],NULL); // black table size
    FLOW_EST_SIZE = 0;       // flow estimator size
    float overTarget = 0.01;    // target overselection rate

    // -----------------------------
    /* initialize random seed: */
    srand (time(NULL));

    // time v.r.
    struct timeval gen_start, gen_end; /* gettimeofday stuff */
    struct timezone tzp;
    long timeInv = 0;

    // ------------------------------
    /*load key file*/
    std::string inFileName = argv[1];
    std::ifstream infile(inFileName.c_str());
    if(!infile)
        std::cout << "Train File Error " << std::endl;

    int flowPrefixInt;
    string flowStr;
    vector<int> keyprefixlengths;
    vector<string> keys;
    while(infile >> flowPrefixInt >> flowStr)
    {
        keyprefixlengths.push_back(flowPrefixInt);
        keys.push_back(flowStr);
    }
    cout<<"* Key size: "<<keys.size()<<endl;
    infile.clear();
    infile.close();

    // --------------------------------
    // get unique prefix length
    vector<int> uniquePrefix;
    vector<int> uniqueAggPrefix;
    prefixNum(keyprefixlengths, uniquePrefix);

    // ------------------------------
    /*assign actions*/
    int actionSize = 4;
    vector<int> keyActions;
    assignAction(keys,keyActions,actionSize);

    // --------------------------------
    /*cuckoo table*/
    // load factor
    // m: key number, f: fingerprint width, bc: slots num per bucket,
    // MaxNumKicks: max kickout number in a cuckoo filter
    cout<<"* Init cuckoo table ... ..."<<endl<<endl;
    float a;
    int f,bc;
    a = 0.9;
    bc = 4;
    long m = long(keys.size()/(a*bc))+1;
    f = 12;
    long MaxNumKicks = 1000;
    cuckooTableKey.ClearTable();
    cuckooTableKey.CuckooTableInit(m,f,bc,MaxNumKicks);

    // --------------------------------
    // Add original key to cuckoo table
    cout<<"* Add key to cuckoo table ... ..."<<endl;
    bool isAddTable;
    for (int i = 0; i < keys.size(); i++)
    {
        isAddTable = cuckooTableKey.AddKeyPrefix(keys[i],(keyprefixlengths[i]),keyActions[i]);

        if(isAddTable == 0)
        {
            cout<<"* Flag_add fail"<<endl;
            cout<<"* Order: "<<i<<"  ";
        }

    }

    // -----------------------------------------------
    // init cuckooFilter for flow estimation
    long flowEstSize = FLOW_EST_SIZE;
    m = flowEstSize/(a*bc)+1;
    f = 14;
    cuckooFilterFlowEst.ClearTable();
    cuckooFilterFlowEst.cuckooFilterInit(m,f,bc,MaxNumKicks);

    // init black table
    m = CUCKOO_BLACK_SIZE/(a*bc)+1;
    cuckooBlackKeyTable.ClearTable();
    cuckooBlackKeyTable.CuckooTableInit(m,f,bc,
                                        MaxNumKicks);
    // -----------------------------------
    // Init blackkey file
    cout<<"* Write blackkey to file!"<<endl;
    ofstream blackKeyFileOut;
    BLACKFILENAME = "../test_no_noise_3/blackkeyfile_" + string(argv[2]) + '_' +
                          string(argv[4])+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6]);
    blackKeyFileOut.open(BLACKFILENAME.c_str());
    blackKeyFileOut.clear();
    blackKeyFileOut.close();
    uint16_t blackBackSize = 0;
    float feedSumPortion = 0.0;

    // -----------------------------------
    // Init aggr file
    cout<<"* Write aggr to file!"<<endl;
    ofstream aggrFileOut;
    AGGRFILENAME = "../test_no_noise_3/aggrfile" + string(argv[2]) + '_' +
                          string(argv[4])+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6]);
    aggrFileOut.open(AGGRFILENAME.c_str());
    aggrFileOut.clear();
    aggrFileOut.close();

    // -----------------------------------------------
    // Two counters for estimator
    int2s bigNonCounts;
    long2s timeCounts;
    bigNonCounts= vector<vector<int> > (m, vector<int>(bc, 0));
    timeCounts = vector<vector<long> > (m, vector<long>(bc, 0));

    // ---------------------------------------------
    /* init cuckoo filter without aggregation*/
    cout<<"* Init cuckoo filter ... ..."<<endl;

    float storage = strtof(argv[2],NULL); // storage size
    int finger = 0;
    long but = 200000/3;
    char mL0[but][4][20];//
    bzero(&mL0,sizeof(mL0));
    //cout<<sizeof(mL0)<<endl;
    //return 0;
    vector<vector<size_t> > keyCountcs = vector<vector<size_t> > (but, vector<size_t>(4, 0));;
    vector<vector<size_t> > keyCountcs0= vector<vector<size_t> > (but, vector<size_t>(4, 0));;
    vector<vector<size_t> > keyCountDiffs= vector<vector<size_t> > (but, vector<size_t>(4, 0));;
    initCuckoo(keys,keyprefixlengths,keyActions,storage, finger, mL0);
    int finger0 = finger;

    /*for(size_t i = 0; i < 100; i++)
    {
        for(int j=0; j <4; j++)
            cout<<mL0[i][j]<<" ";
    }
    cout<<endl;*/

    // ---------------------------------------------------
    /*define mask*/
    vector<size_t> mask;
    size_t maskIP;
    mask.clear();
    cout<<"* Compute prefixes main ... ..."<<endl<<endl;;
    for(int i = 8; i <= 32; i++)
    {
        maskIP = (size_t(pow(2,i))-1)<<(32-i);
        mask.push_back(maskIP);
    }

    // ------------------------------------------------
    // Init aggregation
    bool isInit = 1;
    //initAggregation(keys,keyprefixlengths,keyActions,\
                    mask, actionSize, storage, isInit, finger,uniqueAggPrefix);
                    
    // add keys to actual filter
    // load keys
    vector<string> aggr_keys;
    vector<int> aggr_keyActions;
	char* infile_name = (char*)"../key/key_out_0_3_wo";
    
    load_keys(infile_name, aggr_keys, aggr_keyActions);
    cout<<aggr_keys.size()<<endl;

    add_act_filter(aggr_keys, aggr_keyActions);

    // ----------------------------
    // Compute unique prefix length
    vector<int> aggr_keyPrefixes;
    //ints uniqueAggPrefix;
    for(size_t i = 0; i< aggr_keys.size(); i++)
    {
        size_t found = aggr_keys[i].find('/');
        string str = aggr_keys[i].substr(0,found);
        string prefixStr = aggr_keys[i].substr(found+1, aggr_keys[i].size()-found);
        aggr_keyPrefixes.push_back(str2num(prefixStr));
    }

    prefixNum(aggr_keyPrefixes,uniqueAggPrefix);
    
    // add keys tp prefilter
    vector<string> other_keys;
    vector<int> other_keyActions;

    char* other_infile_name = (char*)"../key/other_key_out_0_3_wo";
    
    load_keys(other_infile_name, other_keys, other_keyActions);

	CuckooTable cuckooPretable;
    add_prefilter1(cuckooPretable, other_keys, other_keyActions);

    ints uni_prefilter_prefix;
    // ----------------------------
    // Compute unique prefix length
    vector<int>other_keyPrefixes;
    //ints uniqueAggPrefix;
    for(size_t i = 0; i< other_keys.size(); i++)
    {
        size_t found = other_keys[i].find('/');
        string str = other_keys[i].substr(0,found);
        string prefixStr = other_keys[i].substr(found+1, other_keys[i].size()-found);
        other_keyPrefixes.push_back(str2num(prefixStr));
    }
    prefixNum(other_keyPrefixes, uni_prefilter_prefix);

    // init_black_table(other_keys.size());

    // add to aggr table
	strings aggregateKeys;
	ints aggr_actions;

	char* aggr_infile_name = (char*)"../key/aggr_key_out_0_3_wo";
    
    load_keys(aggr_infile_name, aggregateKeys, aggr_actions);
	
    add_aggr_table( aggregateKeys);

    // ------------------------------------------------
    // File name for caida trace
    const char * fp[] = {"mapped-caida-1",
                         "mapped-caida-6",
                         "mapped-caida-11",
                         "mapped-caida-16",
                         "mapped-caida-21",
                         "mapped-caida-26",
                         "mapped-caida-31",
                         "mapped-caida-36",
                         "mapped-caida-41",
                         "mapped-caida-46",
                         "mapped-caida-51",
                         "mapped-caida-56"
                        };

    // -------------------------------------------------
    // Open out file, and write result into it
    std::ofstream outfile0;
    string outfileName = ("../test_no_noise_3/outfile0_simple_"+ string(argv[2]) + '_' +
                          string(argv[4]))+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6])+".csv";
    outfile0.open(outfileName.c_str());

    // Open out file, and write result into it
    std::ofstream keyCountfile;
    outfileName = ("../test_no_noise_3/keyCountfile_simple_"+ string(argv[2]) + '_' +
                          string(argv[4]))+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6])+".csv";
    keyCountfile.open(outfileName.c_str());
    strings keycs;
    vector<int> keyActioncs;

    strings keycsM;
    vector<int> keyActioncsM;
    vector<size_t> keyCountcsM;
    vector<size_t> keyCountcs0M;
    vector<size_t> keyCountDiffsM;
    cuckooFilter.returnKey(keycsM,keyActioncsM,mL0, keyCountcs,keyCountcs0,keyCountDiffs,keyCountcsM,
                                       keyCountcs0M,keyCountDiffsM);
    for(size_t i = 0; i < keycsM.size(); i++)
    {
        keyCountfile<<keycsM[i]<<" ";
    }
    keyCountfile<<endl;
    for(size_t i = 0; i < keycsM.size(); i++)
    {
        keyCountfile<<keyActioncsM[i]<<" ";
    }
    keyCountfile<<endl;
    for(size_t i = 0; i < keycsM.size(); i++)
    {
        keyCountfile<<keyCountcsM[i]<<" ";
    }
    keyCountfile<<endl;
    for(size_t i = 0; i < keycsM.size(); i++)
    {
        keyCountfile<<keyCountcs0M[i]<<" ";
    }
    keyCountfile<<endl;
    for(size_t i = 0; i < keycsM.size(); i++)
    {
        keyCountfile<<keyCountDiffsM[i]<<" ";
    }
    keyCountfile<<endl;

    // -----------------------------------------------
    float falsePos = 0;
    float falsePos0 = 0;
    float haoFalsePos = 0;
    float haoFalsePos0 = 0;
    float haoFalsePosTotal;
    float haoFalsePos0Total;
    float overAggr = 0.0f;

    vector<string> overBigKeys;
    vector<size_t> overBigKeyNos;
    vector<string> blackKeys;
    vector<size_t> blackKeyNos;
    vector<float> haoOvers;
    haoOvers.assign(actionSize,0);

    size_t countNum = 0;
    size_t countNum0 = 0;
    double countIP = 0;
    double countIP0 = 0;
    size_t countBlack = 0;

    double countIPTotal = 0.0f;
    double countIP0Total = 0.0f;
    double keySum = 0.0f;
    double pktSum = 0.0f;
    double keySumTotal = 0.0f;
    double pktSumTotal = 0.0f;
    double countIPTotalOff = 0.0f;
    double countIP0TotalOff = 0.0f;
    double keySumTotalOff = 0.0f;
    double pktSumTotalOff = 0.0f;
    double aggrSum= 0.0f;

    floats keySums;
    floats countIPs;
    keySums.assign(actionSize,0);
    countIPs.assign(actionSize,0);

    uint64_t line = 0;

    // ---------------------------
    // Define a trie
    Trie *trie;            // define tree
    trie = new Trie();            // define tree

    // -------------------------
    //std::ifstream infile;

    // ----------------------------
    // load packets from file
    for (int fi = 0; fi < 6; fi++)
    {
        // --------------------------
        // Open trace file
        string pathname = fp[fi];
        std::string inFileName = argv[3];
        inFileName += pathname;
        infile.open(inFileName.c_str());
        cout<<inFileName.c_str()<<endl;
        cout<<inFileName<<endl;
        if(!infile)
            std::cout << "* TestIP File Error " << std::endl;

        // ------------------------------
        countIPTotalOff += countIP;
        countIP0TotalOff += countIP0;
        keySumTotalOff += keySum;
        pktSumTotalOff += pktSum;

        //line = 0;
        countBlack = 0;
        countIP = 0.0f;
        countIP0 = 0.0f;
        keySum = 0.0f;
        pktSum = 0.0f;
        //aggrSum = 0.0f;
        keySums.clear();
        countIPs.clear();
        keySums.assign(actionSize,0);
        countIPs.assign(actionSize,0);

        int nthreads, tid;

        size_t  ei;
        bool isEndFlag = 0;
        size_t updateInvDis = 10000; // interval for display
        size_t readNum = strtof(argv[5],NULL);    // interval for feedback

        while(!isEndFlag )
        {
            // -------------------------------
            // read file
            strings flows;
            size_ts flowNos;
            readFile0(infile, flows, flowNos, readNum, isEndFlag);
            size_t flowNum = flows.size();

            # pragma omp parallel for \
            shared ( infile,isEndFlag, updateInvDis, line, actionSize, mask, keySum, pktSum, aggrSum, keySums, countIPs, countIP, countNum, countIP0, countNum0, \
                     countBlack, overBigKeys,overBigKeyNos,  trie, bigNonCounts, timeCounts, uniquePrefix,uniqueAggPrefix, cuckooFilter, \
                     cuckooFilterInit0,cuckooBlackKeyTable,cuckooTableKey,cuckooAggrKeyTable,cuckooFilterFlowEst) \
            private ( ei )

            for(ei = 0; ei < flowNum; ei++)
            {
                // read 1 flow
                bool readFlag1 = 0;
                size_t flowNoInt = 1;
                uint32_t flowInt;

                #pragma omp critical
                {
                    /*if(!readFlag1)
                    {
                        readFlag1 = (infile>>flowInt>>flowSizeInt);
                    }*/
                    flowInt = parseIPV4string(flows[ei].c_str());
                    flowNoInt = flowNos[ei];
                    readFlag1 = 1;

                }

                if(readFlag1)
                {

                    // ---------------------------
                    // ---------------------------
                    // lookup  key
                    {
                        int prefix;
                        int flag_look,flag_look0, flag_lookkey;
                        bool flag_lookblack = 0;
                        uint32_t subIP;
                        string flowstr;
                        string flowIpv4 = parsedec2IPV4(flowInt);

                        uint32_t ip = flowInt; //convert IP to int
                        string keyType_cur = "0";
                        string flow_action_str;

                        int keySumLocal = 0;
                        int aggrSumLocal = 0;
                        int countIPLocal = 0;
                        int countNumLocal = 0;
                        int countIP0Local = 0;
                        int countNum0Local = 0;
                        vector<int> keySumsLocal;
                        keySumsLocal.assign(actionSize,0);
                        vector<int> countIPsLocal;
                        countIPsLocal.assign(actionSize,0);

                        //pkt_sum += flowNoInt;

                        // -------------------------------
                        // look up blackkey
                        int iflowaction;
                        if(CUCKOO_BLACK_SIZE > 0)
                        {
							
                            int mi = 24;
                            //for(int mi = 0; mi <= 32-8; mi++)
                            {

                                subIP = ip & mask[mi]; // mask
                                flowstr = parsedec2IPV4(ip);//s.str();
                                prefix = mi+8;

                                flag_lookblack = cuckooBlackKeyTable.LookUpKeyAction(flowstr,prefix,iflowaction);
                                if (flag_lookblack)
                                {
                                    countBlack += flowNoInt;
                                    //break;
                                }
                                //cout<<"check black key: "<<flag_lookblack<<endl;
                            }


                        }

                        // ------------------------------
                        flag_look = 0;
                        if (flag_lookblack == 0) // not a blackkey
                        {
                            // -------------------------------
                            // look up key
                            for(int mi = uniquePrefix.size()-1; mi >=0; mi--)
                            {
                                subIP = ip & mask[uniquePrefix[mi]-8];
                                flowstr = parsedec2IPV4(subIP);
                                prefix = uniquePrefix[mi];

                                flag_lookkey = cuckooTableKey.LookUpKeyActionCount(flowstr,prefix,iflowaction,flowNoInt);
                                if (flag_lookkey)
                                {

                                    keyType_cur = "1";
                                    flow_action_str = num2str(iflowaction);
                                    //#pragma omp atomic
                                    keySumLocal += (flowNoInt);
                                    for(int ai = 0; ai <actionSize; ai++)
                                    {
                                        if(iflowaction == ai)
                                        {
                                            //#pragma omp atomic
                                            keySumsLocal[ai] += (flowNoInt);
                                        }

                                    }
                                    break;

                                }
                                
                            }

                            // lookup prefilter
                            bool lookup_prefilter = 0;
                            for(int mi = uni_prefilter_prefix.size()-1; mi >=0; mi--)
                            {
                                subIP = ip & mask[uni_prefilter_prefix[mi]-8];
                                flowstr = parsedec2IPV4(subIP);
                                prefix = uni_prefilter_prefix[mi];

                                lookup_prefilter = cuckooTableKey.LookUpKeyActionCount(flowstr,prefix,iflowaction,flowNoInt);
                                if (lookup_prefilter)
                                {

                                    countNumLocal += lookup_prefilter;
                                    countIPLocal += lookup_prefilter*(flowNoInt);

                                    // each action lookups
                                    for(int ai = 0; ai <actionSize; ai++)
                                    {
 
                                       if(iflowaction == ai)
                                       {
                                           countIPsLocal[ai] += (flowNoInt);
                                       }

                                        
                                    }
                                    break;

                                }
                                
                            }
							if(!lookup_prefilter)
							{
								// -------------------------------
								// look up aggregate key
								bool isAggregatekey = 0;
								//if(cuckooAggrKeyTable.mm > 1)
								{
								   /* for(int mi = uniqueAggPrefix.size()-1; mi >=0; mi--)
									{
										
										
										subIP = ip & mask[uniqueAggPrefix[mi]-8];
										flowstr = parsedec2IPV4(subIP);
										prefix = uniqueAggPrefix[mi];
										cout<<"check filter: "<<isAggregatekey<<" prefix: "<<(uniqueAggPrefix[mi]-8)<<endl;
										isAggregatekey = cuckooAggrKeyTable.LookUpKey(flowstr,prefix);
										 cout<<"check filter: "<<isAggregatekey<<" prefix: "<<(uniqueAggPrefix[mi]-8)<<endl;
										if (isAggregatekey && !flag_lookkey)
										{
											//#pragma omp atomic
											//cout<<"check filter: "<<flag_look<<endl;
											aggrSumLocal += (flowNoInt);
											break;

										}
										cout<<"check filter: "<<isAggregatekey<<" prefix: "<<(uniqueAggPrefix[mi]-8)<<endl;
									}*/
								}

								// -------------------------------
								// lookup key from cuckoo filter
								vector<int> iactions;
								bool isAKey = 0;

								for(int mi = 0; mi < uniqueAggPrefix.size(); mi++)
								{
									bool overbigFlag = 0;
									subIP = ip & mask[uniqueAggPrefix[mi]-8];
									flowstr = parsedec2IPV4(subIP)+"/"+num2str(uniqueAggPrefix[mi]);
									iactions.clear();
									flag_look = cuckooFilter.LookUpKeyActionsCount(flowstr,iactions,flowNoInt, mL0,
																					keyCountcs,  keyCountcs0, keyCountDiffs);
									if (flag_look != 0)
									{
										
										isAKey = 1;
										countNumLocal += flag_look;
										countIPLocal += flag_look*(flowNoInt);

										// each action lookups
										for(int ai = 0; ai <actionSize; ai++)
										{
											for(int li = 0; li < iactions.size(); li++)
											{
												if(iactions[li] == ai)
												{
													countIPsLocal[ai] += (flowNoInt);
												}

											}
										}

										// overselection count
										if(!flag_lookkey && !overbigFlag)
										{
											overbigFlag = 1;
											#pragma omp critical
											{
												trie->addWordCountNum(DecToBin(flowInt),32, isnonkey, 0, flowNoInt);
											}

										}
									}
									//cout<<"check filter: "<<flag_look<<" prefix: "<<(uniqueAggPrefix[mi]-8)<<endl;
								}

								// ---------------------------------------------------
								// add to flow est filter
								if(isAKey == 0)  // not inside cuckooFilter
								{
									// add to filter
									//if(FLOW_EST_SIZE>0)
								   // addFlowEst(flowIpv4, bigNonCounts, timeCounts, overbigKey );

								}
							}

                        }

                        // -------------------------------
                        // lookup key from cuckoo filterInit0
                        vector<int> iactions;
                        for(int mi = 0; mi < uniquePrefix.size(); mi++)
                        {
                            subIP = ip & mask[uniquePrefix[mi]-8];
                            flowstr = parsedec2IPV4(subIP)+"/"+num2str(uniquePrefix[mi]);
                            iactions.clear();
                            flag_look0 =cuckooFilterInit0.LookUpKeyActions(flowstr,iactions);
                            if (flag_look0)
                            {
                                //#pragma omp atomic
                                countNum0Local += flag_look0;
                                //#pragma omp atomic
                                countIP0Local += flag_look0*(flowNoInt);

                            }

                        }
                        // lookup key end
                        // ---------------------------------------------
                        // ---------------------------------------------

                        // ---------------------------------------------
                        // update global variable

                        {
                            line ++;
                            #pragma omp atomic
                            pktSum += flowNoInt;
                            #pragma omp atomic
                            keySum += keySumLocal;
                            #pragma omp atomic
                            aggrSum += aggrSumLocal;
                            #pragma omp atomic
                            countIP += countIPLocal;
                            #pragma omp atomic
                            countNum += countNumLocal;
                            #pragma omp atomic
                            countIP0 += countIP0Local;
                            #pragma omp atomic
                            countNum0 += countNum0Local;
                            for(int ai = 0; ai <actionSize; ai++)
                            {
                                #pragma omp atomic
                                keySums[ai] += keySumsLocal[ai];
                                #pragma omp atomic
                                countIPs[ai] += countIPsLocal[ai];
                            }
                        }

                    }
                }
                else
                {
                    isEndFlag = 1; // end read flag
                }
                //tid = omp_get_thread_num();

                // ---------------------------------
                // update rates and write to file

                #pragma omp critical
                {
                    if(size_t(line)%200 == 0)
                    {
                        // update total overselects value
                        countIPTotal = countIPTotalOff + countIP;
                        countIP0Total = countIP0TotalOff + countIP0;
                        keySumTotal = keySumTotalOff + keySum;
                        pktSumTotal = pktSumTotalOff + pktSum;

                        // --------------------------------------------
                        // overselection rate
                        // ----------------------------------------------
                        if ((pktSum - keySum) != 0)
                        {
                            falsePos = float(countIP-keySum)/float(pktSum-keySum);
                            falsePos0 = float(countIP0-keySum)/float(pktSum-keySum);
                        }


                        if(keySum != 0)
                        {
                            haoFalsePos = float(countIP-keySum)/float(keySum);
                            haoFalsePos0 = float(countIP0-keySum)/float(keySum);
                        }

                        for(int ai = 0; ai < actionSize; ai++)
                        {
                            if(keySums[ai] != 0)
                            {
                                haoOvers[ai] = float(countIPs[ai]-keySums[ai])/float(keySums[ai]);
                            }
                        }

                        if(keySumTotal != 0)
                        {
                            haoFalsePosTotal = (countIPTotal-keySumTotal)/(keySumTotal);
                            haoFalsePos0Total = (countIP0Total-keySumTotal)/(keySumTotal);
                        }

                        overAggr = aggrSum/keySumTotal;
                    }

                    
                }

            }

            #pragma omp barrier

			// ------------------------------------------------
			// Display and write to file
			outfile0<<"line,"<<line<<",total,"<<haoFalsePosTotal<<",total0,"<<haoFalsePos0Total<<",false,"<<
					falsePos<<",over,"<<haoFalsePos<<",false0,"<<falsePos0<<",over0,"<<
					haoFalsePos0<<",overaggr,"<<overAggr<<",overcuckoo,"<<(haoFalsePosTotal-overAggr)<<",";

			outfile0<<"countIP,"<<countIPTotal<<",countIP0,"<<countIP0Total<<",keysum,"<<keySumTotal<<
					",pktSum,"<<pktSumTotal<<",aggrSum,"<<aggrSum<<",blackkey_num,"<<countBlack<<",feedback,"<<blackBackSize<<",feedsumportion,"
					<<feedSumPortion<<",finger0,"<<finger0<<",finger,"<<finger<<",time,"<<runTime<<endl;

			cout<<endl<<"line,"<<line<<",total,"<<haoFalsePosTotal<<",total0,"<<haoFalsePos0Total<<",false,"<<
					falsePos<<",over,"<<haoFalsePos<<",false0,"<<falsePos0<<",over0,"<<
					haoFalsePos0<<",overaggr,"<<overAggr<<",overcuckoo,"<<(haoFalsePosTotal-overAggr)<<endl;

			cout<<"countIP,"<<countIPTotal<<",countIP0,"<<countIP0Total<<",keysum,"<<keySumTotal<<
					",pktSum,"<<pktSumTotal<<",aggrSum,"<<aggrSum<<",blackkey_num,"<<countBlack<<",feedback,"<<blackBackSize<<",feedsumportion,"
					<<feedSumPortion<<",finger0,"<<finger0<<",finger,"<<finger<<",time,"<<runTime<<endl;
		
            // ---------------------------------------------
            // feed back portionFeedBack% overselections and overselection rates for actions
            //if(size_t(line)%2000000 == 0 )
            if(size_t(runTime)%5 == 0)
            {

                // ---------------------------
                // the overselects for feedback
                vector<char> word;
                trie->getLeaf(trie->root,word,blackKeys,blackKeyNos);
                trie->deleteChild(trie->root);
                delete trie;
                trie = new Trie();            // define tree

                int init = 0;
                float overTotalSum = accumulate(blackKeyNos.begin(), blackKeyNos.end(), init);

                // sort balckkeys
                keySort(blackKeys,blackKeyNos);

                float portionFeedBack = strtof(argv[4],NULL);
                blackBackSize = portionFeedBack*blackKeys.size();
                /*if(blackKeys.size()<100)
                {
                    blackBackSize = blackKeys.size();
                }
                else if(blackKeys.size()<200 && portionFeedBack<=0.5)
                {
                    blackBackSize = 0.5*blackKeys.size();
                }
                else if(blackKeys.size()<1000 && portionFeedBack <=0.1)
                {
                    blackBackSize = 0.1*blackKeys.size();
                }*/
                blackKeys.erase(blackKeys.begin()+blackBackSize, blackKeys.end());
                blackKeyNos.erase(blackKeyNos.begin()+blackBackSize, blackKeyNos.end());

                float feedSum = accumulate(blackKeyNos.begin(), blackKeyNos.end(), init);
                feedSumPortion = feedSum/overTotalSum;

                // ------------------------------------
                // feedback
                // -------------------------
                //time before calling function
                gettimeofday(&gen_start, &tzp);
                // ---------------------------------
                // call function
                feedbackBlackkey(blackKeys);

                // ------------------------------
                // time after calling function
                gettimeofday(&gen_end, &tzp);
                // time interval
                timeInv = print_elapsed("Aggr: ", &gen_start,&gen_end, 1);

                // write to file
                strings().swap(keycsM);
                ints().swap(keyActioncsM);
                size_ts().swap(keyCountcsM);
                size_ts().swap(keyCountcs0M);
                size_ts().swap(keyCountDiffsM);

                cuckooFilter.returnKey(keycsM,keyActioncsM,mL0, keyCountcs,keyCountcs0,keyCountDiffs,keyCountcsM,
                                       keyCountcs0M,keyCountDiffsM);
                for(int i = 0; i < 10; i++)
                cout<<keyCountcsM[i]<<" "<<keyCountcs0M[i]<<" "<<keyCountDiffsM[i]<<"     ";
                cout<<endl;
                for(size_t i = 0; i < keycsM.size(); i++)
                {
                    keyCountfile<<keyCountcsM[i]<<" ";
                }
                keyCountfile<<endl;
                for(size_t i = 0; i < keycsM.size(); i++)
                {
                    keyCountfile<<keyCountcs0M[i]<<" ";
                }
                keyCountfile<<endl;
                for(size_t i = 0; i < keycsM.size(); i++)
                {
                    keyCountfile<<keyCountDiffsM[i]<<" ";
                }
                keyCountfile<<endl;

                // ------------------------------------
                // print haoOvers
                /*for(int i = 0; i< actionSize; i++)
                    cout<<"Action: "<<i <<" haoOver: "<<" "<<haoOvers[i]<<" ";
                cout<<endl;*/

                // ------------------------------------
                // aggregation
                //isInit = 1;
                /*if(line >= 20000000 && line%20000000 == 0 && haoFalsePosTotal > overTarget)
                {

                //trie->deleteChild(trie->root);
                //delete trie;
                //trie = new Trie();            // define tree

                    aggregation(keys,keyprefixlengths,keyActions, \
                            mask, actionSize, storage, isInit, finger,uniqueAggPrefix, blackKeys, blackKeyNos, haoOvers, overTarget);
                }*/


                strings().swap(blackKeys);
                vector<size_t>().swap(blackKeyNos);
            }

            strings().swap(flows);
            vector<size_t>().swap(flowNos);

        }
        //ifstream().swap(infile);
        infile.clear();
        infile.close();
    } //for fi

    //-------------------------------------------------------------
    /*clear the data structure*/
    mask.clear();

    cuckooFilter.ClearTable();
    cuckooBlackKeyTable.ClearTable();
    cuckooTableKey.ClearTable();
    cuckooAggrKeyTable.ClearTable();

    (keyCountfile.close());
    (outfile0.close());
}

void feedbackBlackkey(vector<string>& overBigKeys)
{
    // Add blackkeys to cuckooTable
    GLOBAL_BIGNONKEYNUM = 100.0;
    vector<string> blackKeys;
    blackKeys.clear();

    //cout<<"* overBigKeys No: "<<overBigKeys.size()<<endl;

    for(int i = 0; i < overBigKeys.size(); i++)
    {
        //if(overBigKeyNos[i]>GLOBAL_BIGNONKEYNUM)
        {
            blackKeys.push_back(overBigKeys[i]);
        }
    }

    // Add blackkey to cuckooTable
    //cout<<"* Add balckkey to cuckooTable!"<<endl;
    float loadFactor = 0.9f;
    int slotNo = 4;
    size_t blackKeySize = CUCKOO_BLACK_SIZE;
    size_t bucketSize = int(blackKeySize/(loadFactor*slotNo))+1;
    int fingerprintNew = 12;
    long MaxKickoutNum = 1000;
    cuckooBlackKeyTable.ClearTable();
    cuckooBlackKeyTable.CuckooTableInit(bucketSize,fingerprintNew,slotNo, \
                                        MaxKickoutNum);

    for(size_t i = 0; i < blackKeys.size(); i++)
    {
        if(i<CUCKOO_BLACK_SIZE)
        {
            bool addFlag = cuckooBlackKeyTable.AddKeyPrefix(blackKeys[i],32, 4);
            if(!addFlag)
            {
                cout<<"1 add fail..."<<endl;
            }
        }
    }
    size_t newBlackSize = blackKeys.size();

    // Import previous blackkeys
    ifstream blackKeyFileIn;
    blackKeyFileIn.open(BLACKFILENAME.c_str());
    string blackKeyStr;
    int cuckooBlackSize = CUCKOO_BLACK_SIZE;
    int prefix = 32;
    while(blackKeyFileIn>>blackKeyStr>>prefix && blackKeys.size()<cuckooBlackSize)
    {
        bool blackFlag = 0;
        int iflowaction;
        blackFlag = cuckooBlackKeyTable.LookUpKeyAction(blackKeyStr,prefix,iflowaction);

        if(!blackFlag)
        blackKeys.push_back(blackKeyStr);
    }
    blackKeyFileIn.clear();
    blackKeyFileIn.close();

    // -----------------------------------
    // Write blackkey to file
    //cout<<"* Write blackkey to file!"<<endl;
    ofstream blackKeyFileOut;
    blackKeyFileOut.open(BLACKFILENAME.c_str());
    for(int i = 0; i < blackKeys.size(); i++)
        blackKeyFileOut<<blackKeys[i]<<" "<<32<<endl;
    blackKeyFileOut.clear();
    blackKeyFileOut.close();

    // --------------------------------------
    //cout<<"balckkey  size: "<<blackKeys.size()<<endl;
    for(size_t i = newBlackSize; i < blackKeys.size(); i++)
    {
        bool addFlag = cuckooBlackKeyTable.AddKeyPrefix(blackKeys[i],32, 4);
        if(!addFlag)
        {
            cout<<"2 add fail..."<<endl;
        }
    }



}


void feedbackBlackkey1(strings& overBigKeys)
{

    float loadFactor = 0.9f;
    int slotNo = 4;
    int cuckooBlackSize = CUCKOO_BLACK_SIZE;
    int blackKeySize = cuckooBlackSize;
    int bucketSize = int(blackKeySize/(loadFactor*slotNo))+1;
    int fingerprintNew = 12;
    int MaxKickoutNum = 1000;
    cuckooBlackKeyTable.ClearTable();
    cuckooBlackKeyTable.CuckooTableInit(bucketSize,fingerprintNew,slotNo,
                                        MaxKickoutNum);
    int prefix = 32;
    // write to file
    ofstream blackkeyFile;
    blackkeyFile.open("blackkeyFile");

    sort(overBigKeys.begin(),overBigKeys.end());
    vector<string>::iterator it;
    it = unique (overBigKeys.begin(), overBigKeys.end());
    overBigKeys.resize( distance(overBigKeys.begin(),it) );
    cout<<"* overBigKeys No: "<<overBigKeys.size()<<endl;
    int mini = 0;
    if(overBigKeys.size()>cuckooBlackSize)
    {
        mini = overBigKeys.size()-cuckooBlackSize;
        overBigKeys.erase(overBigKeys.begin(),overBigKeys.begin()+mini);
    }
    for(int i = overBigKeys.size()-1; i > 0 ; i--)
    {

        cuckooBlackKeyTable.AddKeyPrefix(overBigKeys[i],prefix, 4);
        blackkeyFile<<overBigKeys[i]<<endl;
    }
    blackkeyFile.close();

}

bool readFile0(ifstream& infile, vector<string> &flow, vector<size_t> &flow_cnt, size_t readNum, bool& isEndFlag)
{
    cout<<"* Read File!"<<endl;
    Trie *bTrie = new Trie();
    uint32_t flowInt;
    size_t flowNo;
    size_t in_num = 0;

    double timeInv = 0.0f;
    double nTime = 0.0f;

    while((infile >> flowInt >>flowNo) && timeInv < 1.0f/*&& (in_num < readNum)*/)
    {

        bTrie->addWordCount(DecToBin(flowInt),32, isnonkey, flowNo);
        in_num ++;

        // generate next time
        nTime = nextTime(rateParameter0);
        timeInv += nTime;

        if(nTime< 0)
        {
            cout<<"* readFile0: timeInv < 0. wrong!"<<endl;
            exit(0);
        }

    }

    cout<<"* Read File end!"<<" timeInv: "<<timeInv<<endl;
    // run time
    runTime += timeInv;
    if(runTime >= INFINITY)
    {
        cout<<"* readFile0:runTime: "<<runTime<<" timeInv: "<<timeInv<<" nTime: "<<nTime<<endl;
        exit(0);
    }
    //cout<<"* timeInv: "<<endl;

    vector<char> word;
    vector<int> flowActions;

    //vector<string>().swap(flow);
    //vector<size_t>().swap(flow_cnt);
    vector<int>().swap(flowActions);

    bTrie->getLeaf(bTrie->root,word,flow,flow_cnt);


    if(timeInv < 1.0f/*in_num < readNum*/)
    {
        isEndFlag = 1;
        cout<<"* readFile0:timeInv: "<<timeInv<<endl;
    }

    //cout<<"* Flow size: "<<flow.size()<<endl;

	bTrie->deleteChild(bTrie->root);
    delete bTrie;
    return true;
}


/*
bool readFile0(ifstream& infile, vector<string> &flow, vector<size_t> &flow_cnt, size_t readNum, bool& isEndFlag)
{
    Trie *bTrie = new Trie();
    uint32_t flowInt;
    size_t flowNo;
    size_t in_num = 0;

    while((infile >> flowInt >>flowNo) && (in_num < readNum))
    {

        bTrie->addWordCount(DecToBin(flowInt),32, isnonkey, flowNo);
        in_num ++;
        //if(in_num%1000000 == 0)
            //cout<<"loading file ...."<<in_num<<"  ";
    }
    vector<char> word;
    bTrie->getLeaf(bTrie->root,word,flow,flow_cnt);


    if(in_num < readNum)
    {
        isEndFlag = 1;
    }
    //cout<<"* Flow size: "<<flow.size()<<endl;
	bTrie->deleteChild(bTrie->root);
    delete bTrie;
    return true;
}*/
