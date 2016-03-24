
/**
trie_cuckoo.cpp
create by: Yunhong
create time: 05/25/2015
*/


#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>
#include <sys/time.h>
#include <arpa/inet.h>
//#include <itpp/itbase.h>
#include "trie.h"
using namespace std;

typedef std::vector<std::string> strings;
typedef std::vector<size_t> size_ts;
typedef std::vector<int> ints;
typedef std::vector<vector<int> > int2s;
typedef std::vector<vector<long> > long2s;
typedef std::vector<long> longs;
typedef std::vector<float> floats;
typedef std::vector<double> doubles;
typedef std::vector<size_t> size_ts;
typedef std::vector<vector<size_t> > size_tss;

extern CuckooFilter cuckooFilter;
extern CuckooFilter cuckooFilterInit0;
extern CuckooTable cuckooBlackKeyTable;
extern CuckooTable cuckooTableKey;
extern CuckooTable cuckooAggrKeyTable;
extern CuckooFilter cuckooFilterFlowEst;


int initCuckoo(vector<string> &key,vector<int> &keyPrefix,
                vector<int> &keyAction,float &storage, int& finger, char mL0[][4][20]);

bool addCuckooFilter(vector<string> &keys, vector<int> &keyPrefixes,
                      vector<int> &keyActions);

bool addCuckooFilter(vector<string> &keys, vector<int> &keyActions);

bool addCuckooFilter0(vector<string> &keys, vector<int> &keyPrefixes,
                      vector<int> &keyActions, char mL0[][4][20]);

bool kickoutBlack(vector<size_t> &mask, vector<string> &flow,
                  vector<string> &flow_cnt, vector<string> &keyprefixlength,
                  vector<string> &testIP, vector<string> &testIP_cnt,
                  vector<string> &flowprelen, double &blackkey_sum);

bool aggregateTrie(Trie *bTrie, int ai, ActionOrder *actionOrder, size_t &countkey,
                   size_t &countaggregatekey, size_t &countblackkey,
                   size_t &countorikey, vector<string> &key,
                   vector<int> &keyaction,vector<string> &blackkey,
                   vector<int> &blackkeyPrefixes, vector<string> &aggregatekey, int& prefixlength,
                   bool isPrint, bool isInit);


bool comWeightKeytype(vector<string> &flow_cnt, vector<string> &keytypestr,
                      double &pkt_sum, double &key_sum,
                      vector<double> &weight,vector<keyType> &keytype);

bool clusterAction(vector<string> &flow, vector<int> &flowaction,
                   vector<size_t> &mask, vector<size_t> &index);

bool sortUniqueFlow(vector<string> &flow_in, vector<size_t> &flow_cnt_in);

bool lookupKey(size_t &line, int actionSize, size_t &flow_int,size_t &flow_cnt_int, vector<size_t> &mask,
               double &key_sum, double &pkt_sum, double &aggrSum, floats& keySums, floats& countIps, double &countIP,
               size_t &countNum, double &countIP0, size_t &countNum0, size_t &countblack_tmp,
               vector<string> &overbigKey, vector<size_t> &overbigKey_cnt,  Trie* trie, int2s& bigNonCounts, long2s& timeCounts,
               vector<int> &UniquePrefix,vector<int> &UniqueAggPrefix);

bool assignAction(vector<string> &key,vector<int> &keyaction,int &actionSize);

size_t initAggregation(vector<string> &keyin,vector<int> &keyprefix,
                        vector<int> &keyactionin,vector<size_t> &mask, int actionSize,
                        float &storage, bool isInit, int& finger,
                        vector<int> &UniqueAggPrefix);

size_t aggregation(vector<string> &keyIns,vector<int> &keyPrefixIns,
                        vector<int> &keyActionIns, vector<size_t> &maskes,
                      int actionSize, float &storage, bool isInit,
                      int &fingerprintOld,vector<int> &uniqueAggKeyprefixes, strings& overKeys, size_ts& overKeyNos,
                      floats& haoOvers, float target);

void assignAction(vector<string> &flow,vector<size_t> &flow_cnt,
                  vector<int> &flowaction,vector<size_t> &pkt_num_action,
                  int actionSize);


bool insertWordTrieSimple(Trie *bTrie, int trieNum, vector<size_t> &index,
                          vector<string> &flow,vector<int> &keyprefixlength,
                          vector<int> &flowaction);

void prefixNum(vector<int> &keyprefix, vector<int> &UniquePrefix);

void keySort(vector<string>& keys, vector<size_t>& keyNos);

void keySortInt(vector<string>& keys, vector<int>& keyNos);

void addFlowEst(string& key, int2s& bigNonCounts, long2s& timeCounts, strings& overBigKeys );

void compressAct(floats& haoOvers, float target, int actionSize, ints& compressActions, ints& decompressActions);

//void keySort(vector<string>& keys, vector<int>& keyNos, vector<int>& keyPrefixes);
