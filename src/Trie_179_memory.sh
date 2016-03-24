#!/bin/bash
# Trie script

keyFile=~/Documents/research/PCAP/caida_trace/BT_hao_Originalkey_prefix_key

ipfileFolder=/home/hao/classbench_tracefiles/0809_gen_mapped_caida_trace/mapped_caida_trace/

memSize=(200 250 300 320 340 360 380 410 450)

feedbackPortion=(0.01 0.05 0.1 0.2 0.4 0.6)

interval=(2000000 5000000 10000000 20000000 50000000 100000000)

blackKeySize=(500 1000 2000 4000 8000)

for i in $(seq 0 8)

do
	export OMP_NUM_THREADS=8

	./trieNoiseMain ${keyFile} ${memSize[i]} ${ipfileFolder} ${feedbackPortion[2]} ${interval[0]} ${blackKeySize[2]}

done






