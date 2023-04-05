
# Date

Update: Nov. 14, 2014  
Author: You Lizhao  

# Note

This folder contains following functions:

* Tools for generating simulation data;  
* PNC Decoder;  
* NCMA Decoder;  
* SIC Decoder for NCMA;  
* (potential) Pt2Pt OFDM Decoder.  

The latter three decoders are based on RawOFDM parameters (eg., pn preambles).

# How to run?

* Signal_Generator_Main.m: entrance function for generating simulated samples;  

* PNC_Relay_RX_Main.m: entrance function for PNC decoder given input samples;  

* NCMA_phy_sim.m: entrance function for NCMA PHY Simulation; then run NCMA_gen_full_maps.m to get full maps;  

* SIC_phy_sim.m: entrance function for SIC PHY Simulation; then run SIC_gen_phy_maps.m to get fininal maps  

===

# Matlab Decoder for Raw samples
* File: debug_pnc_signal_out.m  
* Input: raw samples (rx-sigmix.dat, rx-sync.datb) or (rx-sig.dat)   
* TODO: implement Find_Sync_Index for rx-sig.dat mode  


# Matlab Decoder for FFT signals
* File: debug_pnc_fft_output.m  
* Input: rx-ofdmdemod.dat (fft results); rx-ofdmdemod.datb (flag) 

===

# Folders

* ncma_txdata1/: old transmission data (for NCMA paper);  
* ncma_txdata2/: transmission data (for new experiments);  
* benchmark_rxdata/: benchmark reception data;  
* ncma_demod_data/: data after sampler & fft (for NCMA paper) - too large to put here, please contact me;  
* python-sim/: simulater for NCMA MAC-layer;  

## Difference between ncma_txdata1 and ncma_txdata2

* 16S: amp(txdata1)*10=amp(txdata2); for user A, the LTS of ncma_txdata1 is different from ncma_txdata2  
* 64S: they are the same;  
* 256S: they are the same;  
* 512S: they are the same;  

# Matlab Decoder Functions

* PNC_Packet_Decoder.m: Decoder for single packet;  

* PNC_Data_Decoder.m: Decoder for data symbols;  

* PNC_Subcarrier_Decoder.m: Decoder for single subcarrier;  

# rawofdm.dat
* Training symbols for RawOFDM/PNC
* TS1: time-domain short training symbol for time sync, cfo estimation (same for A and B)
* TS2_freq_a: freq-domain As long training symbol for channel estimation
* TS2_freq_b: freq-domain Bs long training symbol for channel estimation 

