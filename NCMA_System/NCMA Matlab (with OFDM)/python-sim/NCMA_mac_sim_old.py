#!/usr/bin/python
#
# @author: You Lizhao
#
# Try to calculate fourth-equation system (with RS code) throughput
# We only have python implementation of fourth-equation system
#

from ncma_functions import phy_get_maps
from ncma_functions import mac_count_equations, mac_update_seqno, mac_find_last_index

fft_size   = 64
data_tones = 48
nsym       = 64

RANGE_L = [4,8,16,32]
RANGE_L = [16]
nl = len(RANGE_L)

RANGE_RATIO = [1,1.5,2,3,4,5]
RANGE_RATIO = [1]
nratio = len(RANGE_RATIO)

nindex = 10
[A_MUD_INDEX,B_MUD_INDEX,
 A_NCMA1_INDEX,B_NCMA1_INDEX,
 A_NCMA2_INDEX,B_NCMA2_INDEX,
 A_NCMA3_INDEX,B_NCMA3_INDEX,
 A_NCMA4_INDEX,B_NCMA4_INDEX] = range(0,nindex)
print [A_MUD_INDEX,B_MUD_INDEX,
 A_NCMA1_INDEX,B_NCMA1_INDEX,
 A_NCMA2_INDEX,B_NCMA2_INDEX,
 A_NCMA3_INDEX,B_NCMA3_INDEX,
 A_NCMA4_INDEX,B_NCMA4_INDEX]

npktsAll = [[0 for col in range(nindex)] for row in range(nl*nratio)]

nindex = 8
[MUD_INDEX, NCMA_MINUS_INDEX, 
 NCMA_INDEX, NCMA_UPPER_INDEX,
 NCMA_PLUS_INDEX, NCMA_PLUS_UPPER_INDEX, 
 NCMA_PLUS2_INDEX, NCMA_PLUS2_UPPER_INDEX] = range(0,nindex)

npkts = [[0 for col in range(nindex)] for row in range(nl*nratio)]

SLOTS = 1000
[A_INDEX,B_INDEX,X_INDEX]=range(3)

ss='python_data/rmud_phy/0128_011_75dB'
[phyRawMap,phyBridgeMap,phyTwoPhaseMap]=phy_get_maps(ss)
num_pkts=len(phyRawMap)
#print num_pkts
#-----------------------------------------------------------#
#                     Main Funcations                       #
#-----------------------------------------------------------#
for rr in range(nratio):
  ratio = RANGE_RATIO[rr]

  for ll in range(nl):
    L = RANGE_L[ll]
    L_A = round(L*ratio)
    L_B = L
    L_X = max(L_A,L_B)

    print " rr=%d ll=%d L=%d ratio=%f \n" % (rr,ll,L,ratio)

    d_ncma2_X_decodable = 0
    d_ncma3_X_decodable = 0
    d_ncma3_X_decodable = 0
    
    mudMACMap   = [[0]*3 for ii in range(num_pkts)] # MUD
    ncma1MACMap = [[0]*3 for ii in range(num_pkts)] # NCMA-: MUD + PHY Bridging
    ncma2MACMap = [[0]*3 for ii in range(num_pkts)] # NCMA : MUD + PHY/MAC Bridging
    ncma3MACMap = [[0]*3 for ii in range(num_pkts)] # NCMA+: MUD + PHY Bridging + MAC 4-th Equs 
    ncma4MACMap = [[0]*3 for ii in range(num_pkts)] # NCMA++: NCMA+ + PHY 2-pass decoding 

    mudCntA=0; mudCntB=0;
    ncma1CntA=0; ncma1CntB=0;
    ncma2CntA=0; ncma2CntB=0; ncma2CntX=0;
    ncma3CntA=0; ncma3CntB=0; ncma3CntX=0;
    ncma4CntA=0; ncma4CntB=0; ncma4CntX=0;

    mudPktNumA=1; mudPktNumB=1;
    ncma1PktNumA=1; ncma1PktNumB=1;
    ncma2PktNumA=1; ncma2PktNumB=1; ncma2PktNumX=1;
    ncma3PktNumA=1; ncma3PktNumB=1; ncma3PktNumX=1;
    ncma4PktNumA=1; ncma4PktNumB=1; ncma4PktNumX=1;

    N = 255
    L = num_pkts
    #C = rs_code.RSCode(N,L,systematic=0)
    matrix = mu_rs_code.NCMA()
    for pktIndex in range(num_pkts):
      #---------------------------------------------#
      #               MUD MAC Layer                 #
      #---------------------------------------------#
      okA=phyRawMap[pktIndex][A_INDEX]
      okB=phyRawMap[pktIndex][B_INDEX]
      if okA:
        mudCntA += 1
        mudMACMap[pktIndex][A_INDEX] = mudPktNumA
      if okB:
        mudCntB += 1
        mudMACMap[pktIndex][B_INDEX] = mudPktNumB
      if mudCntA == L_A:
        mudPktNumA += 1
        mudCntA = 0
      if mudCntB == L_B:
        mudPktNumB += 1
        mudCntB = 0

      #---------------------------------------------#
      #          NCMA-/NCMA MAC Layer               #
      #---------------------------------------------#
      okA=phyBridgeMap[pktIndex][A_INDEX]
      okB=phyBridgeMap[pktIndex][B_INDEX]
      okX=phyBridgeMap[pktIndex][X_INDEX]
      if okA:
        ncma1CntA += 1
        ncma1MACMap[pktIndex][A_INDEX] = ncma1PktNumA
        ncma2CntA += 1
        ncma2MACMap[pktIndex][A_INDEX] = ncma2PktNumA
      if okB:
        ncma1CntB += 1
        ncma1MACMap[pktIndex][B_INDEX] = ncma1PktNumB
        ncma2CntB += 1
        ncma2MACMap[pktIndex][B_INDEX] = ncma2PktNumB
      if okX:
        ncma2CntX += 1
        ncma2MACMap[pktIndex][X_INDEX] = ncma2PktNumX

      # NCMA-
      if ncma1CntA == L_A:
        ncma1PktNumA += 1
        ncma1CntA = 0
      if ncma1CntB == L_B:
        ncma1PktNumB += 1
        ncma1CntB = 0

      # NCMA
      matrix.PacketReception(pktIndex,okA,okB,okX)

      # Check NCMA equations
      inSeq=[ncma2CntA,ncma2PktNumA,ncma2CntB,ncma2PktNumB,ncma2CntX,ncma2PktNumX]
      #print pktIndex,inSeq,d_ncma2_X_decodable,L_A,L_B,ncma2MACMap
      #print [ncma2CntA,ncma2PktNumA,ncma2CntB,ncma2PktNumB,ncma2CntX,ncma2PktNumX]
      [ncma2CntA,ncma2PktNumA,ncma2CntB,ncma2PktNumB,ncma2CntX,ncma2PktNumX,d_ncma2_X_decodable] = \
          mac_update_seqno(pktIndex,inSeq,matrix,d_ncma2_X_decodable,ncma2MACMap,L_A,L_B,'ncma2 ')
      #print [ncma2CntA,ncma2PktNumA,ncma2CntB,ncma2PktNumB,ncma2CntX,ncma2PktNumX]

  cur = rr*nl+ll
  
  npktsAll[cur][A_MUD_INDEX]   = 1.0*((mudPktNumA-1)*L_A+mudCntA)/SLOTS;
  npktsAll[cur][B_MUD_INDEX]   = 1.0*((mudPktNumB-1)*L_B+mudCntB)/SLOTS;
  npktsAll[cur][A_NCMA1_INDEX] = 1.0*((ncma1PktNumA-1)*L_A+ncma1CntA)/SLOTS;
  npktsAll[cur][B_NCMA1_INDEX] = 1.0*((ncma1PktNumB-1)*L_B+ncma1CntB)/SLOTS;
  npktsAll[cur][A_NCMA2_INDEX] = 1.0*((ncma2PktNumA-1)*L_A+ncma2CntA)/SLOTS;
  npktsAll[cur][B_NCMA2_INDEX] = 1.0*((ncma2PktNumB-1)*L_B+ncma2CntB)/SLOTS;
  npktsAll[cur][A_NCMA3_INDEX] = 1.0*((ncma3PktNumA-1)*L_A+ncma3CntA)/SLOTS;
  npktsAll[cur][B_NCMA3_INDEX] = 1.0*((ncma3PktNumB-1)*L_B+ncma3CntB)/SLOTS;
  npktsAll[cur][A_NCMA4_INDEX] = 1.0*((ncma4PktNumA-1)*L_A+ncma4CntA)/SLOTS;
  npktsAll[cur][B_NCMA4_INDEX] = 1.0*((ncma4PktNumB-1)*L_B+ncma4CntB)/SLOTS;

  npkts[cur][MUD_INDEX] = npktsAll[cur][A_MUD_INDEX]+npktsAll[cur][B_MUD_INDEX];
  npkts[cur][NCMA_MINUS_INDEX] = npktsAll[cur][A_NCMA1_INDEX]+npktsAll[cur][B_NCMA1_INDEX];
  npkts[cur][NCMA_INDEX]       = npktsAll[cur][A_NCMA2_INDEX]+npktsAll[cur][B_NCMA2_INDEX];
  npkts[cur][NCMA_PLUS_INDEX]  = npktsAll[cur][A_NCMA3_INDEX]+npktsAll[cur][B_NCMA3_INDEX];
  npkts[cur][NCMA_PLUS2_INDEX] = npktsAll[cur][A_NCMA4_INDEX]+npktsAll[cur][B_NCMA4_INDEX];

  [two_equations_num,one_equations_num]=mac_count_equations(phyRawMap);
  npkts[cur][NCMA_UPPER_INDEX] = 1.0*(two_equations_num*2+one_equations_num)/SLOTS;

  print "ncma1: pktinfo=",[ncma1PktNumA,ncma1CntA,ncma1PktNumB,ncma1CntB]
  print "ncma2: pktinfo=",[ncma2PktNumA,ncma2CntA,ncma2PktNumB,ncma2CntB]
  print "round: %d %d %d" % (cur, rr, ll)
  print "npktsAll=",npktsAll[cur]
  print "npkts=",npkts[cur]
