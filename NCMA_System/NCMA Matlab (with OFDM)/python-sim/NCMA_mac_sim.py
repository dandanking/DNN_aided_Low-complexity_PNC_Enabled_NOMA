#!/usr/bin/python
#
# @author: You Lizhao
#
# Try to calculate fourth-equation system (with RS code) throughput
# We only have python implementation of fourth-equation system
#

import mu_rs_code, mu_rs_code_tes
from ncma_functions import phy_get_maps
from ncma_functions import mac_count_equations #, mac_update_seqno, mac_find_last_index

fwriter = open('ncma_results.dat','w')

RANGE_L = [8,16,32,64]
RANGE_L = [16]
nl = len(RANGE_L)

RANGE_RATIO = [0.5,0.67,1,1.5,2] #,3,4
RANGE_RATIO = [1]
#RANGE_RATIO = [2,3]
nratio = len(RANGE_RATIO)

nindex = 12
[A_MUD1_INDEX, B_MUD1_INDEX,
 A_MUD2_INDEX, B_MUD2_INDEX,
 A_NCMA_TES1_INDEX, B_NCMA_TES1_INDEX,
 A_NCMA_TES2_INDEX, B_NCMA_TES2_INDEX,
 A_NCMA_UES1_INDEX, B_NCMA_UES1_INDEX,
 A_NCMA_UES2_INDEX, B_NCMA_UES2_INDEX] = range(0,nindex)

npktsAll = [[0 for col in range(nindex)] for row in range(nl*nratio)]

nindex = 9
[MUD0_INDEX, MUD1_INDEX, MUD2_INDEX, 
 NCMA_TES1_INDEX, NCMA_UES1_INDEX, NCMA_UES1_UPPER_INDEX,
 NCMA_TES2_INDEX, NCMA_UES2_INDEX, NCMA_UES2_UPPER_INDEX] = range(0,nindex)

npkts = [[0 for col in range(nindex)] for row in range(nl*nratio)]

SLOTS = 1000
[A_INDEX,B_INDEX,X_INDEX]=range(3)

mac_files=['0128_011_75dB','0126_09_8dB','0125_08_85dB','0128_01_9dB','0126_01_95dB','0126_04_10dB','0126_06_105dB']
mac_files=['0128_11_75dB_75dB','0128_013_75dB_95dB','0130_01_75dB_105dB','0128_14_75dB_135dB']
mac_files=['0130_03_95dB_75dB','0130_04_95dB_85dB','0130_05_95dB_95dB','0130_06_95dB_105dB','0130_08_95dB_135dB']
mac_files=['0131_01_01_209dB_123dB','0131_02_01_90dB_123dB','0131_07_01_70dB_74dB','0131_08_01_70dB_90dB','0131_10_74dB_123dB']
#mac_files=['0131_10_74dB_123dB']
#mac_files=['0126_09_8dB']

for ff in range(len(mac_files)):
    ss = mac_files[ff]    #ss='data/rmud_phy/0128_011_75dB'	
    ss = 'data/rmud_phy/' + ss
    print "\nRead file ", ss
    [phyRawMap,phyBridgeMap,phyTwoPhaseMap]=phy_get_maps(ss)
    num_pkts=len(phyRawMap)
    print num_pkts
    #-----------------------------------------------------------#
    #                      Main Functions                       #
    #-----------------------------------------------------------#
    for rr in range(nratio):
      ratio = RANGE_RATIO[rr]
    
      for ll in range(nl):
        L = RANGE_L[ll]
        L_A = int(round(L*ratio))
        L_B = L
        
        print
        print "rr=%d ll=%d L_A=%d L_B=%d ratio=%f \n" % (rr,ll,L_A,L_B,ratio)
    
        N = 255
        K = 1
        #C = rs_code.RSCode(N,L,systematic=0)
        mudMatrix0 = mu_rs_code.MUDSimpleCode(L_A,L_B)
        mudMatrix1 = mu_rs_code.MUDSimpleCode(L_A,L_B)
        mudMatrix2 = mu_rs_code.MUDSimpleCode(L_A,L_B)
        ncmaMatrix_T1 = mu_rs_code_tes.NCMACode(N,K,L_A,L_B,systematic=0)
        #ncmaMatrix_T1.debug2 = True
        ncmaMatrix_T2 = mu_rs_code_tes.NCMACode(N,K,L_A,L_B,systematic=0)
        #print ncmaMatrix_T1, ncmaMatrix_T2
        ncmaMatrix_U1 = mu_rs_code.NCMAUESCode(N,K,L_A,L_B,systematic=0)
        ncmaMatrix_U2 = mu_rs_code.NCMAUESCode(N,K,L_A,L_B,systematic=0)
        #ncmaMatrix1.debug = True
        
        for pktIndex in range(num_pkts):
            okA=phyRawMap[pktIndex][A_INDEX]
            okB=phyRawMap[pktIndex][B_INDEX]
            okX=phyRawMap[pktIndex][X_INDEX]	
            mudMatrix0.RxPacket(okA,okB)  

            #---------------------------------------------#
            #               NCMA MAC Layer                #
            #---------------------------------------------#
            okA=phyBridgeMap[pktIndex][A_INDEX]
            okB=phyBridgeMap[pktIndex][B_INDEX]
            okX=phyBridgeMap[pktIndex][X_INDEX]   
            if not okA and okB and okX:
                okA = 1
            elif okA and not okB and okX:
                okB = 1
            elif okA and okB and not okX:
                okX = 1
             
            mudMatrix1.RxPacket(okA,okB)
            ncmaMatrix_T1.RxPacket(okA,okB,okX)
            ncmaMatrix_U1.RxPacket(okA,okB,okX)
                    
            #---------------------------------------------#
            #              NCMA+ MAC Layer                #
            #---------------------------------------------#
            okA=phyTwoPhaseMap[pktIndex][A_INDEX]
            okB=phyTwoPhaseMap[pktIndex][B_INDEX]
            okX=phyTwoPhaseMap[pktIndex][X_INDEX]     
            if not okA and okB and okX:
                okA = True
            elif okA and not okB and okX:
                okB = True
            elif okA and okB and not okX:
                okX = True
                 
            mudMatrix2.RxPacket(okA,okB)
            ncmaMatrix_T2.RxPacket(okA,okB,okX)
            ncmaMatrix_U2.RxPacket(okA,okB,okX)
          
        [mud0PktNumA, mud0CntA, mud0PktNumB, mud0CntB] = mudMatrix0.GetState()  
        [mud1PktNumA, mud1CntA, mud1PktNumB, mud1CntB] = mudMatrix1.GetState()
        [mud2PktNumA, mud2CntA, mud2PktNumB, mud2CntB] = mudMatrix2.GetState()
        [ncmaT1PktNumA, ncmaT1CntA, ncmaT1PktNumB, ncmaT1CntB] = ncmaMatrix_T1.GetState()      
        [ncmaU1PktNumA, ncmaU1CntA, ncmaU1PktNumB, ncmaU1CntB] = ncmaMatrix_U1.GetState()
        [ncmaT2PktNumA, ncmaT2CntA, ncmaT2PktNumB, ncmaT2CntB] = ncmaMatrix_T2.GetState()
        [ncmaU2PktNumA, ncmaU2CntA, ncmaU2PktNumB, ncmaU2CntB] = ncmaMatrix_U2.GetState()	  
        
        cur = rr*nl+ll
        
        npktsAll[cur][A_MUD1_INDEX] = 1.0*(mud1PktNumA*L_A+mud1CntA)/SLOTS
        npktsAll[cur][B_MUD1_INDEX] = 1.0*(mud1PktNumB*L_B+mud1CntB)/SLOTS
        npktsAll[cur][A_MUD2_INDEX] = 1.0*(mud2PktNumA*L_A+mud2CntA)/SLOTS
        npktsAll[cur][B_MUD2_INDEX] = 1.0*(mud2PktNumB*L_B+mud2CntB)/SLOTS
        npktsAll[cur][A_NCMA_TES1_INDEX] = 1.0*(ncmaT1PktNumA*L_A+ncmaT1CntA)/SLOTS
        npktsAll[cur][B_NCMA_TES1_INDEX] = 1.0*(ncmaT1PktNumB*L_B+ncmaT1CntB)/SLOTS
        npktsAll[cur][A_NCMA_UES1_INDEX] = 1.0*(ncmaU1PktNumA*L_A+ncmaU1CntA)/SLOTS
        npktsAll[cur][B_NCMA_UES1_INDEX] = 1.0*(ncmaU1PktNumB*L_B+ncmaU1CntB)/SLOTS
        npktsAll[cur][A_NCMA_TES2_INDEX] = 1.0*(ncmaT2PktNumA*L_A+ncmaT2CntA)/SLOTS
        npktsAll[cur][B_NCMA_TES2_INDEX] = 1.0*(ncmaT2PktNumB*L_B+ncmaT2CntB)/SLOTS
        npktsAll[cur][A_NCMA_UES2_INDEX] = 1.0*(ncmaU2PktNumA*L_A+ncmaU2CntA)/SLOTS
        npktsAll[cur][B_NCMA_UES2_INDEX] = 1.0*(ncmaU2PktNumB*L_B+ncmaU2CntB)/SLOTS
    
        npkts[cur][MUD0_INDEX] = 1.0*(mud0PktNumA*L_A+mud0CntA)/SLOTS+1.0*(mud0PktNumB*L_B+mud0CntB)/SLOTS;
        npkts[cur][MUD1_INDEX] = npktsAll[cur][A_MUD1_INDEX]+npktsAll[cur][B_MUD1_INDEX]
        npkts[cur][MUD2_INDEX] = npktsAll[cur][A_MUD2_INDEX]+npktsAll[cur][B_MUD2_INDEX]
        npkts[cur][NCMA_TES1_INDEX] = npktsAll[cur][A_NCMA_TES1_INDEX]+npktsAll[cur][B_NCMA_TES1_INDEX]
        npkts[cur][NCMA_UES1_INDEX] = npktsAll[cur][A_NCMA_UES1_INDEX]+npktsAll[cur][B_NCMA_UES1_INDEX]
        npkts[cur][NCMA_TES2_INDEX] = npktsAll[cur][A_NCMA_TES2_INDEX]+npktsAll[cur][B_NCMA_TES2_INDEX]
        npkts[cur][NCMA_UES2_INDEX] = npktsAll[cur][A_NCMA_UES2_INDEX]+npktsAll[cur][B_NCMA_UES2_INDEX]
    
        [two_equations_num,one_equations_num]=mac_count_equations(phyRawMap)
        npkts[cur][NCMA_UES1_UPPER_INDEX] = 1.0*(two_equations_num*2+one_equations_num)/SLOTS
        [two_equations_num,one_equations_num]=mac_count_equations(phyTwoPhaseMap)
        npkts[cur][NCMA_UES2_UPPER_INDEX] = 1.0*(two_equations_num*2+one_equations_num)/SLOTS
    
        print "ncma-mud1: pktinfo=",[mud1PktNumA,mud1CntA,mud1PktNumB,mud1CntB]
        print "ncma-tes1: pktinfo=",[ncmaT1PktNumA,ncmaT1CntA,ncmaT1PktNumB,ncmaT1CntB]
        print "ncma-ues1: pktinfo=",[ncmaU1PktNumA,ncmaU1CntA,ncmaU1PktNumB,ncmaU1CntB]
        print "ncma-mud2: pktinfo=",[mud2PktNumA,mud2CntA,mud2PktNumB,mud2CntB]
        print "ncma-tes2: pktinfo=",[ncmaT2PktNumA,ncmaT2CntA,ncmaT2PktNumB,ncmaT2CntB]
        print "ncma-ues2: pktinfo=",[ncmaU2PktNumA,ncmaU2CntA,ncmaU2PktNumB,ncmaU2CntB]
        print "round: %d %d %d" % (cur, rr, ll)
        print "npktsAll = [%6f %6f %6f %6f %6f %6f %6f %6f %6f %6f %6f %6f]" % (
                                                                npktsAll[cur][A_MUD1_INDEX], npktsAll[cur][B_MUD1_INDEX],
                                                                npktsAll[cur][A_NCMA_TES1_INDEX], npktsAll[cur][B_NCMA_TES1_INDEX],
		  											            npktsAll[cur][A_NCMA_UES1_INDEX], npktsAll[cur][B_NCMA_UES1_INDEX],
													            npktsAll[cur][A_MUD2_INDEX], npktsAll[cur][B_MUD2_INDEX],
                                                                npktsAll[cur][A_NCMA_TES2_INDEX], npktsAll[cur][B_NCMA_TES2_INDEX],
    										                    npktsAll[cur][A_NCMA_UES2_INDEX], npktsAll[cur][B_NCMA_UES2_INDEX])
        print "npkts = [%6f %6f %6f %6f %6f %6f %6f %6f]" % (npkts[cur][MUD1_INDEX], npkts[cur][NCMA_TES1_INDEX], npkts[cur][NCMA_UES1_INDEX], npkts[cur][NCMA_UES1_UPPER_INDEX], 
                                                             npkts[cur][MUD2_INDEX], npkts[cur][NCMA_TES2_INDEX], npkts[cur][NCMA_UES2_INDEX], npkts[cur][NCMA_UES2_UPPER_INDEX])
        log = "%d %d %d %f %f %f %f %f %f %f %f %f\n"   % (ff, rr, ll, npkts[cur][MUD0_INDEX],
												     npkts[cur][MUD1_INDEX], npkts[cur][NCMA_TES1_INDEX], npkts[cur][NCMA_UES1_INDEX], npkts[cur][NCMA_UES1_UPPER_INDEX], 
    									             npkts[cur][MUD2_INDEX], npkts[cur][NCMA_TES2_INDEX], npkts[cur][NCMA_UES2_INDEX], npkts[cur][NCMA_UES2_UPPER_INDEX])
        fwriter.write(log)
      
fwriter.close()