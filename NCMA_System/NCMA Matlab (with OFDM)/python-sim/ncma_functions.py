#!/usr/bin/python

def phy_get_maps(ss):
  postfix = ['_raw_map.dat','_bridge_map.dat','_twophase_map.dat']
  ret = [[],[],[]]
  for ii in range(len(postfix)):
    ss1=ss+postfix[ii]
    fp=open(ss1)
    while 1:
      r=fp.readline()
      if not r:
        break
      r=r.strip('\n')
      l=[int(n) for n in r.split(',')]
      ret[ii].append(l)
    fp.close()
  return ret

def mac_find_last_index(macMap,cur_index,pktAnum,pktBnum,pktXnum,SELF_RECOVER_SIGNAL,XOR_RECOVER_SIGNAL):
  for pkt_index in range(cur_index,-1,-1):
    [seqA,seqB,seqX] = macMap[pkt_index]
    if seqA == pktAnum or seqA == pktAnum * SELF_RECOVER_SIGNAL or seqA == pktAnum * XOR_RECOVER_SIGNAL:
      break
    if seqB == pktBnum or seqB == pktBnum * SELF_RECOVER_SIGNAL or seqX == pktXnum * XOR_RECOVER_SIGNAL:
      break
    if seqX == pktXnum or seqX == pktXnum * SELF_RECOVER_SIGNAL or seqX == pktXnum * XOR_RECOVER_SIGNAL:
      break
  return pkt_index

def mac_count_equations(phyMap):
  [A_INDEX,B_INDEX,X_INDEX]=range(3)
  two = 0
  one = 0
  for ii in range(len(phyMap)):
    esum = phyMap[ii][A_INDEX]+phyMap[ii][B_INDEX]+phyMap[ii][X_INDEX]
    if esum >= 2: 
      two += 1
    if esum == 1:
      one += 1
  return [two,one]
      
def mac_update_seqno(cur_index,inputSeqInfo,matrix,d_X_decodable,macMap,L_A,L_B,ncma_mode):
  [A_INDEX,B_INDEX,X_INDEX]=range(3)
  [countA,pktAnum,countB,pktBnum,countX,pktXnum]=inputSeqInfo
  NULL_SIGNAL=0; SELF_RECOVER_SIGNAL=-1; XOR_RECOVER_SIGNAL=-10;
  L_X = max(L_A,L_B)

  # XOR is recoverable
  if macMap[cur_index][X_INDEX]==0 and d_X_decodable:
    macMap[cur_index][X_INDEX] = SELF_RECOVER_SIGNAL * pktXnum

  d_updateA_flag = 0; d_updateB_flag = 0; d_count = 0;
  while (countA >= L_A) or (countB >= L_B) or ((d_X_decodable==0) and (countX >= L_X)):
    print "%s: d_count=%d curIndex=%3d countA=%3d countB=%3d countX=%3d d_X_decodable=%d" % (ncma_mode,d_count,cur_index,countA,countB,countX,d_X_decodable)
    last_index = mac_find_last_index(macMap,cur_index,pktAnum-1,pktBnum-1,pktXnum-1,SELF_RECOVER_SIGNAL,XOR_RECOVER_SIGNAL)
    #print "last_index=", last_index

    d_count += 1
    assert(d_count <= 3)

    # check user A
    if countA >= L_A:
      d_updateA_flag = 1
      
      if d_updateB_flag == 0:
        if pktAnum == 1:
          last_index = -1

        # MAC-layer Bridging
        for pkt_index in range(last_index+1,cur_index+1):
          if macMap[pkt_index][A_INDEX] == 0:
            macMap[pkt_index][A_INDEX] = SELF_RECOVER_SIGNAL * pktAnum
          if (macMap[pkt_index][X_INDEX] or d_X_decodable) and (macMap[pkt_index][B_INDEX]==0):
            countB += 1
            macMap[pkt_index][B_INDEX] = XOR_RECOVER_SIGNAL * pktBnum
        
        # update X's state
        d_X_decodable = 0
        pktXnum += 1
        countX = 0

        # update Matrix's state (update Eqn A)
        matrix.UpdateOnEqnSolved(A_INDEX)

      # update self's state
      pktAnum += 1
      countA = 0

    # check user B
    if countB >= L_B:
      d_updateB_flag = 1

      if d_updateA_flag == 0:
        if pktBnum == 1:
          last_index = -1

        # MAC-layer Bridging
        for pkt_index in range(last_index+1,cur_index+1):
          if macMap[pkt_index][B_INDEX] == 0:
            macMap[pkt_index][B_INDEX] = SELF_RECOVER_SIGNAL * pktBnum
          if (d_X_decodable or macMap[pkt_index][X_INDEX]) and (macMap[pkt_index][A_INDEX]==0):
            countA += 1
            macMap[pkt_index][A_INDEX] = XOR_RECOVER_SIGNAL * pktAnum

        # update X's state
        d_X_decodable = 0
        pktXnum += 1
        countX = 0

        # update Matrix's state (update Eqn B)
        matrix.UpdateOnEqnSolved(B_INDEX)

      # update self's state
      pktBnum += 1
      countB = 0
      
    if d_X_decodable == 0 and countX >= L_X:
      #print ncma_mode, cur_index, countX, L_X, d_X_decodable
      d_X_decodable = 1
      if d_updateA_flag == 0 and d_updateB_flag == 0:
        if pktXnum == 1:
          last_index = -1 
        
        # MAC-layer Bridging
        #print last_index,cur_index
        for pkt_index in range(last_index+1,cur_index+1):
          if macMap[pkt_index][X_INDEX] == 0:
            macMap[pkt_index][X_INDEX] = SELF_RECOVER_SIGNAL * pktXnum
            if macMap[pkt_index][B_INDEX] and (macMap[pkt_index][A_INDEX] == 0):
              countA += 1
              macMap[pkt_index][A_INDEX] = XOR_RECOVER_SIGNAL * pktAnum
              matrix.UpdateMatrix(A_INDEX,pkt_index) # update Eqn A
            elif macMap[pkt_index][A_INDEX] and (macMap[pkt_index][B_INDEX] == 0):
              countB += 1
              macMap[pkt_index][B_INDEX] = XOR_RECOVER_SIGNAL * pktBnum
              matrix.UpdateMatrix(B_INDEX,pkt_index) # update Eqn B

  # Check the the 4-th equation
  if matrix.Solveable(pkt_index):
    null = 1

  ret = [countA,pktAnum,countB,pktBnum,countX,pktXnum,d_X_decodable]
  return ret

