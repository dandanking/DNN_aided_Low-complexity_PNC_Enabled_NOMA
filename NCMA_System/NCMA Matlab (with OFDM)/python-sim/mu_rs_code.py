
# Copyright Emin Martinian 2002.  See below for license terms.
# Version Control Info: $Id: rs_code.py,v 1.3 2003/10/28 21:42:56 emin Exp $

"""
This package implements the RSCode class designed to do
Reed-Solomon encoding and (erasure) decoding.  The following
docstrings provide detailed information on various topics.

  RSCode.__doc__   Describes the RSCode class and how to use it.

  license_doc      Describes the license and lack of warranty.

"""

import ffield
import genericmatrix
import math
import random

[A_INDEX,B_INDEX,X_INDEX]=range(3);

class MUDSimpleCode:
    def __init__(self,LA,LB):
        self.ka = LA
        self.kb = LB
        self.msgCntA = 0
        self.msgCntB = 0
        self.pktCntA = 0
        self.pktCntB = 0
        
    def RxPacket(self,okA,okB,okX=None,index=None):  # for compatibility
        if okA and okB:
            self.pktCntA += 1
            self.pktCntB += 1
        elif okA and not okB:
            self.pktCntA += 1
        elif not okA and okB:
            self.pktCntB += 1
        fstr = 'N'
            
        if self.pktCntA == self.ka:
            self.msgCntA += 1
            self.pktCntA = 0
            fstr = 'A'
        if self.pktCntB == self.kb:
            self.msgCntB += 1
            self.pktCntB = 0 
            fstr = 'B'
            
        return [fstr,None]  # for compatibility
    
    def GetState(self):
        return [self.msgCntA, self.pktCntA, self.msgCntB, self.pktCntB]

class MUDCode:
    """
    The MUD class implements a MUD packet reception instance
    """
    def __init__(self,N,K,ka,kb,log2FieldSize=-1,systematic=1,shouldUseLUT=-1):
        if (log2FieldSize < 0):
            log2FieldSize = int(math.ceil(math.log(N)/math.log(2)))
        self.field = ffield.FField(log2FieldSize,useLUT=shouldUseLUT)
        self.n = N
        self.K = K
        self.ka = ka
        self.kb = kb
        self.k = self.ka+self.kb
        self.fieldSize = 1 << log2FieldSize
        
        # TODO: remove these parameters???
        self.curPktIndex = 0  # for getting generator rows. In real implementation, we do not need this parameter.
        self.debug = False    # debug flag
                
        self.CreateEncoderMatrics(systematic)
        self.ResetRxStatus()
        
        self.decoderMatrixA = genericmatrix.GenericMatrix(
            (0,self.ka),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide,fillMode='e')   
        self.decoderMatrixB = genericmatrix.GenericMatrix(
            (0,self.kb),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide,fillMode='e')   
        
    def CreateEncoderMatrics(self,systematic):                 
        self.encoderMatrixA = genericmatrix.GenericMatrix(
            (self.n,self.ka),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide)
        self.encoderMatrixA[0,0] = 1
        shiftA = 3
        for i in range(0,self.n):
            term = 1
            for j in range(0, self.ka):
                index = (i+shiftA) % self.n
                self.encoderMatrixA[index,j] = term
                term = self.field.Multiply(term,i)            

        self.encoderMatrixB = genericmatrix.GenericMatrix(
            (self.n,self.kb),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide)
        self.encoderMatrixB[0,0] = 1
        shiftB = 5
        for i in range(0,self.n):
            term = 1
            for j in range(0, self.kb):
                index = (i+shiftB) % self.n
                self.encoderMatrixB[index,j] = term
                term = self.field.Multiply(term,i)
        
        if (systematic):
            self.encoderMatrixA.Transpose()
            self.encoderMatrixA.LowerGaussianElim()
            self.encoderMatrixA.UpperInverse()
            self.encoderMatrixA.Transpose()

            self.encoderMatrixB.Transpose()
            self.encoderMatrixB.LowerGaussianElim()
            self.encoderMatrixB.UpperInverse()
            self.encoderMatrixB.Transpose()

    def RxPacket(self,okA,okB,okX,index=None):       
        if index is None:
            pkt_index = self.curPktIndex % self.n
        else:
            pkt_index = index % self.n
            
        #print "Index: ", pkt_index
        seq_no = pkt_index
        itemsA = self.encoderMatrixA.GetRow(seq_no)
        itemsB = self.encoderMatrixB.GetRow(seq_no)
        
        if okA and okB:
            self.pktCntA += 1
            self.pktCntB += 1
            self.decoderMatrixA.AppendRow(itemsA)
            self.decoderMatrixB.AppendRow(itemsB)
        elif okA and not okB:
            self.pktCntA += 1
            self.decoderMatrixA.AppendRow(itemsA)
        elif not okA and okB:
            self.pktCntB += 1
            self.decoderMatrixB.AppendRow(itemsB)
            
        fstr = 'N'
        if self.pktCntA == self.ka:
            self.msgCntA += 1
            self.pktCntA = 0
            self.SimulateCalculation(self.decoderMatrixA, self.ka)
            self.decoderMatrixA = genericmatrix.GenericMatrix(
                (0,self.ka),0,1,self.field.Add,self.field.Subtract,
                self.field.Multiply,self.field.Divide,fillMode='e')   
            fstr = 'A'
                   
        if self.pktCntB == self.kb:
            self.msgCntB += 1
            self.pktCntB = 0 
            #print "Matrix: ", self.decoderMatrixB.rows, self.kb
            #print self.decoderMatrixB
            self.SimulateCalculation(self.decoderMatrixB, self.kb)
            self.decoderMatrixB = genericmatrix.GenericMatrix(
                (0,self.kb),0,1,self.field.Add,self.field.Subtract,
                self.field.Multiply,self.field.Divide,fillMode='e')   
            if fstr == 'N':
                fstr = 'B'
            elif fstr == 'A':
                fstr = 'AB'
            
        self.curPktIndex += 1
        return fstr
    
    def SimulateCalculation(self,G,length):
        while G.rows > length:
            G.DeleteLastRow()
        Inv = G.Inverse()
        Gy = [random.randint(0, (1<<self.field.n)-1) for i in range(length)]
        for i in range(self.K):
            decVec = Inv.LeftMulColumnVec(Gy)
    
    def ResetRxStatus(self):
        self.msgCntA = 0      # number of user A messages received so far
        self.msgCntB = 0      # number of user B messages received so far
        self.pktCntA = 0       # number of user A packets within a message received so far
        self.pktCntB = 0       # number of user B packets within a message received so far        
    
    def GetState(self):
        return [self.msgCntA, self.pktCntA, self.msgCntB, self.pktCntB]
    
    def RandomTest(self,numTests,rxProbA=0.1,rxProbB=0.1):
        assert(rxProbA<1)
        assert(rxProbB<1)
        
        for i in range(numTests):
            self.ResetRxStatus()
            
            for j in range(1000):
                r = random.random()
                okA = False; okB = False
                if r < rxProbA:
                    okA = True
                if r < rxProbB:
                    okB = True
                self.RxPacket(okA, okB, j)
            print "[%d] Result: [%d %d %d %d]" % (i, self.msgCntA, self.pktCntA, self.msgCntB, self.pktCntB)
        
    
class NCMAUESCode:
    """
    The NCMA class implements a NCMA (UES) packet reception instance
    """

    def __init__(self,n,K,ka,kb,log2FieldSize=-1,systematic=1,shouldUseLUT=-1):
        """
        Function:   __init__(n,k,log2FieldSize,systematic,shouldUseLUT)
        Purpose:    Create a Reed-Solomon coder for an (n,k) code.
        Notes:      The last parameters, log2FieldSize, systematic
                    and shouldUseLUT are optional.

                    The log2FieldSize parameter 
                    represents the base 2 logarithm of the field size.
                    If it is omitted, the field GF(2^p) is used where
                    p is the smalles integer where 2^p >= n.

                    If systematic is true then a systematic encoder
                    is created (i.e. one where the first k symbols
                    of the encoded result always match the data).

                    If shouldUseLUT = 1 then a lookup table is used for
                    computing finite field multiplies and divides.
                    If shouldUseLUT = 0 then no lookup table is used.
                    If shouldUseLUT = -1 (the default), then the code
                    decides when a lookup table should be used.
        """
        if (log2FieldSize < 0):
            log2FieldSize = int(math.ceil(math.log(n)/math.log(2)))
        self.field = ffield.FField(log2FieldSize,useLUT=shouldUseLUT)
        self.n = n
        self.K = K
        self.ka = ka
        self.kb = kb
        self.k = self.ka+self.kb
        self.fieldSize = 1 << log2FieldSize
                
        # TODO: remove these parameters???
        self.curPktIndex = 0  # for getting generator rows. In real implementation, we do not need this parameter.
        self.rxFullVec = []   # input packet content. In real implementation, we do not need this parameter.
        
        self.rxPktNum = 0     # total number of pkts received so far   
        self.msgCntA = 0      # number of user A messages received so far
        self.msgCntB = 0      # number of user B messages received so far
        self.pktCntA = 0      # number of user A packets within a message received so far
        self.pktCntB = 0      # number of user B packets within a message received so far

        self.maps = []        # map of the received pkts
        self.rxVec = []       # received pkt content
        self.decVec = []      # decoded source pkts

        self.debug = False    # debug flag
        
        self.CreateEncoderMatrics(systematic)
        self.ResetRxStatus()
        
        self.inVecA = range(self.ka); self.inVecB = range(self.kb);
        self.codedVecA = []; self.codedVecB = []; self.codedVecX = [];
        self.InitTxVec()

    def __repr__(self):
        rep = ('<RSCode (n,k) = (' + `self.n` +', ' + `self.k` + ')'
               + '  over GF(2^' + `self.field.n` + ')\n' +
               `self.encoderMatrixA` + '\n' + '>' + `self.encoderMatrixB` + '\n' + '>')
        return rep
     
    def CreateEncoderMatrics(self,systematic):                 
        self.encoderMatrixA = genericmatrix.GenericMatrix(
            (self.n,self.ka),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide)
        self.encoderMatrixA[0,0] = 1
        shiftA = 3
        for i in range(0,self.n):
            term = 1
            for j in range(0, self.ka):
                index = (i+shiftA) % self.n
                self.encoderMatrixA[index,j] = term
                term = self.field.Multiply(term,i)            

        self.encoderMatrixB = genericmatrix.GenericMatrix(
            (self.n,self.kb),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide)
        self.encoderMatrixB[0,0] = 1
        shiftB = 5
        for i in range(0,self.n):
            term = 1
            for j in range(0, self.kb):
                index = (i+shiftB) % self.n
                self.encoderMatrixB[index,j] = term
                term = self.field.Multiply(term,i)
        
        if (systematic):
            self.encoderMatrixA.Transpose()
            self.encoderMatrixA.LowerGaussianElim()
            self.encoderMatrixA.UpperInverse()
            self.encoderMatrixA.Transpose()

            self.encoderMatrixB.Transpose()
            self.encoderMatrixB.LowerGaussianElim()
            self.encoderMatrixB.UpperInverse()
            self.encoderMatrixB.Transpose()
            
        self.encoderMatrix = genericmatrix.GenericMatrix(
        (self.n,self.ka+self.kb),0,1,self.field.Add,self.field.Subtract,
        self.field.Multiply,self.field.Divide)
        
        for i in range(0,self.n):
            termA = self.encoderMatrixA.data[i]
            termB = self.encoderMatrixB.data[i]
            self.encoderMatrix.SetRow(i,termA+termB)
            
        #print "Encoder Matrix:"
        #print self.encoderMatrixA.rows, self.encoderMatrixA.cols
        #print self.encoderMatrixB.rows, self.encoderMatrixB.cols
        #print self.encoderMatrixA
        #print self.encoderMatrixB     
        #print self.encoderMatrix

    def TestRdMatrix(self,L):
        self.randomMatrix = genericmatrix.GenericMatrix(
            (L,L),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide)
        for i in range(L):
            items = [random.randint(0, self.n-1) for ii in range(L)]
            self.randomMatrix.SetRow(i,items)
        return self.randomMatrix

    def TestMURdMatrix(self,lA,lB,lX,sd=None):
        LA = self.ka
        LB = self.kb
        assert lA+lB+lX == LA+LB
        if sd is not None:   
            random.seed(sd)
        else:
            random.seed()
        self.rdMUMatrix = genericmatrix.GenericMatrix(
            (lA+lB+lX,LA+LB),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide)

        emptyA=[0]*LA
        emptyB=[0]*LB
        for i in range(0,lA+lB+lX):
            itemsA = [random.randint(0,self.n-1) for ii in range(LA)]
            itemsB = [random.randint(0,self.n-1) for ii in range(LB)]
            if i<lA:
                row = itemsA+emptyB
            elif i<lA+lX:
                row = itemsA+itemsB
            else:
                row = emptyA+itemsB
            #print row
            self.rdMUMatrix.SetRow(i,row)
        return self.rdMUMatrix

    # Create NCMA Multi-User Erasure Matrix (use for encoding data)    
    def TestMURsMatrix(self,lA,lB,lX,sd=None):
        LA = self.ka
        LB = self.kb
        #assert lA+lB+lX == LA+LB        
        if sd is not None:   
            random.seed(sd)
        else:
            random.seed()
        self.rsMUMatrix = genericmatrix.GenericMatrix(
            (lA+lB+lX,LA+LB),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide)
        
        # use same generator matrix
        rawlistA = self.encoderMatrixA.data
        rawlistB = self.encoderMatrixB.data
        # shuffle the order randomly
        a = random.sample(range(0,len(rawlistA)),lA+lB+lX)
        b = random.sample(range(0,len(rawlistB)),lA+lB+lX)
        #a.sort() # make items in order
    
        emptyA = [0]*LA
        emptyB = [0]*LB
        for i in range(0,lA+lB+lX):
            itemsA = rawlistA[a[i]]
            itemsB = rawlistB[b[i]]
            if i<lA:
                row = itemsA+emptyB
            elif i<lA+lX:
                row = itemsA+itemsB
            else:
                row = emptyA+itemsB
            #print row
            self.rsMUMatrix.SetRow(i,row)
        return self.rsMUMatrix
    
    def ResetRxStatus(self):
        self.InitRxMatrix()
        
        self.curPktIndex = 0             
        self.msgCntA = 0
        self.msgCntB = 0

        self.debug = False        

    
    def SetMUMatrix(self, matrix, tmap, rx=[], fullRx=[]):
        self.maps = tmap
        self.decoderMatrix = genericmatrix.GenericMatrix(
            (0,self.ka+self.kb),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide,fillMode='e')
        
        for i in range(len(matrix)):
            self.decoderMatrix.AppendRow(matrix[i])

        self.rxFullVec = fullRx
        self.rxVec = rx
        self.decVec = []
        self.rxPktNum = len(rx)

    def InitSURxMatrix(self):
        self.decoderMatrixA = genericmatrix.GenericMatrix(
            (self.ka,self.ka),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide,fillMode='i')
        self.decoderMatrixB = genericmatrix.GenericMatrix(
            (self.kb,self.kb),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide,fillMode='i')

    def InitRxMatrix(self):
        self.decoderMatrix = genericmatrix.GenericMatrix(
            (0,self.ka+self.kb),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide,fillMode='e')        
 
        self.pktCntA = 0
        self.pktCntB = 0

        self.maps = []
        self.rxFullVec = []
        self.rxVec = []
        self.decVec = []
        self.rxPktNum = 0

    def RxPacket(self,okA,okB,okX,index=None):
        if index is None:
            pkt_index = self.curPktIndex % self.n
            rxPktVec = [random.randint(0, (1<<self.field.n)-1) for ii in range(3)]
        else:
            pkt_index = index % self.n
            rxPktVec = [self.codedVecA[pkt_index],self.codedVecB[pkt_index],self.codedVecX[pkt_index]]
                 
        if self.debug:
            print "%d %d: okA=%d okB=%d okX=%d PktA=%d cntA=%d PktB=%d cntB=%d" % (self.curPktIndex, index,okA,okB,okX,self.msgCntA,self.pktCntA,self.msgCntB,self.pktCntB)

        seq_no = pkt_index
        items = [0]*(self.ka+self.kb)
        
        if okA and okB:
            items[0:self.ka] = self.encoderMatrixA.GetRow(seq_no)
            self.maps.append('A')
            self.pktCntA += 1
            self.decoderMatrix.AppendRow(items)
            self.rxVec.append(rxPktVec[0])
            self.rxPktNum += 1
            self.rxFullVec.append(rxPktVec)
            
            items = [0]*(self.ka+self.kb)
            items[self.ka:self.ka+self.kb] = self.encoderMatrixB.GetRow(seq_no)
            self.maps.append('B')
            self.pktCntB += 1
            self.decoderMatrix.AppendRow(items)
            self.rxVec.append(rxPktVec[1])
            self.rxPktNum += 1
            self.rxFullVec.append(rxPktVec)
        elif okA and (not okB):
            items[0:self.ka] = self.encoderMatrixA.GetRow(seq_no)
            self.maps.append('A')
            self.pktCntA += 1
            self.decoderMatrix.AppendRow(items)
            self.rxVec.append(rxPktVec[0])
            self.rxPktNum += 1
            self.rxFullVec.append(rxPktVec)
        elif okB and (not okA):
            items[self.ka:self.ka+self.kb] = self.encoderMatrixB.GetRow(seq_no)
            self.maps.append('B')
            self.pktCntB += 1
            self.decoderMatrix.AppendRow(items)
            self.rxVec.append(rxPktVec[1])
            self.rxPktNum += 1
            self.rxFullVec.append(rxPktVec)
        elif (not okA) and (not okB) and okX:
            items[0:self.ka+self.kb] = self.encoderMatrixA.GetRow(seq_no) + self.encoderMatrixB.GetRow(seq_no)
            self.maps.append('X')
            self.decoderMatrix.AppendRow(items)
            self.rxVec.append(rxPktVec[2])
            self.rxPktNum += 1
            self.rxFullVec.append(rxPktVec)        

        fstr = 'N'
        decVec = []
        flag = True
        while flag:
            flag = False
            if self.pktCntA >= self.ka:
                flag = True
                if self.debug:
                    print "A:",self.curPktIndex
                decVec = self.UpdateMatrix('A')
                if fstr == 'N':
                    fstr = 'A'
                elif fstr == 'B':
                    fstr = 'AB'
                if self.debug:
                    print "A:",self.curPktIndex, self.pktCntA, self.pktCntB, self.decoderMatrix.rows, self.decoderMatrix.cols
                
            if self.pktCntB >= self.kb:
                flag = True
                if self.debug:
                    print "B:",self.curPktIndex
                decVec = self.UpdateMatrix('B')
                if fstr == 'N':
                    fstr = 'B'
                elif fstr == 'A':
                    fstr = 'AB'
                if self.debug:
                    print "B:",self.curPktIndex, self.pktCntA, self.pktCntB, self.decoderMatrix.rows, self.decoderMatrix.cols
                
            r = self.FullRank() 
            if r:
                #if self.decoderMatrix.rows >= self.k:
                #    self.debug = True
                if self.debug:
                    print "AB:",self.curPktIndex, self.rxPktNum, self.pktCntA, self.pktCntB, self.decoderMatrix.rows, self.decoderMatrix.cols
                    #print "GeneratorMatrix = ", self.decoderMatrix
                    #print "Rx FullPkt = ", self.rxFullVec, len(self.rxFullVec)
                    #print "Maps: ", self.maps, len(self.maps)
                    #print "Rx Pkt = ", self.rxVec, len(self.rxVec)
    
                decVec = self.UpdateMatrix('AB')
                fstr = 'AB' 

        # at the end, change the state
        self.curPktIndex += 1
        return [fstr,decVec]

    # TODO: remove RxPacket2 and UpdateMatrix2
    # keep them for compatibility
    # Bug: should use while to check conditions (see RxPacket2)
    # It is not used
    def RxPacket2(self,okA,okB,okX,index=None):
        pkt_index = self.curPktIndex
        if self.debug:
            print "%d: okA=%d okB=%d okX=%d PktA=%d cntA=%d PktB=%d cntB=%d" % (pkt_index,okA,okB,okX,self.msgCntA,self.pktCntA,self.msgCntB,self.pktCntB)

        seq_no = pkt_index % (self.n)
        items = [0]*(self.ka+self.kb)
        if okA and okB:
            items[0:self.ka] = self.encoderMatrixA.GetRow(seq_no)
            self.maps.append('A')
            self.pktCntA += 1
            self.decoderMatrix.AppendRow(items)
            items = [0]*(self.ka+self.kb)
            items[self.ka:self.ka+self.kb] = self.encoderMatrixB.GetRow(seq_no)
            self.maps.append('B')
            self.pktCntB += 1
            self.decoderMatrix.AppendRow(items)
        elif okA and (not okB):
            items[0:self.ka] = self.encoderMatrixA.GetRow(seq_no)
            self.maps.append('A')
            self.pktCntA += 1
            self.decoderMatrix.AppendRow(items)
        elif okB and (not okA):
            items[self.ka:self.ka+self.kb] = self.encoderMatrixB.GetRow(seq_no)
            self.maps.append('B')
            self.pktCntB += 1
            self.decoderMatrix.AppendRow(items)
        elif (not okA) and (not okB) and okX:
            items[0:self.ka+self.kb] = self.encoderMatrixA.GetRow(seq_no) + self.encoderMatrixB.GetRow(seq_no)
            self.maps.append('X')
            self.decoderMatrix.AppendRow(items)
        #print self.decoderMatrix
        #print self.decoderMatrix.rows, self.decoderMatrix.cols, self.ka, self.kb
        #print items            

        flag = False
        fstr = ''
        if self.pktCntA >= self.ka:
            if self.debug:
                print "A:",self.curPktIndex
            self.UpdateMatrix('A')
            flag = True
            fstr = 'A'
            
        if self.pktCntB >= self.kb:
            if self.debug:
                print "B:",self.curPktIndex
            self.UpdateMatrix('B')
            flag = True
            fstr = 'B'
        
        r = self.FullRank()  
        if r:
            if self.debug:
                print "AB:",self.curPktIndex, self.pktCntA, self.pktCntB, self.decoderMatrix.rows, self.decoderMatrix.cols
            self.UpdateMatrix('AB')
            flag = True
            fstr = 'AB' 

        # at the end, change the state
        self.curPktIndex += 1
        return [flag,fstr]

    def GetState(self):
        return [self.msgCntA, self.pktCntA, self.msgCntB, self.pktCntB]

    def FullRank(self):
        # if the nrows is less than ncols, it must not be full rank
        if self.decoderMatrix.rows < self.ka+self.kb:
            return False

        r = self.decoderMatrix.Rank()
        #print "FullRank:", self.decoderMatrix.rows, self.ka+self.kb, r
        return r == self.ka+self.kb

    def SimulateRecovery(self,length):
        # FIXME: simulate other user's data
        assert length == self.ka or length == self.kb
        if length == self.ka:
            data = self.inVecA
        else:
            data = self.inVecB
                        
        v = [random.randint(0, (1<<self.field.n)-1) for ii in range(length)]
        d = []
        GT = genericmatrix.GenericMatrix(
                    (0,length),0,1,self.field.Add,self.field.Subtract,
                    self.field.Multiply,self.field.Divide,fillMode='e')
        GT.AppendRow(v)        
        for i in range(self.K):
            d += GT.LeftMulColumnVec(data)

        #xor operation
        for i in range(self.K):
            d[i] ^= 10

    def UpdateMatrix(self,TYPE):
        assert TYPE == 'A' or TYPE == 'B' or TYPE == 'AB'
        matrix = []
        tmap = []
        cnt = 0
        decVec = []
        rx = []
        fullRx = []
        if TYPE == 'A':
            self.msgCntA += 1
            self.pktCntA = 0
            self.pktCntB = 0
            G = genericmatrix.GenericMatrix(
                (0,self.ka),0,1,self.field.Add,self.field.Subtract,
                self.field.Multiply,self.field.Divide,fillMode='e')
            Gy = []
            
            for i in range(self.decoderMatrix.rows):
                x = self.decoderMatrix.GetRow(i)
                if self.maps[i] == 'B' or self.maps[i] == 'X':
                    if self.pktCntB < self.kb:
                        # FIXME
                        fullRx.append(self.rxFullVec[i])
                        y = self.rxFullVec[i][1]
                        if self.maps[i] == 'X':
                            x[0:self.ka] = [0]*self.ka
                            self.SimulateRecovery(self.ka)  # FIXME
                        matrix.append(x)
                        tmap.append('B')
                        rx.append(y)
                        self.pktCntB += 1
                else:
                    G.AppendRow(x[0:self.ka])
                    Gy.append(self.rxFullVec[i][0])
            self.SetMUMatrix(matrix,tmap,rx,fullRx)
            #print "##### Set TYPE A"
            while G.rows > self.ka:
                G.DeleteLastRow()
            Inv = G.Inverse()
            #print "Self:", G.rows, G.cols
            #print "Inv: ", Inv.rows, Inv.cols, len(Gy[0:self.ka])
            for i in range(self.K):
                decVec = Inv.LeftMulColumnVec(Gy[0:self.ka])
                
        elif TYPE == 'B':
            self.msgCntB += 1
            self.pktCntA = 0
            self.pktCntB = 0
            G = genericmatrix.GenericMatrix(
                (0,self.kb),0,1,self.field.Add,self.field.Subtract,
                self.field.Multiply,self.field.Divide,fillMode='e')
            Gy = []
            
            for i in range(self.decoderMatrix.rows):
                x = self.decoderMatrix.GetRow(i)
                if self.maps[i] == 'A' or self.maps[i] == 'X':
                    if self.pktCntA < self.ka:
                        # FIXME
                        fullRx.append(self.rxFullVec[i])
                        y = self.rxFullVec[i][0]
                        if self.maps[i] == 'X':
                            x[self.ka:self.ka+self.kb] = [0]*self.kb
                            self.SimulateRecovery(self.kb)  # FIXME
                        matrix.append(x)
                        tmap.append('A')
                        rx.append(y)
                        self.pktCntA += 1
                else:
                    G.AppendRow(x[self.ka:self.ka+self.kb])
                    Gy.append(self.rxFullVec[i][1])
            self.SetMUMatrix(matrix,tmap,rx,fullRx)
            #print "##### Set TYPE B"
            while G.rows > self.kb:
                G.DeleteLastRow()
            Inv = G.Inverse()
            #print "Self:", G.rows, G.cols
            #print "Inv: ", Inv.rows, Inv.cols, len(Gy[0:self.kb])
            for i in range(self.K):
                decVec = Inv.LeftMulColumnVec(Gy[0:self.kb])
        else:
            #print "Size = ", self.decoderMatrix.rows, self.decoderMatrix.cols
            #print self.decoderMatrix
            
            # Try to decode the message    
            #Inv = self.decoderMatrix.Inverse()
            if self.debug:
                #print "Orginal Matrix:", self.decoderMatrix, self.decoderMatrix.rows, self.decoderMatrix.cols
                #print "GE Matrix: ", geMatrix, geMatrix.rows, geMatrix.cols
                #print "Result = ", result, result.rows, result.cols
                #print "Inverse Matrix:", Inv, Inv.rows, Inv.cols
                pass
                
            #while geMatrix.rows > self.k:
            #    geMatrix.DeleteLastRow()
            #    result.DeleteLastRow()
            #    print "====="
                
            Inv = self.decoderMatrix.Inverse()
            while Inv.rows > self.k:
                Inv.DeleteLastRow()
                
            #print "##### Set TYPE AB"
            #print "Self:", self.decoderMatrix.rows, self.decoderMatrix.cols
            #print "Inv: ", Inv.rows, Inv.cols, len(self.rxVec), self.rxVec[0:self.k]
            for i in range(self.K):
                decVec = Inv.LeftMulColumnVec(self.rxVec[0:self.k])
            
            #self.decoderMatrix = geMatrix.Inverse()
            #r = result.LeftMulColumnVec(self.rxVec[0:self.k])
            #if self.debug:
                #print "Invert Matrix:", self.decoderMatrix
                #print "iV=",self.rxVec[0:self.k]
                #print "oVr=", r
                #pass
            #decVec = self.decoderMatrix.LeftMulColumnVec(r)
            #if self.debug:
                #print "decV1=",decVec
                #print "decV0=",m
            
            self.msgCntA += 1
            self.msgCntB += 1
            self.InitRxMatrix()
            
        return decVec
    
    def UpdateMatrix2(self,TYPE):
        assert TYPE == 'A' or TYPE == 'B' or TYPE == 'AB'
        matrix = []
        tmap = []
        cnt = 0
        if TYPE == 'A':
            self.msgCntA += 1
            for i in range(self.decoderMatrix.rows):
                x = self.decoderMatrix.GetRow(i)
                if self.maps[i] == 'B' or self.maps[i] == 'X':
                    if self.maps[i] == 'X':
                        x[0:self.ka] = [0]*self.ka
                    matrix.append(x)
                    tmap.append('B')
                    cnt += 1
            self.SetMUMatrix(matrix,tmap,'B',cnt)
        elif TYPE == 'B':
            self.msgCntB += 1
            for i in range(self.decoderMatrix.rows):
                x = self.decoderMatrix.GetRow(i)
                if self.maps[i] == 'A' or self.maps[i] == 'X':
                    if self.maps[i] == 'X':
                        x[self.ka:self.ka+self.kb] = [0]*self.kb
                    matrix.append(x)
                    tmap.append('A')
                    cnt += 1
            self.SetMUMatrix(matrix,tmap,'A',cnt)
        else:
            #print "Size = ", self.encoderMatrix.rows, self.encoderMatrix.cols
            #print self.encoderMatrix
            self.msgCntA += 1
            self.msgCntB += 1
            self.InitRxMatrix()
    
    def Encode(self,data):
        """
        Function:       Encode(data)
        Purpose:        Encode a list of length k into length n.
        """
        assert len(data)==self.k, 'Encode: input data must be size k list.'
        
        return self.encoderMatrix.LeftMulColumnVec(data)
    
    def MUEncode(self,dataA,dataB):
        """
        Function:       Encode(data)
        Purpose:        Encode a list of length k into length n.
        """
        assert len(dataA)==self.ka, 'Encode: input data must be size ka list.'
        assert len(dataB)==self.kb, 'Encode: input data must be size kb list.'
        
        A = self.encoderMatrixA.LeftMulColumnVec(dataA)
        B = self.encoderMatrixB.LeftMulColumnVec(dataB)
        X = self.encoderMatrix.LeftMulColumnVec(dataA+dataB)
        return [A,B,X]

    def PrepareDecoder(self,unErasedLocations):
        """
        Function:       PrepareDecoder(erasedTerms)
        Description:    The input unErasedLocations is a list of the first
                        self.k elements of the codeword which were 
                        NOT erased.  For example, if the 0th, 5th,
                        and 7th symbols of a (16,5) code were erased,
                        then PrepareDecoder([1,2,3,4,6]) would
                        properly prepare for decoding.
        """
        if (len(unErasedLocations) != self.k):
            raise ValueError, 'input must be exactly length k'
        
        limitedEncoder = genericmatrix.GenericMatrix(
            (self.k,self.k),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide)
        for i in range(0,self.k):
            limitedEncoder.SetRow(
                i,self.encoderMatrix.GetRow(unErasedLocations[i]))
        self.decoderMatrix = limitedEncoder.Inverse()

    def Decode(self,unErasedTerms):
        """
        Function:       Decode(unErasedTerms)
        Purpose:        Use the
        Description:
        """
        return self.decoderMatrix.LeftMulColumnVec(unErasedTerms)

    def DecodeImmediate(self,data):
        """
        Function:       DecodeImmediate(data)
        Description:    Takes as input a data vector of length self.n
                        where erased symbols are set to None and
                        returns the decoded result provided that
                        at least self.k symbols are not None.

                        For example, for an (n,k) = (6,4) code, a
                        decodable input vector would be
                        [2, 0, None, 1, 2, None].
        """

        if (len(data) != self.n):
            raise ValueError, 'input must be a length n list'

        unErasedLocations = []
        unErasedTerms = []
        for i in range(self.n):
            if (None != data[i]):
                unErasedLocations.append(i)
                unErasedTerms.append(data[i])
        self.PrepareDecoder(unErasedLocations[0:self.k])
        return self.Decode(unErasedTerms[0:self.k])
    
    def DecodeImmediate2(self,data):
        if (len(data) != self.n):
            raise ValueError, 'input must be a length n list'

        unErasedLocations = []
        unErasedTerms = []
        for i in range(self.n):
            if (None != data[i]):
                unErasedLocations.append(i)
                unErasedTerms.append(data[i])
        self.PrepareDecoder(unErasedLocations)
        return self.Decode(unErasedTerms[0:self.k])

    def PrepareDecoder2(self,unErasedLocations):
        if (len(unErasedLocations) < self.k):
            raise ValueError, 'input must be larger than length k'
        
        limitedEncoder = genericmatrix.GenericMatrix(
            (1,self.k),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide)
        for i in range(0,len(unErasedLocations)):
            limitedEncoder.AppendRow(
                self.encoderMatrix.GetRow(unErasedLocations[i]))
            if limitedEncoder.Rank():
                break
        self.decoderMatrix = limitedEncoder.Inverse()
        
    def InitTxVec(self):
        # encoding
        random.seed(100)
        for j in range(self.ka):
            self.inVecA[j] = random.randint(0, (1<<self.field.n)-1)
        for j in range(self.kb):
            self.inVecB[j] = random.randint(0, (1<<self.field.n)-1)
        [self.codedVecA,self.codedVecB,self.codedVecX] = self.MUEncode(self.inVecA,self.inVecB)
        
    def RandomTest3(self,numTests,rxProbA=0.1,rxProbB=0.1,rxProbX=0.1,rxProbAll=0.2):
        # Simulate the decoding procedure packet by packet
        assert(rxProbA+rxProbB+rxProbX+rxProbAll<=1)
        import random

        for i in range(numTests):
            #print "************** %d ******************" % (i)
            self.ResetRxStatus()

            # erasure and decoding
            #random.seed(1)
            for j in range(1000): #self.n):
                # erasure
                r = random.random()
                #index = j % N
                if r < rxProbA:
                    okA = True; okB = False; okX = False;
                elif r < rxProbA+rxProbB:
                    okA = False; okB = True; okX = False;
                elif r < rxProbA+rxProbB+rxProbX:
                    okA = False; okB = False; okX = True;
                elif r < rxProbA+rxProbB+rxProbX+rxProbAll:
                    okA = True; okB = True; okX = True;
                else:
                    okA = False; okB = False; okX = False;
                    continue;
                    
                [s,decVec] = self.RxPacket(okA,okB,okX,j)
                
                # not decodable
                #if len(s) == 0:
                #    continue

                #if s=='AB':
                #    inVec = self.inVecA + self.inVecB
                #if s=='A':
                #    inVec = self.inVecA
                #if s=='B':
                #    inVec = self.inVecB

                #assert decVec == inVec, ('inVec = ' + `inVec`
                #                + '\ndecVec = ' + `decVec`)
            
            print "[%d] Result: [%d %d %d %d]" % (i, self.msgCntA, self.pktCntA, self.msgCntB, self.pktCntB)

    def RandomTest2(self,numTests):
        # This test takes the whole message as a block
        # It is not the real practice 
        maxErasures = self.n-self.k
        for i in range(numTests):
            inVec = range(self.k)
            for j in range(self.k):
                inVec[j] = random.randint(0, (1<<self.field.n)-1)
            codedVec = self.Encode(inVec)
            
            numErasures = random.randint(0,maxErasures)
            for j in range(numErasures):
                j = random.randint(0,self.n-1)
                while(codedVec[j] == None):
                    j = random.randint(0,self.n-1)
                codedVec[j] = None
            decVec = self.DecodeImmediate2(codedVec)
            assert decVec == inVec, ('inVec = ' + `inVec`
                                     + '\ncodedVec = ' + `codedVec`
                                     + '\ndecVec = ' + `decVec`)
        
    def RandomTest(self,numTests):
        # This test takes the whole message as a block
        # It is not the real practice
        import random
        
        maxErasures = self.n-self.k
        for i in range(numTests):
            inVec = range(self.k)
            for j in range(self.k):
                inVec[j] = random.randint(0, (1<<self.field.n)-1)
            codedVec = self.Encode(inVec)
            numErasures = random.randint(0,maxErasures)
            for j in range(numErasures):
                j = random.randint(0,self.n-1)
                while(codedVec[j] == None):
                    j = random.randint(0,self.n-1)
                codedVec[j] = None
            decVec = self.DecodeImmediate(codedVec)
            assert decVec == inVec, ('inVec = ' + `inVec`
                                     + '\ncodedVec = ' + `codedVec`
                                     + '\ndecVec = ' + `decVec`)

class NCMATESCode:
    """
    The NCMA class implements a NCMA packet reception instance
    """

    def __init__(self,n,K,ka,kb,log2FieldSize=-1,systematic=1,shouldUseLUT=-1):
        """
        Function:   __init__(n,k,log2FieldSize,systematic,shouldUseLUT)
        Purpose:    Create a Reed-Solomon coder for an (n,k) code.
        Notes:      The last parameters, log2FieldSize, systematic
                    and shouldUseLUT are optional.

                    The log2FieldSize parameter 
                    represents the base 2 logarithm of the field size.
                    If it is omitted, the field GF(2^p) is used where
                    p is the smalles integer where 2^p >= n.

                    If systematic is true then a systematic encoder
                    is created (i.e. one where the first k symbols
                    of the encoded result always match the data).

                    If shouldUseLUT = 1 then a lookup table is used for
                    computing finite field multiplies and divides.
                    If shouldUseLUT = 0 then no lookup table is used.
                    If shouldUseLUT = -1 (the default), then the code
                    decides when a lookup table should be used.
        """
        if (log2FieldSize < 0):
            log2FieldSize = int(math.ceil(math.log(n)/math.log(2)))
        self.field = ffield.FField(log2FieldSize,useLUT=shouldUseLUT)
        self.n = n
        self.K = K
        self.ka = ka
        self.kb = kb
        self.kx = max(self.ka,self.kb)
        self.fieldSize = 1 << log2FieldSize
                
        # TODO: remove these parameters
        self.curPktIndex = 0  # for getting generator rows. In real implementation, we do not need this parameter.
        self.rxFullVec = []   # input packet content. In real implementation, we do not need this parameter.
        
        self.rxPktNum = 0     # total number of pkts received so far   
        self.msgCntA = 0      # number of user A messages received so far
        self.msgCntB = 0      # number of user B messages received so far
        self.pktCntA = 0      # number of user A packets within a message received so far
        self.pktCntB = 0      # number of user B packets within a message received so far
        self.pktCntX = 0

        self.maps = []        # map of the received pkts
        self.rxVec = []       # received pkt content
        self.decVec = []      # decoded source pkts

        self.debug = False    # debug flag
        self.debug2 = False
        self.XOR_STATE= False
        
        self.CreateMatrics(systematic)
        self.ResetRxStatus()
        
        self.inVecA = range(self.ka); self.inVecB = range(self.kb);
        self.codedVecA = []; self.codedVecB = []; self.codedVecX = [];
        self.InitTxVec()
        
    def __repr__(self):
        rep = ('<RSCode (n,k) = (' + `self.n` +', ' + `self.k` + ')'
               + '  over GF(2^' + `self.field.n` + ')\n' +
               `self.encoderMatrixA` + '\n' + '>' + `self.encoderMatrixB` + '\n' + '>')
        return rep
     
    def ResetRxStatus(self):
        self.InitRxMatrix()
        
        self.curPktIndex = 0             
        self.msgCntA = 0
        self.msgCntB = 0

        self.debug = False 
        
    def InitTxVec(self):
        # encoding
        random.seed(100)
        for j in range(self.ka):
            self.inVecA[j] = random.randint(0, (1<<self.field.n)-1)
        for j in range(self.kb):
            self.inVecB[j] = random.randint(0, (1<<self.field.n)-1)
    
    def InitRxMatrix(self):
        self.decoderMatrix = genericmatrix.GenericMatrix(
            (0,self.ka+self.kb),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide,fillMode='e')        
 
        self.pktCntA = 0
        self.pktCntB = 0

        self.maps = []
        self.rxPktNum = 0
         
    def CreateMatrics(self,systematic):                 
        self.encoderMatrixA = genericmatrix.GenericMatrix(
            (self.n,self.ka),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide)
        self.encoderMatrixA[0,0] = 1
        shiftA = 3
        for i in range(0,self.n):
            term = 1
            for j in range(0, self.ka):
                index = (i+shiftA) % self.n
                self.encoderMatrixA[index,j] = term
                term = self.field.Multiply(term,i)            

        self.encoderMatrixB = genericmatrix.GenericMatrix(
            (self.n,self.kb),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide)
        self.encoderMatrixB[0,0] = 1
        shiftB = shiftA
        for i in range(0,self.n):
            term = 1
            for j in range(0, self.kb):
                index = (i+shiftB) % self.n
                self.encoderMatrixB[index,j] = term
                term = self.field.Multiply(term,i)
        
        if (systematic):
            self.encoderMatrixA.Transpose()
            self.encoderMatrixA.LowerGaussianElim()
            self.encoderMatrixA.UpperInverse()
            self.encoderMatrixA.Transpose()

            self.encoderMatrixB.Transpose()
            self.encoderMatrixB.LowerGaussianElim()
            self.encoderMatrixB.UpperInverse()
            self.encoderMatrixB.Transpose()
            
        self.encoderMatrixX = genericmatrix.GenericMatrix(
        (self.n,self.kx),0,1,self.field.Add,self.field.Subtract,
        self.field.Multiply,self.field.Divide)
        
        for i in range(0,self.n):
            if self.kx == self.ka:
                term = self.encoderMatrixA.data[i]
            if self.kx == self.kb:
                term = self.encoderMatrixB.data[i]
            self.encoderMatrixX.SetRow(i,term)
            
        #print "Encoder Matrix:"
        #print self.encoderMatrixA.rows, self.encoderMatrixA.cols
        #print self.encoderMatrixB.rows, self.encoderMatrixB.cols
        #print self.encoderMatrixA
        #print self.encoderMatrixB     
        #print self.encoderMatrix

    # assume bridging results [okA,okB,okX]
    def RxPacket(self,okA,okB,okX,index=None):      
        """
        The TES.MAP_TYPE: A, B, AB, X
        which is different from UES.MAP_TYPE = A, B, X
        """ 
        if index is None:
            pkt_index = self.curPktIndex % self.n
        else:
            pkt_index = index % self.n
        
        # let self.decodeMatrix contain all generators   
        # use map to differentiate        
        seq_no = pkt_index
        items = [0]*(self.ka+self.kb)
        items[0:self.ka] = self.encoderMatrixA.GetRow(seq_no)
        items[self.ka:self.ka+self.kb] = self.encoderMatrixB.GetRow(seq_no)
        
        if not okX and self.XOR_STATE:
            okX = True
            if not okA and okB:
                okA = True
            if okA and not okB:
                okB = True
        
        pkt_index = self.curPktIndex
        if self.debug2:
            print "%d: okA=%d okB=%d okX=%d msgA=%d msgB=%d pktA=%d pktB=%d pktX=%d %s" % (pkt_index,okA,okB,okX,self.msgCntA,self.msgCntB,self.pktCntA,self.pktCntB,self.pktCntX,self.XOR_STATE)
        
        if okA and okB:
            self.pktCntA += 1
            self.pktCntB += 1
            self.pktCntX += 1
            self.maps.append('AB')
            self.decoderMatrix.AppendRow(items)
            #print "+%d: okA=%d okB=%d okX=%d msgA=%d pktA=%d msgB=%d pktB=%d %s cntX=%d" % (pkt_index,okA,okB,okX,self.msgCntA,self.pktCntA,self.msgCntB,self.pktCntB,self.XOR_STATE,self.pktCntX)
        elif okA and not okB:
            self.pktCntA += 1
            self.maps.append('A')
            self.decoderMatrix.AppendRow(items)
        elif okB and not okA:
            self.pktCntB += 1
            self.maps.append('B')
            self.decoderMatrix.AppendRow(items)
        elif (not okA) and (not okB) and okX:
            self.pktCntX += 1
            if 1: #self.XOR_STATE == False:
                self.maps.append('X')
                self.decoderMatrix.AppendRow(items)

        fstr = 'N'
        decVec = []
        flag = True
        while flag:
            flag = False
            #if self.debug2:
            #    print "D-%d: okA=%d okB=%d okX=%d msgA=%d msgB=%d pktA=%d pktB=%d pktX=%d %s" % (pkt_index,okA,okB,okX,self.msgCntA,self.msgCntB,self.pktCntA,self.pktCntB,self.pktCntX,self.XOR_STATE)
                
            if pkt_index == 29:
                pkt_index = 29
            
            if self.pktCntA >= self.ka:
                #print "A:", self.curPktIndex, self.pktCntA, self.pktCntB, self.pktCntX, self.ka, self.kb, self.kx
                
                flag = True
                if self.debug:
                    print "A:",self.curPktIndex
                decVec = self.UpdateMatrix('A')
                if fstr == 'N' or fstr == 'X':
                    fstr = 'A'
                elif fstr == 'B':
                    fstr = 'AB'
                if self.debug:
                    print "A:",self.curPktIndex, self.pktCntA, self.pktCntB, self.decoderMatrix.rows, self.decoderMatrix.cols
                    print "A%d: okA=%d okB=%d okX=%d msgA=%d msgB=%d pktA=%d pktB=%d pktX=%d %s" % (pkt_index,okA,okB,okX,self.msgCntA,self.msgCntB,self.pktCntA,self.pktCntB,self.pktCntX,self.XOR_STATE)
            
            
            if self.pktCntB >= self.kb:
                #print "B:", self.curPktIndex, self.pktCntA, self.pktCntB, self.pktCntX, self.ka, self.kb, self.kx
                
                flag = True
                if self.debug:
                    print "B:",self.curPktIndex
                    print "B%d: okA=%d okB=%d okX=%d msgA=%d msgB=%d pktA=%d pktB=%d pktX=%d %s" % (pkt_index,okA,okB,okX,self.msgCntA,self.msgCntB,self.pktCntA,self.pktCntB,self.pktCntX,self.XOR_STATE)
                    print self.maps, len(self.maps)
                decVec = self.UpdateMatrix('B')
                if fstr == 'N' or fstr == 'X':
                    fstr = 'B'
                elif fstr == 'A':
                    fstr = 'AB'
                if self.debug:
                    print "B:",self.curPktIndex, self.pktCntA, self.pktCntB, self.decoderMatrix.rows, self.decoderMatrix.cols
                    print "B%d: okA=%d okB=%d okX=%d msgA=%d msgB=%d pktA=%d pktB=%d pktX=%d %s" % (pkt_index,okA,okB,okX,self.msgCntA,self.msgCntB,self.pktCntA,self.pktCntB,self.pktCntX,self.XOR_STATE)
                
            if self.pktCntX >= self.kx and self.XOR_STATE == False: 
                #print "X:", self.curPktIndex, self.pktCntA, self.pktCntB, self.pktCntX, self.ka, self.kb, self.kx
                
                flag = True           
                decVec = self.UpdateMatrix('X')
                fstr = 'X' 
                if self.debug:
                    print "X:",self.curPktIndex, self.rxPktNum, self.pktCntA, self.pktCntB, self.decoderMatrix.rows, self.decoderMatrix.cols
                    print "X%d: okA=%d okB=%d okX=%d msgA=%d msgB=%d pktA=%d pktB=%d pktX=%d %s" % (pkt_index,okA,okB,okX,self.msgCntA,self.msgCntB,self.pktCntA,self.pktCntB,self.pktCntX,self.XOR_STATE)
    

        # at the end, change the state
        self.curPktIndex += 1
        return [fstr,decVec]
    
    def SetMUMatrix(self, matrix, tmap):
        self.maps = tmap
        self.decoderMatrix = genericmatrix.GenericMatrix(
            (0,self.ka+self.kb),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide,fillMode='e')
        
        for i in range(len(matrix)):
            self.decoderMatrix.AppendRow(matrix[i])

        self.decVec = []
        
    def SimulateRecovery(self,length):
        # FIXME: simulate other user's data
        assert length == self.ka or length == self.kb
        if length == self.ka:
            data = self.inVecA
        else:
            data = self.inVecB
                        
        v = [random.randint(0, (1<<self.field.n)-1) for ii in range(length)]
        d = []
        GT = genericmatrix.GenericMatrix(
                    (0,length),0,1,self.field.Add,self.field.Subtract,
                    self.field.Multiply,self.field.Divide,fillMode='e')
        GT.AppendRow(v)
        for i in range(self.K):
            d += GT.LeftMulColumnVec(data)
            
        #xor operation
        for i in range(self.K):
            d[i] ^= 10
            
    def SimulateCalculation(self,G,length):
        while G.rows > length:
            G.DeleteLastRow()
        Inv = G.Inverse()
        Gy = [random.randint(0, (1<<self.field.n)-1) for ii in range(length)]
        for i in range(self.K):
            decVec = Inv.LeftMulColumnVec(Gy)
    
    def UpdateMatrix(self,TYPE):
        assert TYPE == 'A' or TYPE == 'B' or TYPE == 'X'
        
        matrix = []
        tmap = []
        if TYPE == 'A':
            self.XOR_STATE = False
            self.msgCntA += 1     
            self.pktCntX = 0; self.pktCntA = 0; self.pktCntB = 0;
            G = genericmatrix.GenericMatrix(
                (0,self.ka),0,1,self.field.Add,self.field.Subtract,
                self.field.Multiply,self.field.Divide,fillMode='e')                   
            for i in range(len(self.maps)):
                x = self.decoderMatrix.GetRow(i)
                if self.maps[i] == 'AB' or self.maps[i] == 'A':
                    G.AppendRow(x[0:self.ka])
                if self.maps[i] == 'AB' or self.maps[i] == 'B' or self.maps[i] == 'X':
                    if self.pktCntB < self.kb:
                        tmap.append('B')
                        self.pktCntB += 1
                        matrix.append(x)
                        if self.maps[i] == 'X':
                            self.SimulateRecovery(self.ka)  #FIXME
            self.SetMUMatrix(matrix,tmap)
            self.SimulateCalculation(G,self.ka)  #FIXME
        elif TYPE == 'B':
            self.XOR_STATE = False
            self.msgCntB += 1
            self.pktCntX = 0; self.pktCntA = 0; self.pktCntB = 0;
            G = genericmatrix.GenericMatrix(
                (0,self.kb),0,1,self.field.Add,self.field.Subtract,
                self.field.Multiply,self.field.Divide,fillMode='e')
            for i in range(len(self.maps)):
                x = self.decoderMatrix.GetRow(i)
                if self.maps[i] == 'AB' or self.maps[i] == 'B':
                    G.AppendRow(x[self.ka:self.ka+self.kb])
                if self.maps[i] == 'AB' or self.maps[i] == 'A' or self.maps[i] == 'X':
                    if self.pktCntA < self.ka:
                        tmap.append('A')
                        self.pktCntA += 1
                        matrix.append(x)
                        if self.maps[i] == 'X':
                            self.SimulateRecovery(self.kb)  #FIXME
            self.SetMUMatrix(matrix,tmap)
            self.SimulateCalculation(G,self.kb)  #FIXME
        else:
            # FIXME: the current implementation is not optimal
            # We should store everything, even lost pkts, and recover them by XOR equation
            # Our current implementation only cares about received pkts
            self.XOR_STATE = True
            self.pktCntA = 0;  self.pktCntB = 0;
            G = genericmatrix.GenericMatrix(
                (0,self.kx),0,1,self.field.Add,self.field.Subtract,
                self.field.Multiply,self.field.Divide,fillMode='e')             
            for i in range(len(self.maps)):
                x = self.decoderMatrix.GetRow(i)
                if self.maps[i] == 'AB' or self.maps[i] == 'A' or self.maps[i] == 'B':
                    if self.pktCntA < self.ka or self.pktCntB < self.kb:    # FIXME: 
                        self.pktCntA += 1;  self.pktCntB += 1
                        tmap.append('AB')
                        matrix.append(x)
                        if self.maps[i] == 'A' or self.maps[i] == 'B':
                            self.SimulateRecovery(self.kx)  #FIXME
                
                if self.maps[i] == 'AB' or self.maps[i] == 'X':
                    #print x
                    if self.kx == self.ka:
                        y = x[0:self.ka]
                    else:
                        y = x[self.ka:self.ka+self.kb]
                    #if G.rows < self.kx:
                    G.AppendRow(y)
                    
            self.SetMUMatrix(matrix,tmap)
            self.SimulateCalculation(G,self.kx)  #FIXME
            
        decVec = self.Decode(TYPE,G)

    def Decode(self,TYPE,G):
        assert TYPE == 'A' or TYPE == 'B' or TYPE == 'X'
        if TYPE == 'A':
            y = [random.randint(0, (1<<self.field.n)-1) for ii in range(self.ka)]
        elif TYPE == 'B':
            y = [random.randint(0, (1<<self.field.n)-1) for ii in range(self.kb)]
        else:
            y = [random.randint(0, (1<<self.field.n)-1) for ii in range(self.kx)]
        #print "G:", TYPE, G.rows, G.cols, len(y)
        #print G
        Inv = G.Inverse()
        #print "Inv:", TYPE, Inv.rows, Inv.cols, len(y)
        for i in range(self.K):
            decVec = Inv.LeftMulColumnVec(y)
        return decVec

    def GetState(self):
        return [self.msgCntA, self.pktCntA, self.msgCntB, self.pktCntB]
        
    def RandomTest(self,numTests,rxProbA=0.1,rxProbB=0.1,rxProbX=0.1,rxProbAll=0.2):
        # Simulate the decoding procedure packet by packet
        assert(rxProbA+rxProbB+rxProbX+rxProbAll<=1)
        import random

        for i in range(numTests):
            #print "************** %d ******************" % (i)
            self.ResetRxStatus()

            # erasure and decoding
            #random.seed(100)
            for j in range(1000):
                # erasure
                r = random.random()
                #index = j % N
                if r < rxProbA:
                    okA = True; okB = False; okX = False;
                elif r < rxProbA+rxProbB:
                    okA = False; okB = True; okX = False;
                elif r < rxProbA+rxProbB+rxProbX:
                    okA = False; okB = False; okX = True;
                elif r < rxProbA+rxProbB+rxProbX+rxProbAll:
                    okA = True; okB = True; okX = True;
                else:
                    okA = False; okB = False; okX = False;
                    
                self.RxPacket(okA,okB,okX,j)
                            
            print "[%d] Result: [%d %d %d %d]" % (i, self.msgCntA, self.pktCntA, self.msgCntB, self.pktCntB)

license_doc = """
  This code was originally written by Emin Martinian (emin@allegro.mit.edu).
  You may copy, modify, redistribute in source or binary form as long
  as credit is given to the original author.  Specifically, please
  include some kind of comment or docstring saying that Emin Martinian
  was one of the original authors.  Also, if you publish anything based
  on this work, it would be nice to cite the original author and any
  other contributers.

  There is NO WARRANTY for this software just as there is no warranty
  for GNU software (although this is not GNU software).  Specifically
  we adopt the same policy towards warranties as the GNU project:

  BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE PROGRAM 'AS IS' WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
REPAIR OR CORRECTION.

  IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.
"""


# The following code is used to make the doctest package
# check examples in docstrings.

def get_parameters(LA,LB):
    lA = random.randint(0,LA-1)
    lB = random.randint(0,LB-1)
    lX = LA+LB-(lA+lB)
    return [lA,lB,lX]

def _test_su_rd_code():
    N=255; L=10; T=1#0000
    errNum = 0
    C = NCMAUESCode(N,1,1)
    for i in range(T):
        m = C.TestRdMatrix(L)
        c = m.Copy()
        r = m.Rank()
        if r<L:
            errNum += 1
            #print errNum, i, r
    print errNum,T
        
def _test_mu_rd_code():
    #import doctest, mu_rs_code
    #return doctest.testmod(mu_rs_code)
    N=255; LA=10; LB=10; T=100000;
    C = NCMAUESCode(N,LA,LB,systematic=0)

    #[lA,lB,lX] = get_parameters(LA,LB)
    lA=8; lB=8; lX=4;
    errNum = 0    
    for i in range(T):    
        m = C.TestMURdMatrix(lA,lB,lX)
        r = m.Rank()
        if r < LA+LB:
            errNum += 1
            #print errNum, i, r
    print errNum,T

def _test_mu_rs_code():
    #import doctest, mu_rs_code
    #return doctest.testmod(mu_rs_code)
    N=255; LA=10; LB=10;
    C = NCMAUESCode(N,LA,LB,systematic=0)
    for i in range(10000000):
        #[lA,lB,lX] = get_parameters(LA,LB)
        lA=0; lB=0; lX=20;
        
        m = C.TestMURsMatrix(lA,lB,lX)
        c = m.Copy()
        [r1,] = m.FullRank()  # call genericmatrix's function
        r2 = m.Rank()
        if (r1 and r2 < LA+LB) or (not r1 and r2 == LA+LB):
            print "Error: two rank functions are not consistent"
            print "N=%d LA=%d LB=%d lA=%d lB=%d lX=%d" % (N,LA,LB,lA,lB,lX)
            print c, m, r1, r2
            break

def _test_mud():
    N=255; LA=24; LB=16; K=100;
    C = MUDCode(N,K,LA,LB,systematic=0)
    C.RandomTest(10, 0.5, 0.5)

def _test_tes():
    N=255; LA=24; LB=16; K=100; 
    C = NCMATESCode(N,K,LA,LB,systematic=0)
    #C.RandomTest3(1,0.1,0.1,0.7,0)
    C.RandomTest(100)
    
def _test_ues():
    N=255; LA=16; LB=16; K=100; 
    C = NCMAUESCode(N,K,LA,LB,systematic=0)
    #C.RandomTest3(1,0.1,0.1,0.7,0)
    C.RandomTest3(100)
    
if __name__ == "__main__":
    #_test_mu_rs_code()
    #_test_mu_rd_code()
    #_test_su_rd_code()
    #_test_mud()
    _test_tes()
    #_test_ues()
