
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

class NCMACode:
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
                fstr = 'A'
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
                fstr = 'B'
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
    
    def SetMUMatrix(self, matrix, tmap, Type, cnt):
        self.maps = tmap
        self.decoderMatrix = genericmatrix.GenericMatrix(
            (0,self.ka+self.kb),0,1,self.field.Add,self.field.Subtract,
            self.field.Multiply,self.field.Divide,fillMode='e')
        
        for i in range(len(matrix)):
            self.decoderMatrix.AppendRow(matrix[i])
            
        assert Type == 'A' or Type == 'B' or Type == 'X'
        if Type == 'A':
            self.pktCntA=cnt
            self.pktCntB=0
            self.pktCntX=0
        elif Type == 'B':
            self.pktCntA=0
            self.pktCntB=cnt
            self.pktCntX=0
        else:    
            self.pktCntA = cnt[0]
            self.pktCntB = cnt[1]

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
        while G.rows > self.ka:
            G.DeleteLastRow()
        Inv = G.Inverse()
        Gy = [random.randint(0, (1<<self.field.n)-1) for ii in range(length)]
        for i in range(self.K):
            decVec = Inv.LeftMulColumnVec(Gy)
    
    def UpdateMatrix(self,TYPE):
        assert TYPE == 'A' or TYPE == 'B' or TYPE == 'X'
        
        matrix = []
        tmap = []
        cnt = 0
        if TYPE == 'A':
            self.XOR_STATE = False
            self.msgCntA += 1     
            G = genericmatrix.GenericMatrix(
                (0,self.ka),0,1,self.field.Add,self.field.Subtract,
                self.field.Multiply,self.field.Divide,fillMode='e')                   
            for i in range(len(self.maps)):
                x = self.decoderMatrix.GetRow(i)
                if self.maps[i] == 'AB' or self.maps[i] == 'A':
                    G.AppendRow(x[0:self.ka])
                if self.maps[i] == 'AB' or self.maps[i] == 'B' or self.maps[i] == 'X':
                    tmap.append('B')
                    cnt += 1
                    matrix.append(x)
                    self.SimulateRecovery(self.ka)  #FIXME
            self.SetMUMatrix(matrix,tmap,'B',cnt)
            self.SimulateCalculation(G,self.ka)  #FIXME
        elif TYPE == 'B':
            self.XOR_STATE = False
            self.msgCntB += 1
            G = genericmatrix.GenericMatrix(
                (0,self.kb),0,1,self.field.Add,self.field.Subtract,
                self.field.Multiply,self.field.Divide,fillMode='e')
            for i in range(len(self.maps)):
                x = self.decoderMatrix.GetRow(i)
                if self.maps[i] == 'AB' or self.maps[i] == 'B':
                    G.AppendRow(x[self.ka:self.ka+self.kb])
                if self.maps[i] == 'AB' or self.maps[i] == 'A' or self.maps[i] == 'X':
                    tmap.append('A')
                    cnt += 1
                    matrix.append(x)
                    self.SimulateRecovery(self.kb)  #FIXME
            self.SetMUMatrix(matrix,tmap,'A',cnt)
            self.SimulateCalculation(G,self.kb)  #FIXME
        else:
            # FIXME: the current implementation is not optimal
            # We should store everything, even lost pkts, and recover them by XOR equation
            # Our current implementation only cares about received pkts
            self.XOR_STATE = True
            cntA = 0; cntB = 0;
            G = genericmatrix.GenericMatrix(
                (0,self.kx),0,1,self.field.Add,self.field.Subtract,
                self.field.Multiply,self.field.Divide,fillMode='e')             
            for i in range(len(self.maps)):
                x = self.decoderMatrix.GetRow(i)
                if self.maps[i] == 'AB' or self.maps[i] == 'A' or self.maps[i] == 'B':
                    cntA += 1
                    cntB += 1
                    tmap.append('AB')
                    matrix.append(x)
                    self.SimulateRecovery(self.kx)  #FIXME
                
                if self.maps[i] == 'AB' or self.maps[i] == 'X':
                    #print x
                    if self.kx == self.ka:
                        y = x[0:self.ka]
                    else:
                        y = x[self.ka:self.ka+self.kb]
                    #if G.rows < self.kx:
                    G.AppendRow(y)
                    
            cnt = []; cnt.append(cntA); cnt.append(cntB)
            self.SetMUMatrix(matrix,tmap,'X',cnt)
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
        
def _test():
    N=255; LA=24; LB=16; K=100; 
    C = NCMACode(N,K,LA,LB,systematic=0)
    #C.RandomTest3(1,0.1,0.1,0.7,0)
    C.RandomTest(100)

if __name__ == "__main__":
    #_test_mu_rs_code()
    #_test_mu_rd_code()
    #_test_su_rd_code()
    _test()
