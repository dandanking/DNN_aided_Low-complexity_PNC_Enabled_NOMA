#/usr/bin/python
#
# Author: lzyou@ie.cuhk.edu.hk
# Time: March 13, 2014
#
# Note: test (LA+LB,LA+LB) NCMA matrix, RS or Random code
#

import mu_rs_code, random
N = 255

def get_parameters(LA,LB):
    lA = random.randint(0,LA-1)
    lB = random.randint(0,LB-1)
    lX = LA+LB-(lA+lB)
    return [lA,lB,lX]

def get_parameters2(index,LA,LB):
    index += 1
    if index == 1:
        lA=1; lB=9; lX=10;
    elif index == 2:
        lA=5; lB=5; lX=10;
    elif index == 3:
        lA=5; lB=9; lX=6;
    elif index == 4:
        lA=1; lB=1; lX=18;
    elif index == 5:
        lA=9; lB=9; lX=2;
    return [lA,lB,lX]

def get_parameters3(index):
    index += 1
    assert(index<=5)
    if index == 1:
        LA=16; LB=16; lA=8; lB=8; lX=16;
    elif index == 3:
        LA=16; LB=16; lA=4; lB=12; lX=16;
    elif index == 2:
        LA=16; LB=16; lA=13; lB=15; lX=4;
    elif index == 4:
        LA=24; LB=16; lA=13; lB=15; lX=12;
    elif index == 5:
        LA=24; LB=16; lA=4; lB=12; lX=24;
    return [LA,LB,lA,lB,lX]

#lA = 6;  lB = 6; lX = 8

nmatrix = 1000000
nparas = 5
ntimes = 1

for j in range(0,nparas):
    #[lA,lB,lX] = get_parameters2(j,LA,LB)
    [LA,LB,lA,lB,lX] = get_parameters3(j)
    #lA = 9; lB = 9; lX = 3;
    C = mu_rs_code.NCMAUESCode(N,1,LA,LB,systematic=0)
    print "N=%d LA=%d LB=%d lA=%d lB=%d lX=%d" % (N,LA,LB,lA,lB,lX)
    
    rs = []
    rd = []
    for p in range(ntimes):
        errNum = 0
        #"""
        for i in range(nmatrix):
            m = C.TestMURsMatrix(lA,lB,lX)
            #c = m.Copy()
            #print m
            r2 = m.Rank()
            
            if r2 < LA+LB:
                errNum += 1
                #assert(m.Determinant() == 0)
                
            if r2 == LA+LB:
                #x = c.Inverse()
                #print m
                #print x.rows, x.cols
                #print c.rows, c.cols
                pass
            
        print "RS=%3d/%3d" % (errNum,nmatrix),
        rs.append(1-errNum*1.0/nmatrix)
        #"""
        #"""
        errNum = 0
        for i in range(nmatrix):
            m = C.TestMURdMatrix(lA,lB,lX)
            r = m.Rank()
            if r < LA+LB:
                errNum += 1
                
        print "Random=%3d/%3d" % (errNum,nmatrix)
        rd.append(1-errNum*1.0/nmatrix)
        #"""
        
    print 'RS:',rs
    print 'Random:',rd

