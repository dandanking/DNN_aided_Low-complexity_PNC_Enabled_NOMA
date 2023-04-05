#/usr/bin/python
#
# Author: lzyou@ie.cuhk.edu.hk
# Time: March 13, 2014
#
# Note: test (LA+LB,LA+LB) NCMA matrix, RS or Random code
#

import mu_rs_code, random
N = 255
LA = 15
LB = 10
C = mu_rs_code.NCMACode(N,LA,LB,systematic=0)

def get_parameters(LA,LB):
  lA = random.randint(0,LA-1)
  lB = random.randint(0,LB-1)
  lX = LA+LB-(lA+lB)
  return [lA,lB,lX]

def get_parameters2(index,LA,LB):
  index += 1
  if index == 1:
    lA=1; lB=9; lX=10+5;
  elif index == 2:
    lA=5; lB=5; lX=10+5;
  elif index == 3:
    lA=5; lB=9; lX=6+5;
  elif index == 4:
    lA=1; lB=1; lX=18+5;
  elif index == 5:
    lA=9; lB=9; lX=2+5;
  return [lA,lB,lX]

#lA = 6;  lB = 6; lX = 8

nmatrix = 1000000
nparas = 5
ntimes = 1

for j in range(1,nparas):
  [lA,lB,lX] = get_parameters2(j,LA,LB)
  #lA = 8; lB = 9; lX = 3;
  print "N=%d LA=%d LB=%d lA=%d lB=%d lX=%d" % (N,LA,LB,lA,lB,lX)
  
  rs = []
  rd = []
  for p in range(ntimes):
    errNum = 0
    for i in range(nmatrix):
      m = C.TestMURsMatrix(lA,lB,lX)
      c = m.Copy()
      r2 = m.Rank2()

      if r2 < LA+LB:
        errNum += 1
        #print c, m, r2
        #break
        #assert(m.Determinant() == 0)
      
    print "RS=%3d/%3d" % (errNum,nmatrix),
    rs.append(1-errNum*1.0/nmatrix)
            
    errNum = 0
    for i in range(nmatrix):
        m = C.TestMURdMatrix(lA,lB,lX)
        r = m.Rank2()
        if r < LA+LB:
            errNum += 1
        
    print "Random=%3d/%3d" % (errNum,nmatrix)
    rd.append(1-errNum*1.0/nmatrix)


  print 'RS:',rs
  print 'Random:',rd

