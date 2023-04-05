#/usr/bin/python
#
# Calculate the probability that the coupled matrix is full rank
# Req: set [t,LA,LB,lA,lB,lX] manually 
#

t=256
LA=15; LB=10;
lA=9; lB=9; lX=7;
#lA=5; lB=5; lX=10;

r=1.0

for ii in range(LA+1-lA,LA+1):
  r = r*(1-1.0/t**ii)

for ii in range(LB+1-lB,LB+1):
  r = r*(1-1.0/t**ii)

for ii in range(LA+LB+1-lA-lB-lX,LA+LB-lA-lB+1):
  r = r*(1-1.0/t**ii)

print r 
