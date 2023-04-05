import math
import matplotlib.pyplot as plt
import numpy as np

x=np.linspace(0,360,73,endpoint=True)
ym=[]
yp=[]
dPNC=1
snr_lin1=2.511886e-01
snr_lin2=1.584893e-01
snr_lin3=1.000000e-01
ym1=[]
yp1=[]
ym3=[]
yp3=[]

def qfunc(arg):
    return 0.5*math.erfc(arg/math.sqrt(2))

for ph in range(0,181,5):
    p = 0
    p = 8*qfunc(dPNC / math.sqrt(2 * snr_lin2))
    yp.append((1/4)*p)
yp=yp+yp
del yp[37]
for ph in range(0,181,5):
    p = 0
    dmud1=2*math.sin(ph*math.pi/360)
    dmud2=2*math.sin((180-ph)*math.pi/360)
    p = 2*qfunc(dmud1 / math.sqrt(2 * snr_lin2))+2*qfunc(dmud2 / math.sqrt(2 * snr_lin2))
    ym.append((1/2)*p)
ym=ym+ym
del ym[37]
for ph in range(0,181,5):
    p = 0
    p = 8*qfunc(dPNC / math.sqrt(2 * snr_lin1))
    yp1.append((1/4)*p)
yp1=yp1+yp1
del yp1[37]
for ph in range(0,181,5):
    p = 0
    dmud1=2*math.sin(ph*math.pi/360)
    dmud2=2*math.sin((180-ph)*math.pi/360)
    p = 2*qfunc(dmud1 / math.sqrt(2 * snr_lin1))+2*qfunc(dmud2 / math.sqrt(2 * snr_lin1))
    ym1.append((1/2)*p)
ym1=ym1+ym1
del ym1[37]
for ph in range(0,181,5):
    p = 0
    p = 8*qfunc(dPNC / math.sqrt(2 * snr_lin3))
    yp3.append((1/4)*p)
yp3=yp3+yp3
del yp3[37]
for ph in range(0,181,5):
    p = 0
    dmud1=2*math.sin(ph*math.pi/360)
    dmud2=2*math.sin((180-ph)*math.pi/360)
    p = 2*qfunc(dmud1 / math.sqrt(2 * snr_lin3))+2*qfunc(dmud2 / math.sqrt(2 * snr_lin3))
    ym3.append((1/2)*p)
ym3=ym3+ym3
del ym3[37]
plt.xlim(0,360)# 设置横轴的上下限
plt.xticks(np.linspace(0,360,37,endpoint=True))# 设置横轴记号, 是一个 numpy 数组，包含了从  到  等间隔的  个值
plt.ylim(0.001,1)
plt.yscale('log')
plt.grid(linestyle=":")
plt.plot(x, yp1, c='red', marker="o", alpha=0.5, label="PNC SNR=6dB")
plt.plot(x,yp,c='red',marker="|",alpha=0.5,label="PNC SNR=8dB")
plt.plot(x,yp3,c='red',marker="^",alpha=0.5,label="PNC SNR=10dB")
plt.plot(x, ym1, c='blue', marker="o", alpha=0.5, label="MUD SNR=6dB")
plt.plot(x,ym,c='blue',marker="|",alpha=0.5,label="MUD SNR=8dB")
plt.plot(x,ym3,c='blue',marker="^",alpha=0.5,label="MUD SNR=10dB")
plt.legend(loc="lower left")
plt.xlabel('Relative phase offset $\phi$ (degree)')
plt.ylabel('SER')
plt.show()
