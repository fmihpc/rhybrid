import subprocess as sp
import os
import sys
import numpy as np
import matplotlib.pyplot as plt


# check run finished with no errors or warnings
if os.path.isfile('logfile.txt') == False:
 print('FAILED: logfile.txt not found')
 quit()
res = sp.run(['grep -q \'(MAIN) Exiting simulation after successful run.\' logfile.txt'],shell=True)
if res.returncode != 0:
 print('FAILED: run did not finish')
 quit()
res = sp.run(['grep -qi \'error\' logfile.txt'],shell=True)
if res.returncode == 0:
 print('FAILED: run had errors')
 quit()
res = sp.run(['grep -qi \'warning\' logfile.txt'],shell=True)
if res.returncode == 0:
 print('FAILED: run had warning')
 quit()

 # read field log
if os.path.isfile('field.log') == False:
 print('FAILED: field.log not found')
 quit()
flog = np.loadtxt('field.log',comments='%')
if flog.shape[1] != 41:
 print('FAILED: unknown field.log format')
 quit()
if flog.shape[0] < 5:
 print('FAILED: too short field.log')
 quit()

# read first particle log
plogFileName = ''
for f in sorted(os.listdir('./')):
 if (f.startswith('pop001_') and f.endswith('.log')):
  plogFileName = f
if (plogFileName == '') or (os.path.isfile(plogFileName) == False):
 print('FAILED: pop001_*.log not found')
 quit()
plog = np.loadtxt(plogFileName,comments='%')
if plog.shape[1] != 17:
 print('FAILED: unknown pop001_*.log format')
 quit()
if plog.shape[0] < 5:
 print('FAILED: too short pop001_*.log')
 quit()

flog = np.loadtxt('field.log',comments='%')
plog = np.loadtxt(plogFileName,comments='%')

t = flog[:,0]
tmin = t[0]
tmax = t[-1]

#tmiddle = (t[-1] - t[0])/2.0
#ii = (np.argwhere(t >= tmiddle)).flatten()

#t = flog[ii,0]
#D = flog[ii,20]
#D = plog[:,1]

#trend = np.polyfit(t,D,1)
#trendpoly = np.poly1d(trend)
#dt = np.diff(t)
#dt = np.insert(dt,0,dt[0])
#dD = np.fabs(dt*np.gradient(D,t)/D)
#dD = np.fabs(np.gradient(D,t)/D)
#plt.subplot(2,1,1)
#plt.plot(t,D)
#plt.subplot(2,1,2)
#plt.plot(t,trendpoly(t))
#plt.show()

# loop through in N blocks
N = 10
Nb = int(np.floor(len(t)/10))
ii = np.arange(0,Nb)
Nsub = 1
for jj in range(N):
 kk = ii+Nb*jj
 t = flog[kk,0]
 D = flog[kk,9]
 trend = np.polynomial.Polynomial.fit(t,D,1)
 print(trend.coef[1])
 trendpoly = np.polynomial.Polynomial(trend.coef[::-1])
 plt.subplot(N,2,Nsub); Nsub += 1;
 plt.plot(t,D)
 plt.xlim([tmin,tmax])
 plt.subplot(N,2,Nsub); Nsub += 1;
 plt.plot(t,trendpoly(t))
 plt.xlim([tmin,tmax])
plt.show()
