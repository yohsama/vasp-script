import numpy as np
import statsmodels.tsa.api as smt
import ultils
import ase.io
from numba import jit
def velocity_from_newton(dx, a, dt):
    return (dx - 0.5 * a * (dt**2)) / dt


def velocity_from_infinite(dx, dt):
    return dx / dt


def get_velocitys(strcs, dt=1,method='newton'):
    x = np.array([i.get_scaled_positions() for i in strcs])  # frac
    dx = (ultils.get_nearest_distence(x[1:] - x[:-1]) @ strcs[0].cell) * (
        10**-10)  # \AA   -> m
    dt = 10**-15*dt  # fs    -> s
    T = len(strcs)
    if method == "newton":
        f = np.array([i.get_forces() for i in strcs
                  ]) * 1.6021766208 * (10**-19) / (10**-10)  # eV / \AA -> J/m
        m = strcs[0].get_masses() * 1.66053904 * (10**-27)  # a.u   -> kg
        a = f / m[None, :, None]  # J/m/kg = m*s^-2
        v = velocity_from_newton(dx, a[:-1], dt)
    elif method == "back":
        v = velocity_from_infinite(dx, dt)
    return v

@jit(nopython=True) 
def uacf(tdE,tot):
    dE=tdE[len(tdE)-tot:]
    dEbar=dE.mean()
    Cut=np.zeros(tot)
    for i in range(tot):
        temp=0
        temp=(dE[i+1:]-dEbar)*(dE[:-i-1]-dEbar)
        Cut[i]=temp.sum()
    return Cut/tot

#@jit
def get_acf_all(v):
    vsize=v.shape[0]
    v=v.reshape(vsize,-1)
    allacf=np.array([np.correlate(i,i,mode='full')[vsize-1:] for i in v.T],)
    return allacf

#@jit(nopython=True)
#def get_acf_all(v):
#    vsize=v.shape[0]
#    v=v.reshape(vsize,-1)
#    allacf=[uacf(i,vsize) for i in v.T]
#    return allacf


def v2f(vacf,dt=1):
    vacf /= vacf[0]
    T = vacf.shape[0]
    x = np.arange(0, 1000 / dt, 1000 / T * dt) / 3 * 100
    fft = np.abs(np.fft.fft(vacf))
    return x, fft

#if __name__ == "__main__":
#    OUTCAR='OUTCAR'
#    strcs=ase.io.read('OUTCAR',':')
#    v=get_velocitys(strcs)
#    np.save('v.npy',v)
#    allacf=get_acf_all(v).sum(0)
#    x,fft=v2f(allacf)
#    np.savetxt('allacf.dat',allacf)
#    np.savetxt('fftacf.dat',np.array([x,fft]).T)
#    plt.plot(x,fft)
#    plt.show()
