import numpy as np
import matplotlib.pyplot as plt

def next_pow_2(x):
    return np.power(2.,np.ceil(np.log(x)/np.log(2.)));

d=np.power(2,np.arange(13))
N=2048
M=8
c=(d*next_pow_2(M+N/d-1)+next_pow_2(M+N/d-1)*np.log(next_pow_2(M+N/d-1)))/M
plt.plot(np.log(d)/np.log(2),c,label='cost')
plt.plot(np.log(d)/np.log(2),N/d,label='ir part length')
plt.plot(np.log(d)/np.log(2),np.ones((len(d),))*M,label='hopsize')
plt.legend()
plt.show()


