import numpy as np
import math
import time
# import numba as nb
# import cupy as cp
import matplotlib.pyplot as plt


start= time.time()
# ##BTCS
# @nb.njit(nopython =True, parallel=True)
# def BTCS(A,B,U,M,N):
#     for i in range(M-1):
#         B[0]=-U[i,1]-U[i+1,0]
#         B[N-3]=-U[i,N-2]-U[i+1,N-1]
#         for j in range (1,N-3):
#             B[j]= -U[i, j+1]
#         W=np.zeros((N-2,1))
#         W=np.linalg.solve(A,B)
#         for k in range(1,N-1):
#             U[i+1, k]=W[k-1]
#     return U


L=1
T=1
alpha =1
dx= 0.01
dt=0.00001
d= (alpha*dt)/(dx**2)
M= math.ceil(abs(T)/dt)+1 ## number of time nodes
N=math.ceil(L/dx)+1  # no of sapce nodes
print(M,N)
x=np.zeros((N,1))
t=np.zeros((M,1))
for i in range(N):
    x[i]=0+i*dx
for j in range(M):
    t[j]=0+j*dt
U=np.zeros((M,N))
U[:,0]=0
U[:, N-1]=0
for i in range(1,N-1):
    U[0,i]=math.sin(4*math.pi*x[i])

## FTCS
for n in range(M-1):
    for xi in range(1,N-1):
        U[n+1, xi]= d*U[n,xi+1]+(1-2*d)*U[n,xi]+ d*U[n,xi-1]



# A=np.zeros((N-2, N-2))
# B=np.zeros((N-2,1))
# A[0,0]=(1+2*d)
# A[0,1]=-d
# A[N-3, N-4]=-d
# A[N-3, N-3]= (1+2*d)
# for i in range(1,N-3):
#     A[i,i-1]=-d
#     A[i,i]=(1+2*d)
#     A[i,i+1]=-d

# # # U=BTCS(A,B,U,M,N)
    

    
# # BTCS  
# for i in range(M-1):
#     B[0]=U[i,1]+U[i+1,0]
#     B[N-3]=U[i,N-2]+U[i+1,N-1]
#     for j in range (1,N-3):
#         B[j]= U[i, j+1]
#     W=np.linalg.solve(A,B)
#     for k in range(1,N-1):
#         U[i+1, k]=W[k-1]

# Crank Nicolson
# d=d/2
# A=np.zeros((N-2, N-2))
# B= np.zeros((N-2,1))
# A[0,0]=(1+2*d)
# A[0,1]=-d
# A[N-3, N-4]=-d
# A[N-3, N-3]=1+2*d
# for i in range(1,N-3):
#     A[i,i-1]=-d
#     A[i,i]=1+2*d
#     A[i,i+1]=-d

# for i in range(M-1):
#     B[0]=(1-2*d)*U[i,1]+d*U[i,2]+d*(U[i+1,0]+U[i,0])
#     for j in range(N-4):
#         B[j+1]=d*U[i,j+1]+(1-2*d)*U[i,j+2]+d*U[i,j+3]
#     B[N-3]= d*U[i,N-3]+(1-2*d)*U[i,N-2]+d*(U[i+1,N-1] +U[i,N-1])
#     W=np.linalg.solve(A,B)
#     for k in range(1,N-1):
#         U[i+1, k]=W[k-1]
    
print(U)

end=time.time()
print(end-start)

plt.plot(x, U[M-1,:])
plt.title("X vs U at the end time")
plt.xlabel("x")
plt.ylabel("U at end time")

plt.show()

plt.plot(t, U[:, N//2])
plt.title("T vs U at the mid point")
plt.xlabel("t")
plt.ylabel("U at mid point")
plt.show()


