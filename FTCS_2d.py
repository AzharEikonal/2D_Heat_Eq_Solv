import numpy as np
import math
import time
import matplotlib.pyplot as plt

L_x=1
L_y=1
T=1
alpha=1
M=50 
N=50
dx=L_x/(M-1)
dy=L_y/(N-1)

dt=0.0001

L=math.ceil(abs(T/dt))+1
x=np.zeros((M,1))
y=np.zeros((N,1))
t=np.zeros((L,1))
for i in range(M):
    x[i]=0+i*dx
for j in range(N):
    y[j]=0+j*dy
for k in range(L):
    t[k]=0+k*dt
U=np.zeros((M,N))

for i in range(1,M-1):
    for j in range(1,N-1):
        U[i,j]= np.pi*i*dx+ 2*np.pi*j*dy
U[:,0]=0
U[:,-1]=0
U[0,:]=0.1
U[-1,:]=0.1

Dx=alpha*dt/(dx**2)
Dy=alpha*dt/(dy**2)
# FTCS
start_time=time.time()
for n in range(L-1):
    for i in range(1,M-1):
        for j in range(1,N-1):
            U[i,j]=Dx*(U[i+1,j]-2*U[i,j]+U[i-1,j])+Dy*(U[i,j+1]-2*U[i,j]+U[i,j-1])+U[i,j]
print("FTCS time:",time.time()-start_time)
print(U)



# Plot the solution
x = np.linspace(0, L_x, M)
y = np.linspace(0, L_y, N)
X, Y = np.meshgrid(x, y)

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, U, cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Temperature')
ax.set_title('2D Heat Equation Solution using FTCS Scheme')
plt.show()




