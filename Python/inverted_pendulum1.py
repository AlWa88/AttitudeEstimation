# *********************************
#   INVERTED PENDULUM SIMULATION
# *********************************

# Papers:
# "Optimal Control of Nonlinear Inverted Pendulum System Using 
# PID Controller and LQR: Performance Analysis Without and With Disturbance Input"
# https://link.springer.com/content/pdf/10.1007%2Fs11633-014-0818-1.pdf
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# import FuncAnimation, FFMpegFileWriter

# *********************************
#   NONLINEAR MODEL
# *********************************
#
# State vector: [x1,x2,x3,x4]^T = [x,xdot,theta,thetadot]
# Input: u = F
# Constants:
# M = cart mass [kg]
# m = pendulum mass [kg]
# l = pendulum lengt [m]
# I = pendulum inertia

# Constants
M = 2.0
m = 0.25
l = 0.5
I = 0
g = 9.81
# Simulation length and time step
T = 10000
ts = 0.01
cWi = 0.1
cHi = 0.1
xmin = -10.0
xmax = 10.0

# Initial conditions
xFE = np.zeros((4,T))
xRK2 = np.zeros((4,T))
xRK4 = np.zeros((4,T))
u = np.zeros((2,T))
xPid = np.zeros((2,T))
sp = np.ones((1,T))*0.5

#xRK2[:,0] = [0,0,np.pi*0.95,0]
xRK2[:,0] = [0,0.1,0,0]

# Generate some white noise
np.random.seed(10)
Fw = np.random.normal(0, 0.01, size=T)

# Plots / Animation
plotSolvers = False
plotNoise = False
plotPID = True
plotAnimate = True

# Nonlinear state equation function
def nonlinearModel(x,u,Fw):
    return np.array([
    x[1],
    (u[0] + m*l*np.sin(x[2])*(x[3]**2) - m*g*np.cos(x[2])*np.sin(x[3])) / (M + m - m*(np.cos(x[2])**2)),
    x[3],
    (u[0]*np.cos(x[2]) - (M + m)*g*np.sin(x[2]) + m*l*(np.cos(x[2])*np.sin(x[2]))*(x[3]**2)) / (m*l*(np.cos(x[2])**2) - (M + m)*l)
    ])

# Nonlinear state equation function with horizontal wind disturbance on 
def nonlinearModelDisturbance(x,u,Fw):
    return np.array([
    x[1],
    (u + m*l*np.sin(x[2])*(x[3]**2) - m*g*np.cos(x[2])*np.sin(x[2]) + Fw*(np.sin(x[2])**2)) / (M + m - m*(np.cos(x[2])**2)),
    x[3],
    (u*np.cos(x[2]) - (M + m)*g*np.sin(x[2]) + m*l*np.cos(x[2])*np.sin(x[2])*(x[3]**2) - (M/m)*Fw*np.cos(x[2])) / (m*l*(np.cos(x[2])**2) - (M + m)*l)
    ])

# Nonlinear state equation function with "fixed" cart
def nonlinearModelFixedCart(x,u,Fw):
    return np.array([
    x[1],
    u[0] / (M + m),
    x[3],
    (u[0]*np.cos(x[2]) - (M + m)*g*np.sin(x[2]) + m*l*(np.cos(x[2])*np.sin(x[2]))*(x[3]**2)) / (m*l*(np.cos(x[2])**2) - (M + m)*l)
    ])

# General saturatin function
def saturation(x, min, max):
    if x < min:
        return min
    elif x > max:
        return max
    else:
        return x

# General limit function
def limit(x, min, max):
    if x < min:
        return x + max
    elif x > max:
        return x - max
    else:
        return x

# General scale function
# xMin/xMax is the measurement resolution
# pMin/pMax is the scale resolution
def scale(x, xMin, xMax, pMin, pMax):
    return (pMin - pMax) / (xMin - xMax) * x + pMin


# PID controller
class PID():
    # Constructor
    def __init__(self, Kp, Ki, Kd, uMin, uMax):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.integral = 0
        self.errorPrev = 0
        self.u = 0
        self.uMin = uMin
        self.uMax = uMax

    # PID output
    def output(self, x, sp, ts):
        error = x - sp
        # Integral windup protection
        if self.u > self.uMin or self.u < self.uMax:
            self.integral = self.integral + error*ts
        derivative = (error - self.errorPrev)/ts
        self.errorPrev = error
        self.u = self.Kp*error + self.Ki*self.integral + self.Kd*derivative
        # Controller saturation
        self.u = saturation(self.u, self.uMin, self.uMax)
        return self.u

# Create instances of the controllers
pidAngle = PID(100.0,10.0,50.0, -2.0, 2.0)
pidPosition = PID(10.0,2,5.0, -0.05, 0.05)







# *** SIMULATION ****
# Three different solvers to see which is more accurate
for k in range(T-1):
    # Calculate control ouput
    u[0,k] = -pidPosition.output(xRK2[0,k], sp[0,k], ts)
    #positionsOutput = 0

    if xRK2[2,k] > np.pi:
        xPid[0,k] = xRK2[2,k] - 2*np.pi
    else:
        xPid[0,k] = xRK2[2,k]

    u[0,k] = 0.0

    u[1,k] = pidAngle.output(xPid[0,k], u[0,k], ts) 

    # Forward Eulers
    xFE[:,k+1] = xFE[:,k] + ts*nonlinearModelDisturbance(xFE[:,k],u[1,k],Fw[k])

    # 2nd Order Runga-Kutta (RK2) (Heun's method)
    k1 = nonlinearModelDisturbance(xRK2[:,k],u[1,k],Fw[k]) 
    k2 = nonlinearModelDisturbance(xRK2[:,k]+ts*k1,u[1,k],Fw[k]) 
    xRK2[:,k+1] = xRK2[:,k] + ts*(k1/2 + k2/2)

    # 4th Order Runga-Kutta (RK4) (classic)
    k1 = nonlinearModelDisturbance(xRK4[:,k],u[1,k],Fw[k])  
    k2 = nonlinearModelDisturbance(xRK4[:,k]+ts*k1/2,u[1,k],Fw[k])
    k3 = nonlinearModelDisturbance(xRK4[:,k]+ts*k2/2,u[1,k],Fw[k])
    k4 = nonlinearModelDisturbance(xRK4[:,k]+ts*k3,u[1,k],Fw[k])
    xRK4[:,k+1] = xRK4[:,k] + ts*(k1/6+k2/3+k3/3+k4/6)

    # Saturate position within xmin and xmax and theta within 0,2pi
    xFE[0,k+1] = saturation(xFE[0,k+1], xmin+cWi, xmax-cWi)
    xRK2[0,k+1] = saturation(xRK2[0,k+1], xmin+cWi, xmax-cWi)
    xRK4[0,k+1] = saturation(xRK4[0,k+1], xmin+cWi, xmax-cWi)
    xFE[2,k+1] = limit(xFE[2,k+1], 0, 2*np.pi)
    xRK2[2,k+1] = limit(xRK2[2,k+1], 0, 2*np.pi)
    xRK4[2,k+1] = limit(xRK4[2,k+1], 0, 2*np.pi)






# *** ANIMATION ****
if plotAnimate:
    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    #ax2.grid()

    l1,=ax1.plot([],[],'o')
    l2,=ax1.plot([],[],'o')

    #ax1 = plt.axes(xlim=(-1,1), ylim=(-1,1))
    #ax1.set_
    #ax1 = fig.
    #ax2 = plt.plot([], [], lw=2)



    
    #line, = ax1.plot([], [], 'o', lw=2)
    #plt.xlabel('x [m]')
    #plt.ylabel('y [m]')
    #plt.title('Inverted Pendulum')
    

    # # Prepare two plots (pendulum and cart)
    # plotlays, plotcols = [3], ["black","black","black"]
    # lines = []
    # for index in range(1):
    #     if(index<2):
    #         lobj = ax1.plot([],[],lw=1,color=plotcols[index])[0]
    #     else:
    #         lobj = ax1.plot([],[],'o', lw=2,color=plotcols[index])[0]
    #     lines.append(lobj)

    # def init():
    #     for line in lines:
    #         line.set_data([],[])
    #     return lines

    def animate(frame,l1,l2):
        # Pendulum coordinates
        xp = [xRK2[0,frame], l*np.sin(xRK2[2,frame]) + xRK2[0,frame]]
        yp = [0, l*np.cos(xRK2[2,frame]) + 0]
        
        # Car coordinates
        xc = [xRK2[0,frame]-cWi,xRK2[0,frame]+cWi,xRK2[0,frame]+cWi,xRK2[0,frame]-cWi,xRK2[0,frame]-cWi]
        yc = [0, 0, -cHi, -cHi, 0]

        xpm = [l*np.sin(xRK2[2,frame]) + xRK2[0,frame]]
        ypm = [l*np.cos(xRK2[2,frame]) + 0]

        xlist = [xp,xc, xpm]
        ylist = [yp, yc, ypm]

        # l1.set_data([0, 1], [0.5, 0.5])
        l1.set_data(xlist[0], ylist[0])
        l2.set_data(xlist[1], ylist[1])
        #line.set_data(xp[0],yp[0])
        #ax1.set_xlim(xp[0]-1, xp[0]+1)

        ax1limmin, ax1limmax = 0.0,3.0
        ax2.set_xlim(ax1limmin,ax1limmax)

        #for lnum,line in enumerate(lines):
        #    line.set_data(xlist[lnum], ylist[lnum]) # set data for each line separately.
        #    ax1.set_xlim(-2.0,2.0)
            #line.axes.set_xlim(-1,1)

        return fig,


    ani = animation.FuncAnimation(fig, animate, fargs=(l1,l2),
                        frames=T, interval=10, blit=True)

    plt.show()



# *** Plot ****
tvec = np.arange(0,T*ts,ts)

# Solvers
if plotSolvers:
    plt.figure()
    plt.subplot(4, 1, 1)
    #plt.plot(tvec, xFE[0,:], label='FE')
    plt.plot(tvec, xRK2[0,:], label='RK2')
    #plt.plot(tvec, xRK4[0,:], label='RK4')
    plt.legend()
    plt.ylabel('x [m]')
    plt.title('Inverted Pendulum')

    plt.subplot(4, 1, 2)
    #plt.plot(tvec, xFE[1,:], label='FE')
    plt.plot(tvec, xRK2[1,:], label='RK2')
    #plt.plot(tvec, xRK4[1,:], label='RK4')
    plt.legend()
    plt.ylabel('xdot [m/s]')

    plt.subplot(4, 1, 3)
    #plt.plot(tvec, xFE[2,:], label='FE')
    plt.plot(tvec, xRK2[2,:], label='RK2')
    #plt.plot(tvec, xRK4[2,:], label='RK4')
    plt.legend()
    plt.ylabel('theta [rad]')

    plt.subplot(4, 1, 4)
    #plt.plot(tvec, xFE[3,:], label='FE')
    plt.plot(tvec, xRK2[3,:], label='RK2')
    #plt.plot(tvec, xRK4[3,:], label='RK4')
    plt.legend()
    plt.ylabel('thetadot [rad/s]')
    plt.xlabel('Steps [k]')
    plt.show()

# White noise
if plotNoise:
    plt.figure()
    plt.plot(tvec,Fw)
    plt.ylabel('Noise [N]')
    plt.xlabel('Steps [k]')
    plt.title('Horizontal White Noise on Pendulum')
    plt.show()


# Control meas, input and setpoint
if plotPID:
    plt.figure()
    plt.subplot(3, 1, 1)
    plt.plot(tvec, sp[0,:], '--', label='sp [m]')
    plt.plot(tvec, xRK2[0,:], label='x [m]')
    plt.legend()
    plt.ylabel('Position')
    plt.title('Inverted Pendulum - Control meas, input, setpoint')

    plt.subplot(3, 1, 2)
    plt.plot(tvec, u[0,:], '--', label='sp [rad]')
    plt.plot(tvec, xPid[0,:], label='theta scaled [rad]')
    plt.legend()
    plt.ylabel('Angle')

    plt.subplot(3, 1, 3)
    plt.plot(tvec, u[1,:], label='u [N]')
    plt.legend()
    plt.ylabel('Control Input')
    plt.xlabel('Time [s]')
    plt.show()

