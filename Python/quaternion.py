import numpy as np


# **** ATTITUDE TRANSFORMATION MATRIX ****
# Y-axis
def Cy(psi):
    return np.array([[np.cos(psi),0,np.sin(psi)], [0,1,0], [-np.sin(psi),0,np.cos(psi)]])

# Z-axis
def Cz(theta):
    return np.array([[np.cos(theta),-np.sin(theta),0], [np.sin(theta),np.cos(theta),0], [0,0,1]])

# X-axis
def Cx(phi):
    return np.array([[1,0,0],[0,np.cos(phi),-np.sin(phi)],[0,np.sin(phi),np.cos(phi)]])

# Complete represenation body frame to navigation frame
def Cbn(phi,theta,psi):
    return np.dot(np.dot(Cx(psi),Cz(theta)),Cy(psi))



# **** ATTITUDE TRANSFORMATION MATRIX QUATERNION FORM ****
# Complete represenation body frame to navigation frame
def Cbn_q(q):
    a = q[0]**2+q[1]**2-q[2]**2-q[3]**2
    b = 2*q[0]*q[3]+2*q[1]*q[2]
    c = -2*q[0]*q[2]+2*q[1]*q[3]
    d = -2*q[0]*q[3]+2*q[1]*q[2]
    e = q[0]**2-q[1]**2+q[2]**2-q[3]**2
    f = 2*q[0]*q[1]+2*q[2]*q[3]
    g = 2*q[0]*q[2]+2*q[1]*q[3]
    h = -2*q[0]*q[1]+2*q[2]*q[3]
    i = q[0]**2-q[1]**2-q[2]**2+q[3]**2
    return np.array([[a,b,c], [d,e,f], [g,h,i]])



# **** QUATERNION TO EULERS ****
# returns = [psi,theta,phi]
def quat_to_eul(q):
    return np.array([[np.arctan((2*q[0]*q[2]-2*q[1]*q[3])/(q[0]**2+q[1]**2-q[2]**2-q[3]**2))],
    [np.arcsin(2*q[0]*q[3]+2*q[1]*q[2])],
    [np.arctan((2*q[0]*q[1]-2*q[2]*q[3])/(q[0]**2-q[1]**2+q[2]**2-q[3]**2))]])



# **** ACCELEROMETER BASED ATTITUDE  ****
# returns = [theta,phi]
def acc_to_eul(acc):
    g = 9.81
    return np.array([[np.arcsin(acc[0]/g)],[np.arctan(-acc[1]/acc[2])]])

# **** MAGNETOMETER / ACCELEROMETER BASED HEADING   ****
def mag_to_eul(theta,phi,mag,D):
    X = np.array([[np.cos(theta),np.sin(theta),0],[-np.cos(phi)*np.sin(theta),np.cos(phi)*np.cos(theta),np.sin(phi)],[np.sin(theta)*np.sin(phi),-np.sin(phi)*np.cos(theta),np.cos(phi)]])
    Hb = np.array([[mag[0]],[mag[1]],[mag[2]]])
    Hl = np.dot(X,Hb)
    return np.arctan(Hl[0]/Hl[1])+D

# *** Process Model ***
# returns q[k+1]
def process_model(gyr, q, w, t):
    omega_skew = np.array([[0, gyr[2], -gyr[1]],[-gyr[2], 0, gyr[0]],[gyr[1], -gyr[0], 0]])
    omega = 0.5 * np.array([[0,-np.transpose(gyr)],[gyr, omega_skew]])
    return (np.identity(4)+0.5*omega*t)*q+w 