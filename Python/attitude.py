import numpy as np
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import quaternion as qt

# Variables
acc = []
gyr = []
mag = []
acc_cal = []
gyr_cal = []
mag_cal = []

# Read values from data string
def get_float(data):
        val = []
        index = []
        index.append(data.index(':'))
        index.append(data.index(','))
        index.append(data[index[1]+2:].index(',')+index[1]+2)
        index.append(data[index[2]+2:].index(' ')+index[2]+2)
        val.append(float(data[index[0]+2:index[1]]))
        val.append(float(data[index[1]+2:index[2]]))
        val.append(float(data[index[2]+2:index[3]]))
        return val

# *** MAIN ***
# *** OPEN MAIN FILE ***
f_main = open('D:\Program Files (x86)\Microsoft VS Code Workspace\Python\mag.txt', 'r').readlines()
f_calibration = open('D:\Program Files (x86)\Microsoft VS Code Workspace\Python\mag_cal.txt', 'r').readlines()

# Loop through main data and extract floating values
for data in f_main:
        if 'Accel' in data:
                acc.append(get_float(data))
        elif 'Gyro' in data:
                gyr.append(get_float(data))
        elif 'Mag' in data:
                mag.append(get_float(data))
                # mag_x.append(mag[-1][1])
                # mag_y.append(mag[-1][0])
                # mag_z.append(-mag[-1][2])

# Convert lists to arrays
acc = np.asarray(acc)
gyr = np.asarray(gyr)
mag = np.asarray(mag)
mag = mag[:,[1,0,2]]

# *** OPEN CALIBRATION FILE ***
# Loop through main data and extract floating values
for data in f_calibration:
        if 'Accel' in data:
                acc_cal.append(get_float(data))
        elif 'Gyro' in data:
                gyr_cal.append(get_float(data))
        elif 'Mag' in data:
                mag_cal.append(get_float(data))
                # mag_x.append(mag[-1][1])
                # mag_y.append(mag[-1][0])
                # mag_z.append(-mag[-1][2])

# Convert lists to arrays
acc_cal = np.asarray(acc_cal)
gyr_cal = np.asarray(gyr_cal)
mag_cal = np.asarray(mag_cal)
mag_cal = mag_cal[:,[1,0,2]]

# *** SENSOR CALIBRATION ***
# Hard iron distortion
mag_x_offset = (max(mag_cal[:,0]) + min(mag_cal[:,0])) / 2
mag_y_offset = (max(mag_cal[:,1]) + min(mag_cal[:,1])) / 2
mag_z_offset = (max(mag_cal[:,2]) + min(mag_cal[:,2])) / 2

mag_x_corrected = mag_cal[:,0] - mag_x_offset
mag_y_corrected = mag_cal[:,1] - mag_y_offset
mag_z_corrected = mag_cal[:,2] - mag_z_offset

print(mag_x_offset)
print(mag_y_offset)
print(mag_z_offset)

# Raw data after calibration
plt.figure()
mag_xx = plt.scatter(mag_x_corrected,mag_y_corrected)
mag_xy = plt.scatter(mag_x_corrected,mag_z_corrected)
mag_yz = plt.scatter(mag_y_corrected,mag_z_corrected)
plt.legend([mag_xx, mag_xy, mag_yz],['mag_xy', 'mag_xz', 'mag_yz'])
plt.title('Magnetometer Calibration')
plt.show()

# Calibration compensate magnetic sensor values
mag[:,0] -= mag_x_offset
mag[:,1] -= mag_y_offset
mag[:,2] -= mag_z_offset

# Loop through data and test filter functions
att = np.zeros((len(acc),3))
for k in range(len(acc)):
        # map data to filter coordinates
        # acc_x = acc[k][0]
        # acc_y = acc[k][1]
        # acc_z = acc[k][2]

        # scale 1 g * pi
        # acc_x *= 9.81

        # Scale values
        acc[k,2] *= 9.81        # 1G=9.81m/ss

        
        # if (acc_z!=0):
        att[k,0:2] = qt.acc_to_eul(acc[k,:])
        att[k,2:3] = qt.mag_to_eul(att[k,0],att[k,1],mag[k,:],0)
        # theta.append(att[0]*180/np.pi)
        # phi.append(att[1]*180/np.pi)
        # D = np.pi/4
        # psi.append(qt.mag_to_eul(att[0],att[1],mag_x[k]-23.93,mag_y[k]+2.175,mag_z[k]+61.815,D)*180/np.pi)

        

# *** Magnetometer ***
# Raw data without calibration
plt.figure()
mag_xx = plt.scatter(mag[:,0],mag[:,1])
mag_xy = plt.scatter(mag[:,0],mag[:,2])
mag_yz = plt.scatter(mag[:,1],mag[:,2])
plt.legend([mag_xx, mag_xy, mag_yz],['mag_xy', 'mag_xz', 'mag_yz'])
plt.title('Magnetometer Raw')
plt.show()






# *** Plot ***
# Sensors raw data
tvec = np.arange(0,len(acc))
labels = ['x', 'y', 'z']
plt.figure()
plt.subplot(3, 1, 1)
lineObjects = plt.plot(tvec, acc)
plt.ylabel('Acc [g]')
plt.legend(lineObjects, labels)
plt.title('Sensors Raw Data')
plt.subplot(3, 1, 2)
lineObjects = plt.plot(tvec, gyr)
plt.ylabel('Gyr [dps]')
plt.legend(lineObjects, labels)
plt.subplot(3, 1, 3)
lineObjects = plt.plot(tvec, mag)
plt.ylabel('Mag [uT]')
plt.legend(lineObjects, labels)
plt.show()

# Filter results
tvec = np.arange(0,len(theta))
plt.figure()
plt.subplot(3,1,1)
plt.ylabel('Theta [deg]')
plt.title('Accelerometer Based Attitude')
plt.plot(tvec, theta)
plt.subplot(3,1,2)
plt.plot(tvec, phi)
plt.ylabel('Phi [deg]')
plt.subplot(3,1,3)
plt.plot(tvec, psi)
plt.ylabel('Psi [deg]')
plt.xlabel('Sample [k]')
plt.show()

