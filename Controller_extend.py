import math
import numpy as np
from Controller_extend_state import Controller_extend_state


class Controller_extend:
    def __init__(self, r, h, k1, k2):
        self.states = []
        self.r = r
        self.h = h
        self.k1 = k1
        self.k2 = k2



    def s_magnitude(self,velocity,velocity_prec, omega_prec):
        k_prec= omega_prec/velocity_prec
        if(k_prec==0):
            return 0
        else:
            return ((-1 + math.sqrt(1+((k_prec**2)*(self.r + self.h*velocity)**2)))/(k_prec))
    def s(self,theta_prec,velocity_prec,omega_prec,velocity,k):
        s_magnitude_prec = self.s_magnitude(velocity,velocity_prec,omega_prec)
        self.states[k].s_x_prec = s_magnitude_prec * math.sin(theta_prec)
        self.states[k].s_y_prec = -s_magnitude_prec * math.cos(theta_prec)

    def alpha(self,omega_prec,velocity_prec,velocity,k):
        k_prec= omega_prec/velocity_prec
        calculated_alpha = math.atan((k_prec*(self.r + self.h*velocity)))

        if calculated_alpha < -math.pi / 2:
            self.states[k].alpha = -math.pi / 2
        elif calculated_alpha > math.pi / 2:
            self.states[k].alpha = math.pi / 2
        else:
            self.states[k].alpha = calculated_alpha

    def error(self,prec,vehicle,k):
        self.states[k].error_x = prec.states[k].x +self.s_x_prec - vehicle.states[k].x - (
                (self.r + (self.h * vehicle.states[k].velocity)) * math.cos(vehicle.states[k].theta))

        self.states[k].error_y = prec.states[k].y + self.s_y_prec - vehicle.states[k].y - (
                (self.r + (self.h * vehicle.states[k].velocity)) * math.sin(vehicle.states[k].theta))

        self.states[k].error_velocity_x = prec.states[k].velocity * math.cos(prec.states[k].theta) - vehicle.states[k].velocity * math.cos(vehicle.states[k].theta + self.alpha(prec.states[k].omega,prec.states[k].velocity,vehicle.states[k].velocity,k))

        self.states[k].error_velocity_y = prec.states[k].velocity * math.sin(prec.states[k].theta) - vehicle.states[k].velocity * math.sin(vehicle.states[k].theta + self.alpha(prec.states[k].omega,prec.states[k].velocity,vehicle.states[k].velocity,k))


    def update_state_initi(self,k,prec, vehicle,a,w):
        self.states.append(Controller_extend_state())
        if(vehicle.first):
            self.states[k].acceleration = a
            self.states[k].omega = w
            self.states[k].error_x = 0
            self.states[k].error_y = 0
            self.states[k].error_velocity_x = 0
            self.states[k].error_velocity_y = 0
            self.states[k].alpha = 0
        else:
            self.s(prec.states[k].theta,prec.states[k].velocity,prec.states[k].omega,vehicle.states[k].velocity,k)
            self.alpha(prec.states[k].omega,prec.states[k].velocity,vehicle.states[k].velocity,k)
            self.error(prec,vehicle,k)
            self.get_acceleration_omega(prec,vehicle,k)

    def get_acceleration_omega(self, prec, vehicle, k):
        theta = vehicle.states[k].theta
        theta_prec = prec.states[k].theta
        velocity = vehicle.states[k].velocity


        T= np.array([
            [self.h * math.cos(theta) - self.h*math.sin(self.states[k].alpha)*math.sin(theta_prec), -(self.r + self.h * velocity) * math.sin(theta)],
            [self.h * math.sin(theta) + self.h*math.sin(self.states[k].alpha)*math.cos(theta_prec), (self.r + self.h * velocity) * math.cos(theta)]
        ])

        T_inv = np.linalg.inv(T)

        K= np.array([
            [self.k1 * self.states[k].error_x],
            [self.k2 * self.states[k].error_y]
        ])

        Z= np.array([
            [self.states[k].error_velocity_x/math.cos(self.states[k].alpha)],
            [self.states[k].error_velocity_y/math.cos(self.states[k].alpha)]
        ])

        O= np.array([
            [-math.sin(theta)],
            [math.cos(theta)]
        ])

        R= np.array([
            [math.cos(theta_prec), -math.sin(theta_prec)],
            [math.sin(theta_prec), math.cos(theta_prec)]
        ])

        radius = prec.states[k].velocity / prec.states[k].omega
        acceleration_angular = prec.states[k].velocity / radius
        S= np.array([
            [self.s_magnitude(velocity,prec.states[k].velocity,prec.states[k].omega) * prec.states[k].omega],
            [(-((1-math.cos(self.states[k].alpha))/(prec.states[k].omega/prec.states[k].velocity)**2))*((acceleration_angular)/(prec.states[k].acceleration))]
        ])

        O_prec = np.array([
            [math.cos(theta_prec)],
            [math.sin(theta_prec)]
        ])

        B= O * (velocity * math.tan(self.states[k].alpha)) + np.dot(R,S) + (1-(1/math.cos(self.states[k].alpha)))*O_prec * prec.states[k].velocity

        A = np.dot(T_inv, (K + Z + B))

        self.states[k].acceleration = A[0,0]
        self.states[k].omega = A[1,0]
        return A