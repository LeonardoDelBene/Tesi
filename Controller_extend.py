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
        if(velocity_prec==0):
            k_prec=0
        else:
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
        return calculated_alpha

    def error(self,prec,vehicle,k):
        alpha=self.alpha(prec.controller.states[k].omega,prec.states[k].velocity,vehicle.states[k].velocity,k)
        self.states[k].error_x = prec.states[k].x +self.states[k].s_x_prec- vehicle.states[k].x - (
                (self.r + (self.h * vehicle.states[k].velocity)) * math.cos(vehicle.states[k].theta))

        self.states[k].error_y = prec.states[k].y + self.states[k].s_y_prec - vehicle.states[k].y - (
                (self.r + (self.h * vehicle.states[k].velocity)) * math.sin(vehicle.states[k].theta))

        self.states[k].error_velocity_x = prec.states[k].velocity * math.cos(prec.states[k].theta) - vehicle.states[k].velocity * math.cos(vehicle.states[k].theta + alpha)

        self.states[k].error_velocity_y = prec.states[k].velocity * math.sin(prec.states[k].theta) - vehicle.states[k].velocity * math.sin(vehicle.states[k].theta + alpha)

    def update_state_init(self,k,prec, vehicle,a,w):
        self.states.append(Controller_extend_state())
        if(vehicle.first):
            self.states[k].acceleration = a
            self.states[k].omega = w
            self.states[k].error_x = 0
            self.states[k].error_y = 0
            self.states[k].error_velocity_x = 0
            self.states[k].error_velocity_y = 0
            self.states[k].alpha = 0
            self.states[k].s_x_prec = 0
            self.states[k].s_y_prec = 0
        else:
            self.s(prec.states[k].theta,prec.states[k].velocity,prec.controller.states[k].omega,vehicle.states[k].velocity,k)
            self.error(prec,vehicle,k)
            self.states[k].alpha = self.alpha(prec.controller.states[k].omega,prec.states[k].velocity,vehicle.states[k].velocity,k)
            self.get_acceleration_omega(prec,vehicle,k)

    def get_acceleration_omega(self, prec, vehicle, k):
        if(k!=0):
            kk=k+1
        else:
            kk=0
        theta = vehicle.states[k].theta
        theta_prec = prec.states[k].theta
        velocity = vehicle.states[k].velocity
        alpha = self.states[k].alpha

        u= self.h * (self.r + (self.h*velocity))*(1-math.sin(alpha)*math.sin(theta_prec - theta))
        T_inv = np.array([
            [((self.r + (self.h * velocity))*math.cos(theta))/u , ((self.r + (self.h * velocity))*math.sin(theta))/u],
            [((-self.h * math.sin(theta))-((self.h * math.sin(alpha))*math.cos(theta_prec)))/u , ((self.h * math.cos(theta))-((self.h * math.sin(alpha))*math.sin(theta_prec)))/u]
        ])

        K= np.array([
            [self.k1 * self.states[k].error_x],
            [self.k2 * self.states[k].error_y]
        ])

        Z= np.array([
            [self.states[k].error_velocity_x/math.cos(alpha)],
            [self.states[k].error_velocity_y/math.cos(alpha)]
        ])

        O= np.array([
            [-math.sin(theta)],
            [math.cos(theta)]
        ])

        R= np.array([
            [math.cos(theta_prec), -math.sin(theta_prec)],
            [math.sin(theta_prec), math.cos(theta_prec)]
        ])

        if(prec.controller.states[k].omega==0):
            acceleration_angular = 0
        else:
            radius = prec.states[k].velocity / prec.controller.states[k].omega
            acceleration_angular = prec.states[k].velocity / radius

        k_i= ((acceleration_angular * prec.states[k].velocity) - (prec.controller.states[k].omega * prec.controller.states[k].acceleration))/(prec.states[k].velocity ** 2)

        if(prec.controller.states[k].omega==0):
            s_k=0
        else:
            s_k = (((1 - math.cos(alpha)) / (
                        prec.controller.states[k].omega / prec.states[k].velocity) ** 2))

        S= np.array([
            [self.s_magnitude(velocity,prec.states[k].velocity,prec.controller.states[k].omega) * prec.controller.states[k].omega],
            [-s_k * k_i]
        ])

        O_prec = np.array([
            [math.cos(theta_prec)],
            [math.sin(theta_prec)]
        ])

        B= (O * (velocity * math.tan(alpha))) + (np.dot(R,S)) + (((1-(1/math.cos(alpha)))*O_prec) * prec.states[k].velocity)

        M = K + Z + B
        A = np.dot(T_inv,M)

        self.states[kk].acceleration = A[0,0]
        self.states[kk].omega = A[1,0]
        return A

    def updateState(self, k, T, prec, vehicle, a, w):
        self.states.append(Controller_extend_state())
        if (vehicle.first):
            self.states[k].error_x = 0
            self.states[k].error_y = 0
            self.states[k].error_velocity_x = 0
            self.states[k].error_velocity_y = 0
            self.states[k].acceleration = a
            self.states[k].omega = w
            self.states[k].alpha = 0
            self.states[k].s_x_prec = 0
            self.states[k].s_y_prec = 0

        else:
            self.states[k].error_x = self.states[k - 1].error_x + (T * (-self.k1 * self.states[k - 1].error_x))
            self.states[k].error_y = self.states[k - 1].error_y + (T * (-self.k2 * self.states[k - 1].error_y))

            self.states[k].alpha= self.alpha(prec.controller.states[k-1].omega, prec.states[k-1].velocity, vehicle.states[k-1].velocity, k)
            alpha=self.states[k-1].alpha


            H = np.array([
                [math.cos(prec.states[k-1].theta),
                 -(prec.states[k-1].velocity * math.sin(prec.states[k-1].theta))],
                [math.sin(prec.states[k-1].theta),
                 (prec.states[k-1].velocity * math.cos(prec.states[k-1].theta))]
            ])

            A = np.array([
                [prec.controller.states[k-1].acceleration],
                [prec.controller.states[k-1].omega]
            ])

            O = np.array([
                [math.sin(vehicle.states[k - 1].theta + alpha)],
                [-math.cos(vehicle.states[k - 1].theta + alpha)]
            ])

            if (prec.controller.states[k-1].omega == 0):
                acceleration_angular = 0
            else:
                radius = prec.states[k-1].velocity / prec.controller.states[k-1].omega
                acceleration_angular = prec.states[k-1].velocity / radius

            if (prec.states[k-1].velocity == 0):
                k_i = 0
            else:
                k_i = ((acceleration_angular * prec.states[k-1].velocity) - (
                        prec.controller.states[k-1].omega * prec.controller.states[k-1].acceleration)) / (
                              prec.states[k-1].velocity ** 2)

            b=vehicle.states[k - 1].velocity * (self.r + (self.h * vehicle.states[k - 1].velocity)) * (
                    math.cos(alpha) ** 2)

            B2 = b * O * k_i

            u = self.h * (self.r + (self.h * vehicle.states[k-1].velocity)) * (
                        1 - math.sin(alpha) * math.sin(prec.states[k-1].theta - vehicle.states[k-1].theta))
            T12_inv = np.array([
                [((self.r + (self.h * vehicle.states[k-1].velocity)) * math.cos(vehicle.states[k-1].theta)) / u,
                 ((self.r + (self.h * vehicle.states[k-1].velocity)) * math.sin(vehicle.states[k-1].theta)) / u],
                [((-self.h * math.sin(vehicle.states[k-1].theta)) - ((self.h * math.sin(alpha)) * math.cos(prec.states[k-1].theta))) / u,
                 ((self.h * math.cos(vehicle.states[k-1].theta)) - ((self.h * math.sin(alpha)) * math.sin(prec.states[k-1].theta))) / u]
            ])

            d = vehicle.states[k - 1].velocity * self.h * (
                        prec.controller.states[k-1].omega / prec.states[k-1].velocity) * (
                        math.cos(alpha) ** 2)

            T34 = np.array([
                [math.cos(vehicle.states[k - 1].theta + alpha) - d * math.sin(
                    vehicle.states[k - 1].theta + alpha),
                 -vehicle.states[k - 1].velocity * math.sin(vehicle.states[k - 1].theta + alpha)],
                [math.sin(vehicle.states[k - 1].theta + alpha) + d * math.cos(
                    vehicle.states[k - 1].theta + alpha),
                 vehicle.states[k - 1].velocity * math.cos(vehicle.states[k - 1].theta + alpha)]
            ])

            O_1 = np.array([
                [-math.sin(vehicle.states[k - 1].theta)],
                [math.cos(vehicle.states[k - 1].theta)]
            ])

            R = np.array([
                [math.cos(prec.states[k-1].theta), -math.sin(prec.states[k-1].theta)],
                [math.sin(prec.states[k-1].theta), math.cos(prec.states[k-1].theta)]
            ])

            if (prec.controller.states[k-1].omega == 0):
                s_k = 0
            else:
                s_k = (-((1 - math.cos(alpha)) / (
                        prec.controller.states[k-1].omega / prec.states[k-1].velocity) ** 2))

            if(prec.states[k-1].velocity==0):
                k_i=0
            else:
                k_i = ((acceleration_angular * prec.states[k-1].velocity) - (
                        prec.controller.states[k-1].omega * prec.controller.states[k-1].acceleration)) / (
                              prec.states[k-1].velocity ** 2)

            S = np.array([
                [self.s_magnitude(vehicle.states[k - 1].velocity, prec.states[k-1].velocity,
                                  prec.controller.states[k-1].omega) *
                 prec.controller.states[k-1].omega],
                [s_k * k_i]
            ])
            O_prec = np.array([
                [math.cos(prec.states[k-1].theta)],
                [math.sin(prec.states[k-1].theta)]
            ])

            B_res1 =O_1 * (vehicle.states[k - 1].velocity * math.tan(alpha))
            B_res2 = np.dot(R, S)
            B_res3 = (1 - (1 / math.cos(alpha))) * O_prec * prec.states[k-1].velocity

            B1 = B_res1 + B_res2 + B_res3

            K = np.array([
                [self.k1 * self.states[k - 1].error_x],
                [self.k2 * self.states[k - 1].error_y]
            ])

            Z = np.array([
                [self.states[k - 1].error_velocity_x],
                [self.states[k - 1].error_velocity_y]
            ])

            L = np.dot(T34, T12_inv)

            Res0=np.dot(H,A)
            Res1=B2
            Res2=np.dot(L,K)
            Res3=(1 / math.cos(alpha) * Z)+B1
            Res4=np.dot(L,Res3)

            Result = Res0 + Res1 -Res2 - Res4

            self.states[k].error_velocity_x = self.states[k - 1].error_velocity_x + (T * Result[0, 0])
            self.states[k].error_velocity_y = self.states[k - 1].error_velocity_y + (T * Result[1, 0])









