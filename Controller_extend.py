import math
import numpy as np

from Controller import Controller
from Controller_extend_state import Controller_extend_state


class Controller_extend(Controller):
    def __init__(self, r, h, k1, k2):
        super().__init__(r, h, k1, k2)

    def s(self, k, vehicle):
        k_i = vehicle.controller.states[k].omega / vehicle.states[k].velocity
        if k_i == 0:
            s_magnitude = 0
        else:
            s_magnitude = (-1 + np.sqrt(1 + (k_i ** 2) * (self.r + self.h * vehicle.states[k].velocity) ** 2)) / k_i
        self.states[k].s_x = s_magnitude * math.sin(vehicle.states[k].theta)
        self.states[k].s_y = -s_magnitude * math.cos(vehicle.states[k].theta)

    def alpha(self, k, prec, vehicle):
        k_i_1 = prec.controller.states[k].omega / prec.states[k].velocity
        res = k_i_1 * (self.r + self.h * vehicle.states[k].velocity)
        self.states[k].alpha = math.atan(res)

    def error(self, k, prec, vehicle):
        self.states[k].error_x = (
                prec.states[k].x + prec.states[k].s_x
                - vehicle.states[k].x
                - (self.r + self.h * vehicle.states[k].velocity) * math.cos(vehicle.states[k].theta)
        )
        self.states[k].error_y = (
                prec.states[k].y + prec.states[k].s_y
                - vehicle.states[k].y
                - (self.r + self.h * vehicle.states[k].velocity) * math.sin(vehicle.states[k].theta)
        )
        self.states[k].error_velocity_x = (
                prec.states[k].velocity * math.cos(prec.states[k].theta)
                - vehicle.states[k].velocity * math.cos(vehicle.states[k].theta + vehicle.controller.states[k].alpha)
        )
        self.states[k].error_velocity_y = (
                prec.states[k].velocity * math.sin(prec.states[k].theta)
                - vehicle.states[k].velocity * math.sin(vehicle.states[k].theta + vehicle.controller.states[k].alpha)
        )

    def updateState(self, k, T, prec, vehicle, a, w):
        pass

    def update_state_init(self, k, prec, vehicle, a, w, T):
        self.states.append(Controller_extend_state())
        if vehicle.first:
            self.states[k].error_x = 0
            self.states[k].error_y = 0
            self.states[k].error_velocity_x = 0
            self.states[k].error_velocity_y = 0
            self.states[k].acceleration = a
            self.states[k].omega = w
            self.states[k].alpha = 0
            self.s(k, vehicle)
        else:
            self.s(k, vehicle)
            self.alpha(k, prec, vehicle)
            self.error(k, prec, vehicle)
            self.get_acceleration_omega(prec, vehicle, k, T)

    def T_12_inv(self, prec, vehicle, k):
        u_i = self.h * (self.r + self.h * vehicle.states[k].velocity) * (
                1 - math.sin(vehicle.controller.states[k].alpha)
                * math.sin(prec.states[k].theta - vehicle.states[k].theta)
        )
        s_a_i = self.h * math.sin(vehicle.controller.states[k].alpha)
        T_12_inv = np.array([
            [(self.r + self.h * vehicle.states[k].velocity) * math.cos(vehicle.states[k].theta) / u_i,
             (self.r + self.h * vehicle.states[k].velocity) * math.sin(vehicle.states[k].theta) / u_i],
            [((-self.h * math.sin(vehicle.states[k].theta)) - (s_a_i * math.cos(prec.states[k].theta))) / u_i,
             ((self.h * math.cos(vehicle.states[k].theta)) - (s_a_i * math.sin(prec.states[k].theta))) / u_i]
        ])
        return T_12_inv

    def B_1(self,prec,vehicle,k,T):
        O=np.array([
            [-math.sin(vehicle.states[k].theta)],
            [math.cos(vehicle.states[k].theta)]
        ])
        R= np.array([
            [math.cos(prec.states[k].theta), -math.sin(prec.states[k].theta)],
            [math.sin(prec.states[k].theta), math.cos(prec.states[k].theta)]
        ])
        k_i = prec.controller.states[k].omega / prec.states[k].velocity
        if k_i == 0:
            s_magnitude = 0
        else:
            s_magnitude = (-1 + np.sqrt(1 + (k_i ** 2) * (self.r + self.h * prec.states[k].velocity) ** 2)) / k_i
        s_k = (1 - math.cos(vehicle.controller.states[k].alpha)) / (k_i ** 2)

        der_k = (prec.controller.states[k].omega / prec.states[k].velocity) -(prec.controller.states[k-1].omega / prec.states[k-1].velocity) / T

        S = np.array([
            [s_magnitude * prec.controller.states[k].omega],
            [-s_k * der_k]
        ])
        B_1 = O *vehicle.states[k].velocity * math.tan(vehicle.controller.states[k].alpha) + np.dot(R,S)
        return B_1

    def get_acceleration_omega(self, prec, vehicle, k, T):
        if k == 0:
            kk = 0
        else:
            kk = k + 1

        Z_12 = np.array([
            [self.k1 * vehicle.controller.states[k].error_x],
            [self.k2 * vehicle.controller.states[k].error_y]
        ])

        Z_34 = np.array([
            [vehicle.controller.states[k].error_velocity_x / math.cos(vehicle.controller.states[k].alpha)],
            [vehicle.controller.states[k].error_velocity_y / math.cos(vehicle.controller.states[k].alpha)]
        ])

        T_12_inv = self.T_12_inv(prec, vehicle, k)
        B_1 = self.B_1(prec,vehicle,k,T)

        Res=Z_12 + Z_34 + B_1
        A=np.dot(T_12_inv,Res)
        self.states[kk].acceleration = A[0, 0]
        self.states[kk].omega = A[1, 0]




