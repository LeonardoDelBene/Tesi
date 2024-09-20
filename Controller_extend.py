import math
import numpy as np

from Controller import Controller
from Controller_extend_state import Controller_extend_state


class Controller_extend(Controller):
    def __init__(self, r, h, k1, k2):
        super().__init__(r, h, k1, k2)


    def s_magnitude(self, k, vehicle, prec):
        kk = max(k - 1, 0)  # Safeguard for k-1 when k=0
        velocity = vehicle.states[kk].velocity
        velocity_1 = prec.states[kk].velocity
        k_i = prec.controller.states[kk].omega / velocity_1

        if k_i == 0:
            s_magnitude = 0
        else:
            sqrt_val = max(1 + ((k_i ** 2) * (self.r + self.h * velocity) ** 2), 0)  # Ensure non-negative for sqrt
            s_magnitude = (-1 + math.sqrt(sqrt_val)) / k_i
        return s_magnitude

    def Alpha(self, k, prec, vehicle):
        kk = max(k - 1, 0)
        velocity = prec.states[kk].velocity
        k_i_1 = prec.controller.states[kk].omega / velocity
        res = k_i_1 * (self.r + (self.h * vehicle.states[kk].velocity))
        self.states[k].alpha = math.atan(res)

    def error(self, k, prec, vehicle):
        s_magnitude = self.s_magnitude(k,vehicle, prec)
        self.states[k].error_x = (
            prec.states[k].x + (s_magnitude * math.sin(prec.states[k].theta))
            - vehicle.states[k].x
            - (self.r + (self.h * vehicle.states[k].velocity)) * math.cos(vehicle.states[k].theta)
        )
        self.states[k].error_y = (
            prec.states[k].y + (-s_magnitude * math.cos(prec.states[k].theta))
            - vehicle.states[k].y
            - (self.r + (self.h * vehicle.states[k].velocity)) * math.sin(vehicle.states[k].theta)
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
        self.states.append(Controller_extend_state())
        if vehicle.first:
            self.states[k].error_x = 0
            self.states[k].error_y = 0
            self.states[k].error_velocity_x = 0
            self.states[k].error_velocity_y = 0
            self.states[k].acceleration = a
            self.states[k].omega = w
            self.states[k].alpha = 0
        else:
            if (vehicle.id == 1 and k < 1510):
                print("vehicle: ", vehicle.id, " - Iteration: ", k)
            self.Alpha(k, prec, vehicle)
            if (vehicle.id == 1 and k < 1510):
                print(f"Alpha: {self.states[k].alpha}")

            self.states[k].error_x = self.states[k - 1].error_x + (T * (-self.k1 * self.states[k - 1].error_x))
            self.states[k].error_y = self.states[k - 1].error_y + (T * (-self.k2 * self.states[k - 1].error_y))

            kk = k - 1
            H_i_1 = np.array([
                [math.cos(prec.states[kk].theta), -prec.states[kk].velocity * math.sin(prec.states[kk].theta)],
                [math.sin(prec.states[kk].theta), prec.states[kk].velocity * math.cos(prec.states[kk].theta)]
            ])
            if (vehicle.id == 1 and k < 1510):
                print(f"H_i_1: {H_i_1}")

            A_i_1 = np.array([
                [prec.controller.states[kk].acceleration],
                [prec.controller.states[kk].omega]
            ])
            if (vehicle.id == 1 and k < 1510):
                print(f"A_i_1: {A_i_1}")


            k_i_1 = prec.controller.states[kk].omega / prec.states[kk].velocity
            delta_i = vehicle.states[kk].velocity * self.h * (math.cos(vehicle.controller.states[kk].alpha) ** 2) * k_i_1
            if (vehicle.id == 1 and k < 1510):
                print(f"Delta_i: {delta_i}")
            theta_alpha = vehicle.states[kk].theta + vehicle.controller.states[kk].alpha

            T_34_i = np.array([
                [math.cos(theta_alpha) - (delta_i * math.sin(theta_alpha)),
                 -vehicle.states[kk].velocity * math.sin(theta_alpha)],
                [math.sin(theta_alpha) + (delta_i * math.cos(theta_alpha)),
                 vehicle.states[kk].velocity * math.cos(theta_alpha)]
            ])
            if (vehicle.id == 1 and k < 1510):
                print(f"T_34_i: {T_34_i}")

            K = np.array([
                [self.k1 * self.states[kk].error_x],
                [self.k2 * self.states[kk].error_y]
            ])
            if (vehicle.id == 1 and k < 1510):
                print(f"K: {K}")

            cos_alpha = math.cos(vehicle.controller.states[kk].alpha)
            if (vehicle.id == 1 and k < 1510):
                print(f"Cos_alpha: {cos_alpha}")


            Z_34 = np.array([
                [self.states[kk].error_velocity_x / cos_alpha],
                [self.states[kk].error_velocity_y / cos_alpha]
            ])
            if (vehicle.id == 1 and k < 1510):
                print(f"Z_34: {Z_34}")

            B_1 = self.B_1(prec, vehicle, kk, T)
            if (vehicle.id == 1 and k < 1510):
                print(f"B_1: {B_1}")
            T_12_inv = self.T_12_inv(prec, vehicle, kk)
            if (vehicle.id == 1 and k < 1510):
                print(f"T_12_inv: {T_12_inv}")

            der_k = ((prec.controller.states[kk].omega / prec.states[kk].velocity) -
                     (prec.controller.states[kk - 1].omega / prec.states[kk - 1].velocity)) / T
            #der_k = np.clip(der_k, -5,5)
            if (vehicle.id == 1 and k < 1510):
                print(f"Derivative of k in B2: {der_k}")


            B_2_factor = self.r + (self.h * vehicle.states[kk].velocity)
            O_i = np.array([
                [math.sin(theta_alpha)],
                [-math.cos(theta_alpha)]
            ])
            if (vehicle.id == 1 and k < 1510):
                print(f"O_i: {O_i}")
            B_2 = B_2_factor * vehicle.states[kk].velocity * (cos_alpha ** 2)* O_i * der_k
            if (vehicle.id == 1 and k < 1510):
                print(f"B_2: {B_2}")

            Res = np.dot(T_34_i, T_12_inv)
            if (vehicle.id == 1 and k < 1510):
                print(f"T_34 * T_12: {Res}")
            Res1 = (Z_34 + B_1)
            if (vehicle.id == 1 and k < 1510):
                print(f"Z_34 + B_1: {Res1}")
            Res2 = np.dot(Res, K)
            if (vehicle.id == 1 and k < 1510):
                print(f"T_34 * T_12 * K: {Res2}")
            Res3 = np.dot(H_i_1, A_i_1)
            if (vehicle.id == 1 and k < 1510):
                print(f"H_i_1 * A_i_1: {Res3}")

            Result = Res3 + B_2 - Res2 - (np.dot(Res, Res1))
            if (vehicle.id == 1 and k < 1510):
                print(f"Result: {Result}")

            self.states[k].error_velocity_x = self.states[k - 1].error_velocity_x + (T * Result[0, 0])
            self.states[k].error_velocity_y = self.states[k - 1].error_velocity_y + (T * Result[1, 0])
            if (vehicle.id == 1 and k < 1510):
                print(f"Error_x: {self.states[k].error_x}, Error_y: {self.states[k].error_y}, Error_velocity_x: {self.states[k].error_velocity_x}, Error_velocity_y: {self.states[k].error_velocity_y}")

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
        else:
            self.states[k].alpha = 0
            self.error(k, prec, vehicle)
            self.get_acceleration_omega(prec, vehicle, k, T)

    def T_12_inv(self, prec, vehicle, k):
        u_i_denom = (1 - (math.sin(vehicle.controller.states[k].alpha)
                     * math.sin(prec.states[k].theta - vehicle.states[k].theta)))
        u_i = self.h * (self.r + (self.h * vehicle.states[k].velocity)) * u_i_denom
        #if (vehicle.id != 2 and k < 1510):
         #   print(f"U_i: {u_i}")


        s_a_i = self.h * math.sin(vehicle.controller.states[k].alpha)
        #if (vehicle.id != 2 and k < 1510):
         #   print(f"S_a_i: {s_a_i}")
        T_12_inv = np.array([
            [(self.r + (self.h * vehicle.states[k].velocity)) * math.cos(vehicle.states[k].theta) / u_i,
             (self.r + (self.h * vehicle.states[k].velocity)) * math.sin(vehicle.states[k].theta) / u_i],
            [((-self.h * math.sin(vehicle.states[k].theta)) - (s_a_i * math.cos(prec.states[k].theta))) / u_i,
             ((self.h * math.cos(vehicle.states[k].theta)) - (s_a_i * math.sin(prec.states[k].theta))) / u_i]
        ])
        #print("T_12_inv: ", T_12_inv)
        return T_12_inv

    def B_1(self, prec, vehicle, k, T):

        O = np.array([
            [-math.sin(vehicle.states[k].theta)],
            [math.cos(vehicle.states[k].theta)]
        ])

        R = np.array([
            [math.cos(prec.states[k].theta), -math.sin(prec.states[k].theta)],
            [math.sin(prec.states[k].theta), math.cos(prec.states[k].theta)]
        ])

        velocity_1 = prec.states[k].velocity
        k_i = prec.controller.states[k].omega / velocity_1
        if (vehicle.id == 1 and k < 1510):
           print(f"K_i: {k_i}")
        der_k = ((prec.controller.states[k].omega / velocity_1) - (prec.controller.states[k - 1].omega / prec.states[k-1].velocity)) / T
        der_k = np.clip(der_k, -5, 5)
        if (vehicle.id == 1 and k < 1510):
            print(f"Derivative of k in B1: {der_k}")

        s_magnitude = self.s_magnitude(k, vehicle, prec)
        if (vehicle.id == 1 and k < 1510):
            print(f"S_magnitude: {s_magnitude}")
        if k_i == 0:
            s_k = 0
        else:
            cos_alpha = math.cos(vehicle.controller.states[k].alpha)
            s_k = -(1 - cos_alpha) / (k_i ** 2)

        S = np.array([
            [s_magnitude * prec.controller.states[k].omega],
            [s_k * der_k]
        ])
        if (vehicle.id == 1 and k < 1510):
            print(f"S: {S}")

        O_1 = np.array([
            [math.cos(prec.states[k].theta)],
            [math.sin(prec.states[k].theta)]
        ])


        alpha = vehicle.controller.states[k].alpha
        cos_alpha =  math.cos(alpha)
        B_1 = (O * (vehicle.states[k].velocity * math.tan(alpha)) +
               np.dot(R, S) +
               ((1 - cos_alpha) * O_1 * prec.states[k].velocity))
        #if (vehicle.id != 2 and k < 1510):
         #   print(f"B_1: {B_1}")
        return B_1

    def get_acceleration_omega(self, prec, vehicle, k, T):
        # Se k == 0, utilizza kk = 0, altrimenti kk = k + 1
        kk = 0 if k == 0 else k + 1

        Z_12 = np.array([
            [self.k1 * vehicle.controller.states[k].error_x],
            [self.k2 * vehicle.controller.states[k].error_y]
        ])

        cos_alpha = math.cos(vehicle.controller.states[k].alpha)


        Z_34 = np.array([
            [vehicle.controller.states[k].error_velocity_x / cos_alpha],
            [vehicle.controller.states[k].error_velocity_y / cos_alpha]
        ])


        T_12_inv = self.T_12_inv(prec, vehicle, k)

        B_1 = self.B_1(prec, vehicle, k, T)

        Res = Z_12 + Z_34 + B_1


        A = np.dot(T_12_inv, Res)
        acc= np.clip(A[0, 0], -3, 3)
        ome = np.clip(A[1, 0], -  math.pi/4, math.pi/4)
        #acc = A[0, 0]
        #ome = A[1, 0]

        self.states[kk].acceleration = acc
        self.states[kk].omega = ome
        if (vehicle.id == 1 and k < 1510):
            print(" Acceleration:", {self.states[kk].acceleration}," Omega:", {self.states[kk].omega})

        return A





