import math
import numpy as np

from Controller import Controller
from Controller_extend_state import Controller_extend_state


class Controller_extend(Controller):
    def __init__(self, r, h, k1, k2):
        super().__init__(r, h, k1, k2)
        self.epsilon = 1e-8  # Small value to prevent division by zero

    def s_magnitude(self, k, vehicle, prec):
        kk = max(k - 1, 0)  # Safeguard for k-1 when k=0
        velocity = max(vehicle.states[kk].velocity, self.epsilon)
        velocity_1 = max(prec.states[kk].velocity, self.epsilon)
        k_i = prec.controller.states[kk].omega / velocity_1


        if k_i == 0:
            s_magnitude = 0
            sqrt_val = 0
        else:
            sqrt_val = max(1 + ((k_i ** 2) * (self.r + self.h * velocity) ** 2), 0)  # Ensure non-negative for sqrt
            s_magnitude = (-1 + math.sqrt(sqrt_val)) / k_i
        return s_magnitude

    def Alpha(self, k, prec, vehicle):
        kk = max(k - 1, 0)
        velocity = max(prec.states[kk].velocity, self.epsilon)  # Prevent zero velocity
        k_i_1 = prec.controller.states[kk].omega / velocity
        res = k_i_1 * (self.r + (self.h * vehicle.states[kk].velocity))
        self.states[k].alpha = math.atan(res)


    def error(self, k, prec, vehicle):
        s_magnitude=self.s_magnitude(k,vehicle, prec)
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
            self.Alpha(k, prec, vehicle)

            self.states[k].error_x = self.states[k - 1].error_x + (T * (-self.k1 * self.states[k - 1].error_x))
            self.states[k].error_y = self.states[k - 1].error_y + (T * (-self.k2 * self.states[k - 1].error_y))

            kk = k - 1
            H_i_1 = np.array([
                [math.cos(prec.states[kk].theta), -prec.states[kk].velocity * math.sin(prec.states[kk].theta)],
                [math.sin(prec.states[kk].theta), prec.states[kk].velocity * math.cos(prec.states[kk].theta)]
            ])

            A_i_1 = np.array([
                [prec.controller.states[kk].acceleration],
                [prec.controller.states[kk].omega]
            ])

            # Controllo per k_i_1 per evitare divisioni per zero
            k_i_1 = prec.controller.states[kk].omega / max(prec.states[kk].velocity, self.epsilon)
            delta_i = vehicle.states[kk].velocity * self.h * (math.cos(vehicle.controller.states[kk].alpha) ** 2) * k_i_1
            theta_alpha = vehicle.states[kk].theta + vehicle.controller.states[kk].alpha

            T_34_i = np.array([
                [math.cos(theta_alpha) - (delta_i * math.sin(theta_alpha)),
                 -vehicle.states[kk].velocity * math.sin(theta_alpha)],
                [math.sin(theta_alpha) + (delta_i * math.cos(theta_alpha)),
                 vehicle.states[kk].velocity * math.cos(theta_alpha)]
            ])

            K = np.array([
                [self.k1 * self.states[kk].error_x],
                [self.k2 * self.states[kk].error_y]
            ])

            cos_alpha = max(math.cos(vehicle.controller.states[kk].alpha), self.epsilon)


            Z_34 = np.array([
                [self.states[kk].error_velocity_x / cos_alpha],
                [self.states[kk].error_velocity_y / cos_alpha]
            ])

            B_1 = self.B_1(prec, vehicle, kk, T)
            T_12_inv = self.T_12_inv(prec, vehicle, kk)

            der_k = ((prec.controller.states[kk].omega / max(prec.states[kk].velocity, self.epsilon)) -
                     (prec.controller.states[kk - 1].omega / max(prec.states[kk - 1].velocity, self.epsilon))) / T

            # Prevenzione dell'overflow in B_2
            B_2_factor = self.r + (self.h * vehicle.states[kk].velocity)
            O_i = np.array([
                [math.sin(theta_alpha)],
                [-math.cos(theta_alpha)]
            ])
            B_2 = np.clip(B_2_factor * vehicle.states[kk].velocity * (cos_alpha ** 2)* O_i * der_k, -1e10, 1e10)

            Res = np.dot(T_34_i, T_12_inv)
            Res1 = (Z_34 + B_1)
            Res2 = np.dot(Res, K)
            Res3 = np.dot(H_i_1, A_i_1)

            # Prevenzione di valori invalidi o overflow nel risultato finale
            Result = np.clip(Res3 + B_2 - Res2 - (np.dot(Res, Res1)), -1e10, 1e10)

            self.states[k].error_velocity_x = self.states[k - 1].error_velocity_x + (T * Result[0, 0])
            self.states[k].error_velocity_y = self.states[k - 1].error_velocity_y + (T * Result[1, 0])

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
        u_i = self.h * (self.r + (self.h * vehicle.states[k].velocity)) * max(u_i_denom, self.epsilon)

        s_a_i = self.h * math.sin(vehicle.controller.states[k].alpha)
        T_12_inv = np.array([
            [(self.r + (self.h * vehicle.states[k].velocity)) * math.cos(vehicle.states[k].theta) / u_i,
             (self.r + (self.h * vehicle.states[k].velocity)) * math.sin(vehicle.states[k].theta) / u_i],
            [((-self.h * math.sin(vehicle.states[k].theta)) - (s_a_i * math.cos(prec.states[k].theta))) / u_i,
             ((self.h * math.cos(vehicle.states[k].theta)) - (s_a_i * math.sin(prec.states[k].theta))) / u_i]
        ])
        return T_12_inv

    def B_1(self, prec, vehicle, k, T):
        epsilon = self.epsilon  # Costante per prevenire divisioni per zero
        max_angle = np.pi / 2 - 1e-5  # Limite per evitare overflow nella tangente

        O = np.array([
            [-math.sin(vehicle.states[k].theta)],
            [math.cos(vehicle.states[k].theta)]
        ])

        R = np.array([
            [math.cos(prec.states[k].theta), -math.sin(prec.states[k].theta)],
            [math.sin(prec.states[k].theta), math.cos(prec.states[k].theta)]
        ])

        velocity_1 = max(prec.states[k].velocity, epsilon)

        k_i = prec.controller.states[k].omega / velocity_1
        der_k = ((prec.controller.states[k].omega / velocity_1) - (prec.controller.states[k - 1].omega / max(prec.states[k - 1].velocity, epsilon))) / T

        s_magnitude = self.s_magnitude(k, vehicle, prec)

        if abs(k_i) < epsilon:
            s_k = 0
        else:
            cos_alpha = max(min(math.cos(vehicle.controller.states[k].alpha), 1.0), -1.0)
            s_k = -(1 - cos_alpha) / (k_i ** 2)

        S = np.array([
            [s_magnitude * prec.controller.states[k].omega],
            [s_k * der_k]
        ])

        O_1 = np.array([
            [math.cos(prec.states[k].theta)],
            [math.sin(prec.states[k].theta)]
        ])

        # Limita l'angolo per prevenire valori troppo grandi nella tangente
        alpha = vehicle.controller.states[k].alpha
        alpha = max(min(alpha, max_angle), -max_angle)

        cos_alpha_safe = max(min(1 / math.cos(alpha), 1e10), -1e10)
        B_1 = (O * (vehicle.states[k].velocity * math.tan(alpha)) +
               np.dot(R, S) +
               ((1 - cos_alpha_safe) * O_1 * prec.states[k].velocity))

        return B_1

    def get_acceleration_omega(self, prec, vehicle, k, T):
        # Se k == 0, utilizza kk = 0, altrimenti kk = k + 1
        kk = 0 if k == 0 else k + 1

        Z_12 = np.array([
            [self.k1 * vehicle.controller.states[k].error_x],
            [self.k2 * vehicle.controller.states[k].error_y]
        ])

        # Utilizzo di epsilon per evitare divisioni per zero nel calcolo di Z_34
        cos_alpha = max(math.cos(vehicle.controller.states[k].alpha), self.epsilon)  # Prevenzione di cos(alpha) troppo piccolo

        Z_34 = np.array([
            [vehicle.controller.states[k].error_velocity_x / cos_alpha],
            [vehicle.controller.states[k].error_velocity_y / cos_alpha]
        ])

        T_12_inv = self.T_12_inv(prec, vehicle, k)

        B_1 = self.B_1(prec, vehicle, k, T)

        Res = Z_12 + Z_34 + B_1

        A = np.dot(T_12_inv, Res)

        self.states[kk].acceleration = A[0, 0]
        self.states[kk].omega = A[1, 0]

        return A





