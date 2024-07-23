import math
import numpy as np

from Controller_State import ControllerState


class Controller:
    def __init__(self, r, h, k1, k2):
        self.states = []
        self.r = r
        self.h = h
        self.k1 = k1
        self.k2 = k2

    def z1(self, prec, k, vehicle):
        self.states[k].error_x = prec.states[k].x - vehicle.states[k].x - (
                (self.r + (self.h * vehicle.states[k].velocity)) * math.cos(vehicle.states[k].theta))

    def z2(self, prec, k, vehicle):
        self.states[k].error_y = prec.states[k].y - vehicle.states[k].y - (
                (self.r + (self.h * vehicle.states[k].velocity)) * math.sin(vehicle.states[k].theta))

    def z3(self, prec, k, vehicle):
        self.states[k].error_velocity_x = (prec.states[k].velocity * math.cos(prec.states[k].theta)) - (
                vehicle.states[k].velocity * math.cos(vehicle.states[k].theta))
        return self.states[k].error_velocity_x

    def z4(self, prec, k, vehicle):
        self.states[k].error_velocity_y = (prec.states[k].velocity * math.sin(prec.states[k].theta)) - (
                vehicle.states[k].velocity * math.sin(vehicle.states[k].theta))
        return self.states[k].error_velocity_y

    def get_acceleration_omega(self, prec, vehicle, k):
        theta = vehicle.states[k].theta
        velocity = vehicle.states[k].velocity

        # Costruisci la matrice F
        F = np.array([
            [self.h * math.cos(theta),
             -((self.r + (self.h * velocity)) * math.sin(theta))],
            [self.h * math.sin(theta),
             ((self.r + (self.h * velocity)) * math.cos(theta))]
        ])

        # Inverso della matrice F
        F_inv = np.linalg.inv(F)

        # Calcola il vettore degli errori
        E = np.array([
            [self.states[k].error_velocity_x + (self.k1 * self.states[k].error_x)],
            [self.states[k].error_velocity_y + (self.k2 * self.states[k].error_y)]
        ])

        # Calcola U
        U = np.dot(F_inv, E)
        vehicle.controller.states[k].acceleration = U[0, 0]
        vehicle.controller.states[k].omega = U[1, 0]
        return U

    def updateState_init(self, k, prec, vehicle,a,w):
        self.states.append(ControllerState())
        if (prec != None):
            self.z1(prec, k, vehicle)
            self.z2(prec, k, vehicle)
            self.z3(prec, k, vehicle)
            self.z4(prec, k, vehicle)
        if(vehicle.first):
            self.states[k].acceleration = a
            self.states[k].omega = w

    def updateState(self, k, T, prec, vehicle,a,w):
        self.states.append(ControllerState())
        if (prec != None):
            self.states[k].error_x = self.states[k - 1].error_x + (T * (-self.k1 * self.states[k - 1].error_x))
            self.states[k].error_y = self.states[k - 1].error_y + (T * (-self.k2 * self.states[k - 1].error_y))

            G = np.array([
                [(self.h * vehicle.states[k - 1].velocity + self.r * (math.cos(vehicle.states[k - 1].theta) ** 2)) / (
                        self.h * (self.r + self.h * vehicle.states[k - 1].velocity)),
                 (self.r * math.sin(vehicle.states[k - 1].theta) * math.cos(vehicle.states[k - 1].theta)) / (
                         self.h * (self.r + self.h * vehicle.states[k - 1].velocity))],
                [(self.r * math.sin(vehicle.states[k - 1].theta) * math.cos(vehicle.states[k - 1].theta)) / (
                        self.h * (self.r + self.h * vehicle.states[k - 1].velocity)),
                 (self.h * vehicle.states[k - 1].velocity + self.r * (
                         math.sin(vehicle.states[k - 1].theta) ** 2)) / (
                         self.h * (self.r + self.h * vehicle.states[k - 1].velocity))]
            ])

            Z = np.array([
                [self.states[k - 1].error_velocity_x],
                [self.states[k - 1].error_velocity_y]
            ])

            H = np.array([
                [math.cos(prec.states[k - 1].theta),
                 -(prec.states[k - 1].velocity * math.sin(prec.states[k - 1].theta))],
                [math.sin(prec.states[k - 1].theta),
                 (prec.states[k - 1].velocity * math.cos(prec.states[k - 1].theta))]
            ])

            A = np.array([
                [prec.controller.states[k - 1].acceleration],
                [prec.controller.states[k - 1].omega]
            ])

            K = np.array([
                [self.k1 * self.states[k - 1].error_x],
                [self.k2 * self.states[k - 1].error_y]
            ])

            Result = (np.dot(-G, Z)) + ((np.dot(H, A) - np.dot(G, K)))

            self.states[k].error_velocity_x = self.states[k - 1].error_velocity_x + (T * Result[0, 0])
            self.states[k].error_velocity_y = self.states[k - 1].error_velocity_y + (T * Result[1, 0])

        else:
            self.states[k].error_x = 0
            self.states[k].error_y = 0
            self.states[k].error_velocity_x = 0
            self.states[k].error_velocity_y = 0
            self.states[k].acceleration = a
            self.states[k].omega = w
