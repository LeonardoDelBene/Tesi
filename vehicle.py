
import math

import numpy as np

from Vehicle_State import VehicleState


class Vehicle:
    def __init__(self, first, controller,id):
        self.first = first
        self.controller = controller
        self.states = []
        self.id=id

    def updateState(self, k, T, prec):
        self.states.append(VehicleState(0, 0, 0, 0))
        if (prec != None):
            A = self.controller.get_acceleration_omega(prec, self, k-1,T)
            self.states[k].x = self.states[k - 1].x + (T * self.states[k - 1].velocity * math.cos(self.states[k - 1].theta))
            self.states[k].y = self.states[k - 1].y + (T * self.states[k - 1].velocity * math.sin(self.states[k - 1].theta))
            vel = self.states[k - 1].velocity + (T * A[0, 0])
            vel =np.clip(vel,-3,15)
            self.states[k].velocity = vel

            theta = self.states[k - 1].theta + (T * A[1, 0])
            theta = np.clip(theta, -math.pi * 2, math.pi * 2)
            self.states[k].theta = theta

        else:
            self.states[k].x = self.states[k - 1].x + (
                    T * self.states[k - 1].velocity * math.cos(self.states[k - 1].theta))
            self.states[k].y = self.states[k - 1].y + (
                    T * self.states[k - 1].velocity * math.sin(self.states[k - 1].theta))
            self.states[k].velocity = self.states[k - 1].velocity + (T * self.controller.states[k - 1].acceleration)
            self.states[k].theta = self.states[k - 1].theta + (T * self.controller.states[k - 1].omega)

