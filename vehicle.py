
import math

from Vehicle_State import VehicleState


class Vehicle:
    def __init__(self, first, controller):
        self.first = first
        self.controller = controller
        self.states = []

    def updateState(self, k, T, prec):
        self.states.append(VehicleState(0, 0, 0, 0))
        if (prec != None):
            A = self.controller.get_acceleration_omega(prec, self, k - 1)
            self.states[k].x = self.states[k - 1].x + (
                    T * self.states[k - 1].velocity * math.cos(self.states[k - 1].theta))
            self.states[k].y = self.states[k - 1].y + (
                    T * self.states[k - 1].velocity * math.sin(self.states[k - 1].theta))
            self.states[k].velocity = self.states[k - 1].velocity + (T * A[0, 0])
            self.states[k].theta = self.states[k - 1].theta + (T * A[1, 0])

        else:
            self.states[k].x = self.states[k - 1].x + (
                    T * self.states[k - 1].velocity * math.cos(self.states[k - 1].theta))
            self.states[k].y = self.states[k - 1].y + (
                    T * self.states[k - 1].velocity * math.sin(self.states[k - 1].theta))
            self.states[k].velocity = self.states[k - 1].velocity + (T * self.controller.states[k - 1].acceleration)
            self.states[k].theta = self.states[k - 1].theta + (T * self.controller.states[k - 1].omega)


