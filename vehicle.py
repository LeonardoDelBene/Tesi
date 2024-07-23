
import math

from Vehicle_State import VehicleState


class Vehicle:
    def __init__(self, first, controller):
        self.first = first
        self.controller = controller
        self.states = []

    def updateState_first(self,k,T,a,w):
        self.states.append(VehicleState(0,0,0,0))
        self.states[k].x = self.states[k-1].x + (T * self.states[k-1].velocity * math.cos(self.states[k-1].theta))
        self.states[k].y = self.states[k-1].y + (T * self.states[k-1].velocity * math.sin(self.states[k-1].theta))
        self.states[k].theta = self.states[k-1].theta + (T * w)
        self.states[k].velocity = self.states[k-1].velocity + (T * a)

    def updateState(self,k,T,prec):
        self.states.append(VehicleState(0, 0, 0, 0))
        A =self.controller.get_acceleration_omega(prec,self,k-1)
        self.states[k].x = self.states[k - 1].x + (T * self.states[k - 1].velocity * math.cos(self.states[k - 1].theta))
        self.states[k].y = self.states[k - 1].y + (T * self.states[k - 1].velocity * math.sin(self.states[k - 1].theta))
        self.states[k].theta = self.states[k - 1].theta + (T *A[1,0] )
        self.states[k].velocity = self.states[k - 1].velocity + (T * A[0,0])
        print("Vehicle " + str(self.first) + " at time " + str(k))
        print("x:" + str(self.states[k].x) + " y:" + str(self.states[k].y) + " theta:" + str(self.states[k].theta) + " velocity:" + str(self.states[k].velocity))

