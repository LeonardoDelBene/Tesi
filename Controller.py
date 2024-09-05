from abc import ABC, abstractmethod

class Controller(ABC):
    def __init__(self, r, h, k1, k2):
        self.states = []
        self.r = r
        self.h = h
        self.k1 = k1
        self.k2 = k2

    @abstractmethod
    def get_acceleration_omega(self, prec, vehicle, k,T):
        pass

    @abstractmethod
    def update_state_init(self, k, prec, vehicle, a, w,T):
        pass

    @abstractmethod
    def updateState(self, k, T, prec, vehicle, a, w):
        pass