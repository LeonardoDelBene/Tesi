import math
import numpy as np

class ControllerState:
    def __init__(self):
        self.error_x = 0
        self.error_y = 0
        self.error_velocity_x = 0
        self.error_velocity_y = 0
        self.acceleration = 0
        self.omega = 0
