import math

from Controller_State import ControllerState


class Controller_extend:
    def __init__(self, r, h, k1, k2):
        self.states = []
        self.r = r
        self.h = h
        self.k1 = k1
        self.k2 = k2
        self.s_x_prec=0
        self.s_y_prec=0
        self.alpha=0


    def s_magnitude(self,velocity,velocity_prec, omega_prec):
        k_prec= omega_prec/velocity_prec
        if(k_prec==0):
            return 0
        else:
            return ((-1 + math.sqrt(1+((k_prec**2)*(self.r + self.h*velocity)**2)))/(k_prec))
    def s(self,theta_prec,velocity_prec,omega_prec,velocity):
        s_magnitude_prec = self.s_magnitude(velocity,velocity_prec,omega_prec)
        self.s_x_prec = s_magnitude_prec * math.sin(theta_prec)
        self.s_y_prec = -s_magnitude_prec * math.cos(theta_prec)

    def alpha(self,omega_prec,velocity_prec,velocity):
        k_prec= omega_prec/velocity_prec
        calculated_alpha = math.atan((k_prec*(self.r + self.h*velocity)))

        if calculated_alpha < -math.pi / 2:
            self.alpha = -math.pi / 2
        elif calculated_alpha > math.pi / 2:
            self.alpha = math.pi / 2
        else:
            self.alpha = calculated_alpha

    def error(self,prec,vehicle,k):
        self.states[k].error_x = prec.states[k].x +self.s_x_prec - vehicle.states[k].x - (
                (self.r + (self.h * vehicle.states[k].velocity)) * math.cos(vehicle.states[k].theta))

        self.states[k].error_y = prec.states[k].y + self.s_y_prec - vehicle.states[k].y - (
                (self.r + (self.h * vehicle.states[k].velocity)) * math.sin(vehicle.states[k].theta))

        self.states[k].error_velocity_x = prec.states[k].velocity * math.cos(prec.states[k].theta) - vehicle.states[k].velocity * math.cos(vehicle.states[k].theta + self.alpha(prec.states[k].omega,prec.states[k].velocity,vehicle.states[k].velocity))

        self.states[k].error_velocity_y = prec.states[k].velocity * math.sin(prec.states[k].theta) - vehicle.states[k].velocity * math.sin(vehicle.states[k].theta + self.alpha(prec.states[k].omega,prec.states[k].velocity,vehicle.states[k].velocity))


    def update_state_initi(self,k,prec, vehicle,a,w):
        self.states.append(ControllerState())
        if(vehicle.first):
            self.states[k].acceleration = a
            self.states[k].omega = w
            self.states[k].error_x = 0
            self.states[k].error_y = 0
            self.states[k].error_velocity_x = 0
            self.states[k].error_velocity_y = 0
        else:
            self.s(prec.states[k].theta,prec.states[k].velocity,prec.states[k].omega,vehicle.states[k].velocity)
            self.alpha(prec.states[k].omega,prec.states[k].velocity,vehicle.states[k].velocity)
            self.error(prec,vehicle,k)

