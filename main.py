from vehicle import *
from Controller import *
from Controller_State import *
from Vehicle_State import *
import matplotlib.pyplot as plt

def main():
    k1 = 3.5
    k2 = 3.5
    r = 1
    h = 0.2
    T = 0.1
    N = 3 # Numero di veicoli
    num_steps = 500
    vehicles = []
    velocity = 5

    # Arrays to store the trajectory for each vehicle
    x_positions = [[] for _ in range(N)]
    y_positions = [[] for _ in range(N)]

    acceleration = []
    omega = []
    for j in range(num_steps):
        acceleration.append(0)
        if j < 59:
            omega.append(0)
        else:
            omega.append(0.5)

    for i in range(N):
        if i == 0:
            vehicles.append(Vehicle(True, Controller(r, h, k1, k2)))
            vehicles[i].states.append(VehicleState(i, i, 0, velocity))
            vehicles[i].controller.updateState_init(0, vehicles[i - 1], vehicles[i])
        else:
            vehicles.append(Vehicle(False, Controller(r, h, k1, k2)))
            vehicles[i].states.append(VehicleState(-i, i, 0, velocity))
            vehicles[i].controller.updateState_init(0, vehicles[i-1], vehicles[i])


    for k in range(1, num_steps):
        for i in range(N):
            if(i == 0):
                vehicles[i].updateState_first(k, T, acceleration[k], omega[k])
            else:
                vehicles[i].updateState(k, T, vehicles[i-1])
            vehicles[i].controller.updateState(k, T, vehicles[i - 1] if i > 0 else None, vehicles[i])
            state = vehicles[i].states[-1]
            x_positions[i].append(state.x)
            y_positions[i].append(state.y)

    # Plotting the trajectories
    plt.figure(figsize=(10, 6))
    for i in range(N):
        plt.plot(x_positions[i], y_positions[i], label=f'Vehicle {i+1} Trajectory')
    plt.title('Vehicle Trajectories')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
