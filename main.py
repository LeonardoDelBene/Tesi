from vehicle import *
from Controller import *
from Controller_State import *
from Vehicle_State import *
import matplotlib.pyplot as plt

def main():
    k1 = 0.25
    k2 = 0.25
    r = 1
    h = 0.2
    T = 0.1
    N = 2 # Numero di veicoli
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
        else:
            vehicles.append(Vehicle(False, Controller(r, h, k1, k2)))
        vehicles[i].states.append(VehicleState(-i, i, 0, velocity))
        if i != 0:
            previous_vehicle = vehicles[i - 1]
        else:
            previous_vehicle = None
        vehicles[i].controller.updateState_init(0, previous_vehicle, vehicles[i], acceleration[0], omega[0])


    for k in range(1, num_steps):
        for i in range(N):
            vehicles[i].controller.updateState(k, T, vehicles[i - 1] if i > 0 else None, vehicles[i], acceleration[k], omega[k])
            vehicles[i].updateState(k, T, vehicles[i - 1] if i > 0 else None)
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
