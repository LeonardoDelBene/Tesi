from vehicle import *
from Controller_standard import *
from Vehicle_State import *
import matplotlib.pyplot as plt


def main():
    k1 = 3
    k2 = 3
    r = 1
    h = 0.2
    T = 0.01 # Passo di campionameneto
    N = 2 # Numero di veicoli
    num_steps = 2000 # Numero di passi 20 secondi
    vehicles = []
    velocity = 5

    x_positions = [[] for _ in range(N)]
    y_positions = [[] for _ in range(N)]

    acceleration = []
    omega = []
    for j in range(num_steps):
        acceleration.append(0)
        if j < 1/T:
            omega.append(0)
        elif j >= 1/T and j < 4/T:
            omega.append(0.5)
        elif j >= 4/T and j < 6/T:
            omega.append(0)
        elif j >= 6/T and j < 8/T:
            omega.append(-0.5)
        elif j >= 8/T and j < 10/T:
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
