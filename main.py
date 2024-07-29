import math

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Controller_extend import Controller_extend
from vehicle import *
from Controller_standard import *
from Vehicle_State import *



def Error(vehicle, prec, T_max, T, h, v, r,vehicle_number):
    error = []
    times = np.arange(0, T_max, T)
    num_steps = len(times)

    # Prepare lists for positions of the current vehicle
    Pos_vehicle = [(vehicle.states[i].x, vehicle.states[i].y) for i in range(num_steps)]

    # Calculate positions of `prec` at time t - h - (r / v)
    Pos_prec = []
    for i in range(num_steps):
        # Time index for t - h - (r/v)
        t_minus_h_minus_v_r = i - int(h + (r/v[i])/T)
        # Ensure index is within bounds
        if 0 <= t_minus_h_minus_v_r < len(prec.states):
            Pos_prec.append((prec.states[t_minus_h_minus_v_r].x, prec.states[t_minus_h_minus_v_r].y))
        else:
            # If the index is out of bounds, append None or handle accordingly
            Pos_prec.append((None, None))

    # Compute the error
    for i in range(num_steps):
        if Pos_prec[i] == (None, None):
            error.append(np.nan)  # Use NaN to indicate an invalid error value
        else:
            diff = np.array(Pos_vehicle[i]) - np.array(Pos_prec[i])
            error.append(np.linalg.norm(diff))

    # Print debug information
    print("Times:", times)
    print("Error values:", error)

    # Plot the error vs time
    plt.figure(figsize=(10, 6))
    plt.plot(times, error, label=f'Error between Vehicle {vehicle_number} and Vehicle {vehicle_number - 1}')
    plt.title(f'Error between Vehicle {vehicle_number} and Vehicle {vehicle_number -1} Positions over Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Error (L2 Norm)')
    plt.legend()
    plt.grid(True)
    plt.show()

    return error

def main():
    k1 = 2.5
    k2 = 2.5
    r = 1
    h = 0.2
    T = 0.01  # Passo di campionamento
    N = 5  # Numero di veicoli
    T_max = 20  # Tempo massimo
    num_steps = int(T_max/T)  # Numero di passi
    vehicles = []
    velocity = 5

    x_positions = [[] for _ in range(N)]
    y_positions = [[] for _ in range(N)]
    times = [k * T for k in range(num_steps)]

    acceleration = []
    omega = []
    for j in range(num_steps):
        acceleration.append(0)
        if j < 7 / T:
            omega.append(0)
        else:
            omega.append(0.5)

    contr = input("Inserire 1 per Controller_standard, 2 per Controller_extend: ")
    if contr == "1":
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
    elif contr == "2":
        for i in range(N):
            if i == 0:
                vehicles.append(Vehicle(True, Controller_extend(r, h, k1, k2)))
            else:
                vehicles.append(Vehicle(False, Controller_extend(r, h, k1, k2)))
            vehicles[i].states.append(VehicleState(-i, i, 0, velocity))
            if i != 0:
                previous_vehicle = vehicles[i - 1]
            else:
                previous_vehicle = None
            vehicles[i].controller.update_state_init(0, previous_vehicle, vehicles[i], acceleration[0], omega[0])

    # Append initial positions
    for i in range(N):
        state = vehicles[i].states[-1]
        x_positions[i].append(state.x)
        y_positions[i].append(state.y)

    for k in range(1, num_steps):
        for i in range(N):
            vehicles[i].controller.updateState(k, T, vehicles[i - 1] if i > 0 else None, vehicles[i], acceleration[k],
                                               omega[k])
            vehicles[i].updateState(k, T, vehicles[i - 1] if i > 0 else None)
            state = vehicles[i].states[-1]
            x_positions[i].append(state.x)
            y_positions[i].append(state.y)

    # 2D Plots
    plt.figure(figsize=(15, 10))

    # Plot vehicle trajectories
    plt.subplot(3, 1, 1)
    for i in range(N):
        plt.plot(x_positions[i], y_positions[i], label=f'Vehicle {i + 1} Trajectory')
    plt.title('Vehicle Trajectories')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.legend()
    plt.grid(True)

    # Plot X positions vs. time
    plt.subplot(3, 1, 2)
    for i in range(N):
        plt.plot(times, x_positions[i], label=f'Vehicle {i + 1} X Position')
    plt.title('X Position vs Time')
    plt.xlabel('Time (s)')
    plt.ylabel('X Position')
    plt.legend()
    plt.grid(True)

    # Plot Y positions vs. time
    plt.subplot(3, 1, 3)
    for i in range(N):
        plt.plot(times, y_positions[i], label=f'Vehicle {i + 1} Y Position')
    plt.title('Y Position vs Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Y Position')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

    # 3D Plot: X, Y positions vs. time
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    for i in range(N):
        ax.plot(times, x_positions[i], y_positions[i], label=f'Vehicle {i + 1}')
    ax.set_title('3D Trajectory of Vehicles')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('X Position')
    ax.set_zlabel('Y Position')
    ax.legend()
    ax.grid(True)

    plt.show()


    # Calculate and plot errors between vehicles
    v = []
    Error_total = []
    for i in range(1, N):
        for j in range(num_steps):
            v.append(vehicles[i - 1].states[j].velocity)
        e = Error(vehicles[i], vehicles[i - 1], T_max, T, h, v, r, i)
        Error_total.append(e)

    Error_medio = []
    for j in range(1, N):
        # Inizializza la lista 'sum' con zeri
        sum_errors = [0] * num_steps

        # Calcola la somma dei quadrati degli errori
        for i in range(num_steps):
            if np.isnan(Error_total[j - 1][i]):
                sum_errors[i] = sum_errors[i - 1] if i > 0 else 0
            else:
                if i == 0:
                    sum_errors[i] = Error_total[j - 1][i] ** 2
                else:
                    sum_errors[i] = sum_errors[i - 1] + (Error_total[j - 1][i] ** 2)

        average_sum = sum_errors[num_steps - 1] / num_steps
        Error_medio.append(math.sqrt(average_sum))

        print(f"Errore medio tra i veicoli {j - 1 + 1} e {j + 1} : {Error_medio[j - 1]}")










if __name__ == "__main__":
    main()
