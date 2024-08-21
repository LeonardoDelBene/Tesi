from Controller_extend import Controller_extend
from vehicle import *
from Controller_standard import *
from Vehicle_State import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.patches as patches



def Error(vehicle, prec, T_max, T, h, v, r,vehicle_number):
    error = []
    times = np.arange(0, T_max, T)
    num_steps = len(times)

    Pos_vehicle = [(vehicle.states[i].x, vehicle.states[i].y) for i in range(num_steps)]

    # Calcolo della posizione del veicolo precedente al tempo t-h-r/v(t)
    Pos_prec = []
    for i in range(num_steps):
        t_minus_h_minus_v_r = i - int(h + (r/v[i])/T)
        if 0 <= t_minus_h_minus_v_r < len(prec.states):
            Pos_prec.append((prec.states[t_minus_h_minus_v_r].x, prec.states[t_minus_h_minus_v_r].y))
        else:
            Pos_prec.append((None, None))

    # Calcolo dell'errore tra le posizioni dei veicoli
    for i in range(num_steps):
        if Pos_prec[i] == (None, None):
            error.append(np.nan)  # Uso di nan per indicare valori non validi
        else:
            diff = np.array(Pos_vehicle[i]) - np.array(Pos_prec[i])
            error.append(np.linalg.norm(diff))


    # Plot dell'errore tra le posizioni dei veicoli
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
    k1 = 0.25
    k2 = 0.25
    r = 1
    h = 0.2
    T = 0.1  # Passo di campionamento
    N = 1  # Numero di veicoli
    T_max = 20  # Tempo massimo
    num_steps = int(T_max/T)  # Numero di passi
    vehicles = []
    velocity = 5

    x_positions = [[] for _ in range(N)]
    y_positions = [[] for _ in range(N)]
    times = [k * T for k in range(num_steps)]

    acceleration = []
    omega = []
    
    traiettoria=input("Inserire 1 per traiettoria circolare, 2 per traiettoria curvilinea: ")
    for i in range(num_steps):
        acceleration.append(0)
    if(traiettoria=="1"):
        for i in range(num_steps):
            if (i < 7/T):
                omega.append(0)
            else:
                omega.append(0.5)
    elif(traiettoria=="2"):
        for i in range(num_steps):
            if(i< 3/T):
                omega.append(0)
            elif (i>=3/T and i<8/T):
                omega.append(0.5)
            elif (i>=8/T and i<12/T):
                omega.append(-0.5)
            elif (i>=12/T and i<15/T):
                omega.append(0)
            else:
                omega.append(0.5)
    else:
        print("Scelta non valida")

    contr = input("Inserire 1 per Controller_standard, 2 per Controller_extend: ")
    if contr == "1":
        for i in range(N):
            if i == 0:
                vehicles.append(Vehicle(True, Controller_standard(r, h, k1, k2)))
            else:
                vehicles.append(Vehicle(False, Controller_standard(r, h, k1, k2)))
            vehicles[i].states.append(VehicleState(-i, i, 0, velocity))
            if i != 0:
                previous_vehicle = vehicles[i - 1]
            else:
                previous_vehicle = None
            vehicles[i].controller.update_state_init(0, previous_vehicle, vehicles[i], acceleration[0], omega[0])
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

    else:
        print("Scelta non valida")

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

    # Plot traiettorie
    plt.subplot(3,1,1)
    for i in range(N):
        plt.plot(x_positions[i], y_positions[i], label=f'Vehicle {i + 1} Trajectory')
    plt.title('Vehicle Trajectories')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.legend()
    plt.grid(True)

    # Plot posizione lungo x vs tempo
    plt.subplot(3,1,2)
    for i in range(N):
        plt.plot(times, x_positions[i], label=f'Vehicle {i + 1} X Position')
    plt.title('X Position vs Time')
    plt.xlabel('Time (s)')
    plt.ylabel('X Position')
    plt.legend()
    plt.grid(True)

    # Plot posizione lungo y vs tempo
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

    # 3D Plot
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

    # Animazione delle traiettorie
    fig, ax = plt.subplots()
    lines = [ax.plot([], [], label=f'Vehicle {i + 1}')[0] for i in range(N)]
    x_min, x_max = min(map(min, x_positions)), max(map(max, x_positions))
    y_min, y_max = min(map(min, y_positions)), max(map(max, y_positions))
    margin_x = (x_max - x_min) * 0.1
    margin_y = (y_max - y_min) * 0.1
    ax.set_xlim(x_min - margin_x, x_max + margin_x)
    ax.set_ylim(y_min - margin_y, y_max + margin_y)
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_title('Vehicle Trajectories')
    ax.legend()
    ax.grid(True)

    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    def update(frame):
        for i, line in enumerate(lines):
            line.set_data(x_positions[i][:frame], y_positions[i][:frame])
        return lines

    ani = FuncAnimation(fig, update, frames=num_steps, init_func=init, blit=False)
    ani.save('vehicle_trajectories.gif', writer=PillowWriter(fps=30))

    plt.show()

    fig, ax = plt.subplots()
    rects = [patches.Rectangle((x_positions[i][0], y_positions[i][0]), 2, 0.5, angle=0, color=np.random.rand(3, ),
                               label=f'Vehicle {i + 1}') for i in range(N)]
    for rect in rects:
        ax.add_patch(rect)

    x_min, x_max = min(map(min, x_positions)), max(map(max, x_positions))
    y_min, y_max = min(map(min, y_positions)), max(map(max, y_positions))
    margin_x = (x_max - x_min) * 0.1
    margin_y = (y_max - y_min) * 0.1
    ax.set_xlim(x_min - margin_x, x_max + margin_x)
    ax.set_ylim(y_min - margin_y, y_max + margin_y)
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_title('Vehicle Trajectories Rectangles')
    ax.legend()
    ax.grid(True)

    def init():
        for rect in rects:
            rect.set_xy((0, 0))
        return rects

    def update(frame):
        for i, rect in enumerate(rects):
            rect.set_xy((x_positions[i][frame], y_positions[i][frame]))
        return rects

    ani = FuncAnimation(fig, update, frames=num_steps, init_func=init, blit=False)
    ani.save('vehicle_trajectories_rectangles.gif', writer=PillowWriter(fps=30))

    plt.show()


    # Calcolo dell'errore tra i veicoli
    v = []
    Error_total = []
    for i in range(1, N):
        for j in range(num_steps):
            v.append(vehicles[i - 1].states[j].velocity)
        e = Error(vehicles[i], vehicles[i - 1], T_max, T, h, v, r, i)
        Error_total.append(e)

    #Calcolo dell'errore medio di ogni veicolo
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

        print(f"Errore medio tra i veicoli {j} e {j-1} : {Error_medio[j - 1]}")

if __name__ == "__main__":
    main()




