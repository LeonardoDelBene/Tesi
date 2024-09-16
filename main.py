from Controller_extend import Controller_extend
from vehicle import *
from Controller_standard import *
from Vehicle_State import *
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.transforms as transforms
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
    k1 = 2.5
    k2 = 2.5
    r = 1
    h = 0.2
    T = 0.01  # Passo di campionamento
    N = 2 # Numero di veicoli
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
            vehicles[i].controller.update_state_init(0, previous_vehicle, vehicles[i], acceleration[0], omega[0],T)
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
            vehicles[i].controller.update_state_init(0, previous_vehicle, vehicles[i], acceleration[0], omega[0],T)
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
    plt.figure(figsize=(10, 6))

    # Plot traiettorie
    plt.subplot(3,1,1)
    #plt.figure(figsize=(15, 10))
    for i in range(N):
        plt.plot(x_positions[i], y_positions[i], label=f'Vehicle {i + 1} Trajectory')
    plt.title('Vehicle Trajectories')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.legend()
    plt.grid(True)
    #plt.show()


    # Plot posizione lungo x vs tempo
    plt.subplot(3,1,2)
    #plt.figure(figsize=(15, 10))
    for i in range(N):
        plt.plot(times, x_positions[i], label=f'Vehicle {i + 1} X Position')
    plt.title('X Position vs Time')
    plt.xlabel('Time (s)')
    plt.ylabel('X Position')
    plt.legend()
    plt.grid(True)
    #plt.show()


    # Plot posizione lungo y vs tempo
    plt.subplot(3, 1, 3)
    #plt.figure(figsize=(15, 10))
    for i in range(N):
        plt.plot(times, y_positions[i], label=f'Vehicle {i + 1} Y Position')
    plt.title('Y Position vs Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Y Position')
    plt.legend()
    plt.grid(True)
    #plt.show()


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

    # Calcolo dell'errore tra i veicoli
    v = []
    Error_total = []
    for i in range(1, N):
        for j in range(num_steps):
            v.append(vehicles[i - 1].states[j].velocity)
        e = Error(vehicles[i], vehicles[i - 1], T_max, T, h, v, r, i)
        Error_total.append(e)

    # Calcolo dell'errore medio di ogni veicolo
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

        print(f"Errore medio tra i veicoli {j} e {j - 1} : {Error_medio[j - 1]}")

    # Animazione delle traiettorie
    fig, ax = plt.subplots()

    # Creare le linee per la traiettoria e i pallini per le posizioni attuali
    lines = [ax.plot([], [], label=f'Vehicle {i + 1}')[0] for i in range(N)]
    scatters = [ax.scatter([], [], s=20, color=lines[i].get_color()) for i in
                range(N)]  # Pallini alle posizioni correnti

    # Impostazione dei limiti degli assi
    x_min, x_max = min(map(min, x_positions)), max(map(max, x_positions))
    y_min, y_max = min(map(min, y_positions)), max(map(max, y_positions))
    margin_x = (x_max - x_min) * 0.1
    margin_y = (y_max - y_min) * 0.1
    ax.set_xlim(x_min - margin_x, x_max + margin_x)
    ax.set_ylim(y_min - margin_y, y_max + margin_y)

    # Label e titoli
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_title('Vehicle Trajectories')
    ax.legend()
    ax.grid(True)

    # Funzione di inizializzazione
    def init():
        for line, scatter in zip(lines, scatters):
            line.set_data([], [])  # Linea vuota
            scatter.set_offsets(np.array([[0, 0]]))  # Pallino iniziale fittizio
        return lines + scatters

    # Funzione di aggiornamento per ogni frame
    def update(frame):
        for i, (line, scatter) in enumerate(zip(lines, scatters)):
            # Aggiornare la linea con la traiettoria fino al frame attuale
            line.set_data(x_positions[i][:frame], y_positions[i][:frame])

            # Aggiornare il pallino alla posizione corrente
            scatter.set_offsets(np.array([[x_positions[i][frame - 1], y_positions[i][frame - 1]]]))

        return lines + scatters

    # Creare l'animazione
    ani = FuncAnimation(fig, update, frames=num_steps, init_func=init, blit=False)

    # Salvare l'animazione come GIF
    ani.save('vehicle_trajectories.gif', writer=PillowWriter(fps=30))
    plt.show()

    # Animazione delle traiettorie con rettangoli
    thetas = [[state.theta for state in vehicle.states] for vehicle in vehicles]  # Theta in radianti
    angles = [np.degrees(theta) for theta in thetas]  # Theta in gradi
    fig, ax = plt.subplots()
    # Creare una lista di rettangoli iniziali
    rects = [
        patches.Rectangle(
            (x_positions[i][0], y_positions[i][0]), 1, 0.5, color=np.random.rand(3, ), label=f'Vehicle {i + 1}'
        )
        for i in range(N)
    ]

    # Aggiungere i rettangoli all'asse
    for rect in rects:
        ax.add_patch(rect)

    # Calcolare i limiti degli assi
    x_min, x_max = min(map(min, x_positions)), max(map(max, x_positions))
    y_min, y_max = min(map(min, y_positions)), max(map(max, y_positions))

    # Calcola un margine attorno ai dati
    margin_x = (x_max - x_min) * 0.1
    margin_y = (y_max - y_min) * 0.1

    # Impostare i limiti degli assi con margini
    ax.set_xlim(x_min - margin_x, x_max + margin_x)
    ax.set_ylim(y_min - margin_y, y_max + margin_y)

    # Impostare un aspect ratio uguale per mantenere le proporzioni corrette
    ax.set_aspect('equal')

    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_title('Vehicle Trajectories Rectangles')
    ax.legend()
    ax.grid(True)

    # Funzione di inizializzazione
    def init():
        for rect in rects:
            rect.set_xy((0, 0))
        return rects

    # Funzione di aggiornamento per ogni frame
    def update(frame):
        for i, rect in enumerate(rects):
            # Aggiornare la posizione del rettangolo
            rect.set_xy((x_positions[i][frame], y_positions[i][frame]))

            # Creare la trasformazione con rotazione attorno al centro del rettangolo
            trans = transforms.Affine2D().rotate_deg_around(
                x_positions[i][frame], y_positions[i][frame], angles[i][frame]  # Angolo in gradi
            ) + ax.transData

            # Applicare la trasformazione al rettangolo
            rect.set_transform(trans)

        return rects

    # Creare l'animazione
    ani = FuncAnimation(fig, update, frames=num_steps, init_func=init, blit=False)

    # Salvare l'animazione come GIF
    ani.save('vehicle_trajectories_rectangles.gif', writer=PillowWriter(fps=30))

    plt.show()

if __name__ == "__main__":
    main()




