import models.models as mod
import math
import matplotlib.pyplot as plt
import numpy as np


def f_x(t, x):  # 1 скорость

    return -((math.exp(t) + math.exp(-t))/2) * x  # ch(t) * x


def f_y(t, y): # 2 скорость

    return ((math.exp(t) - math.exp(-t))/2) * y  # sh(t) * y


def create_material_body(x_c, y_c, h):  # Функция создания материального тела
    t = 0
    m = 0
    material_points = []  # Создание массива
    for i in range(int(2/h) + 1):
        for j in range(int(2/h) + 1):
            x = x_c - j * h
            y = y_c + i * h
            material_points.append(mod.MaterialPoint(m, x, y, f_x(t, x), f_y(t, y), x, y, t)) # Добавление точек в массив
            m += 1
    material_body = mod.MaterialBody(material_points)  # Формирование тела из точек
    return material_body


def runge_kutta(t, h, x_k, y_k):  # Численное интегрирование методом Рунге-Кутты

    c2, a21, b1, b2 = 1/2, 1/2, 0, 1  # Ввод количества стадий и коэффиценты из таблицы Бутчера (2.1)

    kx1 = f_x(t, x_k)
    kx2 = f_x((t + h * c2), (x_k + a21 * kx1 * h))
    x_t = x_k + h * (kx1 * b1 + kx2 * b2)

    ky1 = f_y(t, y_k)
    ky2 = f_y((t + h * c2), (y_k + a21 * kx1 * h))
    y_t = y_k + h * (ky1 * b1 + ky2 * b2)

    return x_t, y_t


def move_material_body(time, h, mb):  # Функция для движения материального тела
    point_trajectories = []  # Создание массива траекторий точек
    for i in range(len(mb.material_points)):
        t = 0
        x_0 = mb.material_points[i].x_0
        y_0 = mb.material_points[i].y_0
        x_t = [x_0]
        y_t = [y_0]
        for n in range(int(time / h) + 1):
            x_k = x_t[n]
            y_k = y_t[n]
            x_t.append(runge_kutta(t, h, x_k, y_k)[0])
            y_t.append(runge_kutta(t, h, x_k, y_k)[1])

            t += h

        point_trajectories.append(mod.PointTrajectory(mb.material_points[i], x_t, y_t))  # Добавление точек траектории в
                                                                                         # массив
    body_trajectory = mod.BodyTrajectory(point_trajectories, mb)
    return body_trajectory


def plot_trajectory(mb, tr): # Функция построения графика траектории

    for i in range(len(mb.material_points)):
        plt.plot(mb.material_points[i].coord_x, mb.material_points[i].coord_y, 'r.')
    for i in range(len(mb.material_points)):
        plt.plot(tr.point_trajectories[i].x, tr.point_trajectories[i].y, 'b', linewidth=0.5)
    for i in range(len(mb.material_points)):
        time = len(tr.point_trajectories[i].x) - 1
        plt.plot(tr.point_trajectories[i].x[time], tr.point_trajectories[i].y[time], 'g.')
    plt.axis('equal')
    plt.grid()
    # plt.show()
    plt.savefig('assets/plot_trajectory.png', format='png', dpi=300)


def move_through_space(time, h):  # Функция перемещения
    t = h
    m = 0
    a = np.linspace(-4, 4, 9)
    x_s, y_s = np.meshgrid(a, a)
    velocity_fields = []
    for n in range(int(time / h)):
        space_points = []
        for i in range(9):
            for j in range(9):
                x = x_s[i, j]
                y = y_s[i, j]
                space_points.append(mod.SpacePoint(m, x, y, f_x(t, x), f_y(t, y), t))
                m += 1
        velocity_fields.append(mod.SpaceGrid(space_points))
        t += h
    return velocity_fields


def plot_velocity_fields(vf):  # Функция построения графика скорости
    h = vf[0].space_points[0].t
    t = h
    for n in range(len(vf)):
        plt.figure(n)
        plt.suptitle('t = ' + str(t))
        m = 0
        coord_x = []
        coord_y = []
        v_x = []
        v_y = []
        for i in range(9):
            for j in range(9):
                coord_x.append(vf[n].space_points[m].coord_x)
                coord_y.append(vf[n].space_points[m].coord_y)
                v_x.append(vf[n].space_points[m].velocity_x)
                v_y.append(vf[n].space_points[m].velocity_y)
                m += 1
        plt.subplot(1, 2, 1)
        plt.quiver(coord_x, coord_y, v_x, v_y)
        for p in range(1, 4):
            for q in range(1, 4):
                x = np.linspace(-4, -0.1, 100)
                d = math.exp(t) / math.log(t + 1)
                c = q * (p ** d)
                y = c * ((-x) ** (-d))
                plt.subplot(1, 2, 2)
                plt.axis([-4, 4, -4, 4])
                plt.plot(x, y)
        t += h
        # plt.show()
        plt.savefig('assets/velocity_fields' + str(n) + '.png', format='png', dpi=300)


