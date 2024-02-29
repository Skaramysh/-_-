import managers.managers as man

if __name__ == "__main__":
    x_corner = -3
    y_corner = -3
    h_s = 1

    time = 1
    h = 0.01
    rect = man.create_material_body(x_corner, y_corner, h_s)
    tr = man.move_material_body(time, h, rect)

    man.plot_trajectory(rect, tr)


    vf = man.move_through_space(1, 0.1)
    man.plot_velocity_fields(vf)
