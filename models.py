

class MaterialPoint:
    '''
    Иниализация класса материальной точки
    '''
    def __init__(self, i, coord_x, coord_y, velocity_x, velocity_y, x_0, y_0, t):

        self.i = i
        self.coord_x = coord_x
        self.coord_y = coord_y
        self.velocity_x = velocity_x
        self.velocity_y = velocity_y
        self.x_0 = x_0
        self.y_0 = y_0
        self.t = t


class MaterialBody:

    '''
    Инициализция тела
    '''
    def __init__(self, material_points):
        self.material_points = material_points


class PointTrajectory:
    '''
    Класс траектории точки'
    '''
    def __init__(self, material_point, x, y):
        self.material_point = material_point
        self.x = x
        self.y = y


class BodyTrajectory:
    '''
    Класс для траектории тела
    '''
    def __init__(self, point_trajectories, material_body):
        self.point_trajectories = point_trajectories
        self.material_body = material_body


class SpacePoint:
    '''
    Класс точек в пространстве для распределения скоростей
    '''

    def __init__(self, i, coord_x, coord_y, velocity_x, velocity_y, t):
        self.i = i
        self.coord_x = coord_x
        self.coord_y = coord_y
        self.velocity_x = velocity_x
        self.velocity_y = velocity_y
        self.t = t


class SpaceGrid:
    '''
    Класс для хранения данных класса SpacePoint
    '''

    def __init__(self, space_points):
        self.space_points = space_points
