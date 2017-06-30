from math import sqrt, pi, acos

class Vector(object):
    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'
    NO_UNIQUE_PARALLEL_COMPONENT_MSG = 'Zero vector has no unique parallel components'

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple(coordinates)
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __len__(self):
        return len(self.coordinates)

    def __getitem__(self, i):
        return self.coordinates[i]

    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)

    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def __iter__(self):
        self.current = 0
        return self

    def __next__(self):
        if self.current >= len(self.coordinates):
            raise StopIteration
        else:
            current_value = self.coordinates[self.current]
            self.current += 1
            return current_value


    def add(self, v):
        coord = []
        for i in range(len(self.coordinates)):
            coord.append(self.coordinates[i] + v.coordinates[i])
        return Vector(coord)

    def subtract(self, v):
        coord = [x-y for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(coord)

    def scale(self, c):
        coord = [x*c for x in self.coordinates]
        return Vector(coord)

    def magnitude(self):
        coord = [x**2 for x in self.coordinates]
        return sqrt(sum(coord))

    def normalize(self):
        try:
            magnitude = self.magnitude()
            return self.scale(1./magnitude)
        except ZeroDivisionError:
            raise Exception(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)

    def inner_product(self, v):
        coord = [x*y for x,y, in zip(self.coordinates, v.coordinates)]
        return sum(coord)

    def angle_with(self, v, in_degrees=False):
        try:
            u1 = self.normalize()
            u2 = v.normalize()
            # acos has domain of [-1,1], so use rounding to prevent any errors due to precision
            angle_in_radians = acos(round(u1.inner_product(u2), 3))
            if in_degrees:
                degrees_per_radian = 180 / pi
                return angle_in_radians * degrees_per_radian
            else:
                return angle_in_radians
        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('Cannot compute an angle with a zero vector')
            else:
                raise e

    def is_parallel_to(self, v):
        return (self.is_zero() or 
                v.is_zero() or
                self.angle_with(v) == 0 or
                self.angle_with(v) == pi)

    def is_zero(self, tolerance=1e-10):
        return self.magnitude() < tolerance

    def is_orthogonal_to(self, v, tolerance=1e-10):
        return abs(self.inner_product(v)) < tolerance

    def component_parallel_to(self, basis):
        try:
            unit_basis = basis.normalize()
            return unit_basis.scale(self.inner_product(unit_basis))
        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Except(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            else:
                raise e

    def component_orthogonal_to(self, basis):
        try:
            return self.subtract(self.component_parallel_to(basis))
        except Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            else:
                raise e

    def cross_product(self, v):
        try:
            i = self.coordinates[1] * v.coordinates[2] - self.coordinates[2] * v.coordinates[1]
            j = -1 * (self.coordinates[0] * v.coordinates[2] - self.coordinates[2] * v.coordinates[0])
            k = self.coordinates[0] * v.coordinates[1] - self.coordinates[1] * v.coordinates[0]
            return Vector([i, j, k])
        except ValueError as e:
            raise e

    def area_parallelogram(self, v):
        return 0

    def area_triangle(self, v):
        return 0        




# # Quiz add, substract, scale
# v1 = Vector([8.218, -9.341])
# v2 = Vector([-1.129, 2.111])
# v3 = Vector([7.119, 8.215])
# v4 = Vector([-8.223, 0.878])
# v5 = Vector([1.671, -1.012, -0.318])
# print(v1.add(v2))
# print(v3.subtract(v4))
# print(v5.scale(7.41))

# # Quiz 2 magnitude, direction
# v6 = Vector([-0.221, 7.437])
# v7 = Vector([1.996, 3.108, -4.554])
# print(v6.magnitude())
# print(v7.normalize())


# # Quiz 3 inner product, angle
# v8 = Vector([7.887, 4.138])
# v9 = Vector([-8.802, 6.776])
# v10 = Vector([3.183, -7.627])
# v11 = Vector([-2.668, 5.319])
# print(v8.inner_product(v9))
# print(v10.angle_with(v11))

# # Quiz 4 parallel, orthogonal
# v12 = Vector([-7.579, -7.88])
# v13 = Vector([22.737, 23.64])
# print(str(v12.is_parallel_to(v13)))
# print(str(v12.is_orthogonal_to(v13)))

# # Quiz 5 orthogonal, parallel component
# v14 = Vector([3.039, 1.879])
# v15 = Vector([0.825, 2.036]) 
# v16 = Vector([-9.88, -3.264, -8.159])
# v17 = Vector([-2.155, -9.353, -9.473])
# v18 = Vector([3.009, -6.172, 3.692, -2.51])
# v19 = Vector([6.404, -9.144, 2.759, 8.718])
# print(v14.component_parallel_to(v15))
# print(v16.component_orthogonal_to(v17))
# print(v18.component_parallel_to(v19))
# print(v18.component_orthogonal_to(v19))

# Quiz 6 cross product
# v20 = Vector([8.462, 7.893, -8.187])
# v21 = Vector([6.984, -5.975, 4.778])
# print(v20.cross_product(v21))









