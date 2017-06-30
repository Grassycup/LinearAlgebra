from vector import Vector

class Plane(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 3

        if not normal_vector:
            all_zeros = [0]*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = 0
        self.constant_term = constant_term

        self.set_basepoint()


    def set_basepoint(self):
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = [0]*self.dimension

            initial_index = self.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e


    def __str__(self):

        num_float_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_float_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '

            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector

        try:
            initial_index = self.first_nonzero_index(n)
            terms = [write_coefficient(n[i], is_initial_term=(i==initial_index)) + 'x_{}'.format(i+1)
                     for i in range(self.dimension) if round(n[i], num_float_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_float_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output


    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyFloat(item).is_near_zero():
                return k
        raise Exception(Plane.NO_NONZERO_ELTS_FOUND_MSG)
        
    @staticmethod
    def contains_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyFloat(item).is_near_zero():
                return True
        return False
        
    def is_parallel_to(self, plane):
        # check if they have the same normal vector. two planes are parallel if their normal vectors are parallel
        return self.normal_vector.is_parallel_to(plane.normal_vector) 

    def __eq__(self, plane):
        if self.normal_vector.is_zero():
            if not plane.normal_vector.is_zero():
                return False
            else:
                diff = self.constant_term - plane.constant_term
                return MyFloat(diff).is_near_zero()
        elif plane.normal_vector.is_zero():
            return False

        if not self.is_parallel_to(plane):
            return False
        # call parallel check. vector connecting each plane is orthogonal to the normal vector
        conVec = plane.basepoint.subtract(self.basepoint)
        return conVec.is_orthogonal_to(self.normal_vector)

    def intersection(self, plane):
        # check if equal
        if self.__eq__(plane):
            return self

        # check if parallel
        if self.is_parallel_to(plane):
            return None

        # calculate intersection
        A, B = self.normal_vector.coordinates
        C, D = plane.normal_vector.coordinates
        k1 = self.constant_term
        k2 = plane.constant_term

        x_numerator = D*k1 - B*k2
        y_numerator = -C*k1 + A*k2
        one_over_denom = 1/(A*D - B*C)

        return Vector([x_numerator, y_numerator]).scale(one_over_denom)



class MyFloat(float):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps



# plane1 = Plane(Vector([-0.412, 3.806, 0.728]), -3.46)
# plane2 = Plane(Vector([1.03, -9.515, -1.82]), 8.65)
# print('1 is parallel: {}'.format(plane1.is_parallel_to(plane2)))
# print('1 is equal: {}'.format(plane1 == plane2))

# plane3 = Plane(Vector([2.611, 5.528, 0.283]), 4.6)
# plane4 = Plane(Vector([7.715, 8.306, 5.342]), 3.76)
# print('2 is parallel: {}'.format(plane3.is_parallel_to(plane4)))
# print('2 is equal: {}'.format(plane3 == plane4))

# plane5 = Plane(Vector([-7.926, 8.625, -7.212]), -7.952)
# plane6 = Plane(Vector([-2.642, 2.875, -2.404]), -2.443)
# print('3 is parallel: {}'.format(plane5.is_parallel_to(plane6)))
# print('3 is equal: {}'.format(plane5 == plane6))




