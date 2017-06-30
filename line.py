from vector import Vector

class Line(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 2

        if not normal_vector:
            all_zeros = [0]*self.dimension
            normal_vector = all_zeros
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

            initial_index = Line.first_nonzero_index(n)
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
            initial_index = Line.first_nonzero_index(n)
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
        raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)

    def is_parallel_to(self, line):
        # check if they have the same normal vector. two lines are parallel if their normal vectors are parallel
        return self.normal_vector.is_parallel_to(line.normal_vector) 

    def __eq__(self, line):
        if self.normal_vector.is_zero():
            if not line.normal_vector.is_zero():
                return False
            else:
                diff = self.constant_term - line.constant_term
                return MyFloat(diff).is_near_zero()
        elif line.normal_vector.is_zero():
            return False

        if not self.is_parallel_to(line):
            return False
        # call parallel check. vector connecting each line is orthogonal to the normal vector
        conVec = line.basepoint.subtract(self.basepoint)
        return conVec.is_orthogonal_to(self.normal_vector)

    def intersection(self, line):
        # check if equal
        if self.__eq__(line):
            return self

        # check if parallel
        if self.is_parallel_to(line):
            return None

        # calculate intersection
        A, B = self.normal_vector.coordinates
        C, D = line.normal_vector.coordinates
        k1 = self.constant_term
        k2 = line.constant_term

        x_numerator = D*k1 - B*k2
        y_numerator = -C*k1 + A*k2
        one_over_denom = 1/(A*D - B*C)

        return Vector([x_numerator, y_numerator]).scale(one_over_denom)



class MyFloat(float):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps



# # Quiz #1 Intersection
# ell1 = Line(normal_vector=Vector([4.046,2.836]), constant_term=1.21)
# ell2 = Line(normal_vector=Vector([10.115,7.09]), constant_term=3.025)
# print('intersection 1:', ell1.intersection(ell2))

# ell1 = Line(normal_vector=Vector([7.204,3.182]), constant_term=8.68)
# ell2 = Line(normal_vector=Vector([8.172,4.114]), constant_term=9.883)
# print('intersection 2:', ell1.intersection(ell2))

# ell1 = Line(normal_vector=Vector([1.182,5.562]), constant_term=6.744)
# ell2 = Line(normal_vector=Vector([1.773,8.343]), constant_term=9.525)
# print('intersection 3:', ell1.intersection(ell2))




