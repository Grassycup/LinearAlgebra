from copy import deepcopy
from vector import Vector
from plane import Plane
from parametrization import Parametrization

class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def swap_rows(self, row1, row2):
        temp_plane = self[row1]
        self[row1] = self[row2]
        self[row2] = temp_plane


    def multiply_coefficient_and_row(self, coefficient, row):
        self[row].normal_vector = self[row].normal_vector.scale(coefficient)
        self[row].constant_term*=coefficient
        self[row].set_basepoint()

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        row_to_add_plane = self[row_to_add]
        scaled_vector = row_to_add_plane.normal_vector.scale(coefficient)
        scaled_coeff = row_to_add_plane.constant_term*coefficient
        row_to_be_added_plane = self[row_to_be_added_to]
        row_to_be_added_plane.normal_vector = row_to_be_added_plane.normal_vector.add(scaled_vector)
        row_to_be_added_plane.constant_term += scaled_coeff
        row_to_be_added_plane.set_basepoint()

    def compute_triangular_form(self):
        system = deepcopy(self)
        num_eq = len(system)
        num_var = system.dimension
        col = 0
        for i in range(num_eq):
            while col < num_var:
                coeff = system[i].normal_vector[col]
                if(MyFloat(coeff).is_near_zero):
                    non_zero_coeff_col = system.contains_nonzero_coeff_in_nth_var(i, col)
                    if non_zero_coeff_col is not -1: #there's a row under row i with nonzero coeff for var j
                        system.swap_rows(non_zero_coeff_col, i) #swap that row with row i
                    else:
                        col+=1
                        continue
                system.clear_terms_blow(i, col)# clear all terms with var j below row i
                col+=1
                break;
        return system

    def compute_rref(self):
        tf = self.compute_triangular_form()
        try:
            # find the first non zero term for each row
            first_non_zero_term_array = tf.indices_of_first_nonzero_terms_in_each_row()
            # process row by row
            for i,j in enumerate(first_non_zero_term_array):
                # ignore if there's no non zero term
                if j == -1:
                    continue
                # scale term to make it 1
                tf.scale_row_to_make_coefficient_equal_one(i, j)
                # clear terms above 
                tf.clear_terms_above(i, j)
            return tf
        except Exception as e:
            if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                return tf
            else:
                raise e

    def find_solutions(self):
        try:
            rref = self.compute_rref()
            coefficient_list = rref.indices_of_first_nonzero_terms_in_each_row()

            for p in rref.planes:
                if not p.contains_nonzero_index(p.normal_vector) and not MyFloat(p.constant_term).is_near_zero():
                    return rref.NO_SOLUTIONS_MSG

            # check for infinitely many solution
            pivot_num = [ x for x in coefficient_list if x != -1 ]
            if len(pivot_num) != rref.dimension:
                return rref.paramatrize_infinite_solutions()

            return [ rref.planes[x].constant_term for x in range(rref.dimension) ]
        except Exception as e:
            if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                return rref.NO_SOLUTIONS_MSG
            else:
                raise e

    def paramatrize_infinite_solutions(self):
        dimension = self.dimension
        pivot_indicies = self.indices_of_first_nonzero_terms_in_each_row()
        free_indicies = (set(range(dimension)) - set(pivot_indicies))
        
        direction_vectors = []
        basepoint = [0]*dimension
        for index,plane in enumerate(self.planes):
            if index in free_indicies:
                free_vector = [0]*dimension
                for i,p in enumerate(self.planes):
                    free_vector[index] = 1
                    pivot_val = pivot_indicies[i]
                    if pivot_val < 0:
                        break
                    free_vector[pivot_val] = -p.normal_vector[index]
                direction_vectors.append(Vector(free_vector))
            
            pivot_var = pivot_indicies[index]
            if pivot_var < 0:
                continue
            basepoint[pivot_var] = plane.constant_term
        param = Parametrization(Vector(basepoint), direction_vectors)
        return param


    def contains_nonzero_coeff_in_nth_var(self, row, col):
        result = -1
        for i in range(row, len(self.planes)):
            # check if current row's nth(col) variable is not zero
            if not MyFloat(self[i].normal_vector[col]).is_near_zero():
                return i
        return result

    def clear_terms_blow(self, row, col):
        top_row = self[row]
        for i in range(row + 1, len(self.planes)):
            # check if current row's nth(col) variable is not zero
            if not MyFloat(self[i].normal_vector[col]).is_near_zero():
                # clear out that term
                self.add_multiple_times_row_to_row(-1*self[i].normal_vector[col] / top_row.normal_vector[col], row, i)

    def clear_terms_above(self, row, col):
        top_row = self[row]
        for i in range(row):
            # check if current row's nth(col) variable is not zero
            if not MyFloat(self[i].normal_vector[col]).is_near_zero():
                # clear out that term
                self.add_multiple_times_row_to_row(-1*self[i].normal_vector[col], row, i)

    def scale_row_to_make_coefficient_equal_one(self, row, col):
        self.multiply_coefficient_and_row(1/self[row].normal_vector[col], row)

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices


    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyFloat(float):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

# gaussian elimination row operations
p0 = Plane(normal_vector=Vector([1,1,1]), constant_term=1)
p1 = Plane(normal_vector=Vector([0,1,0]), constant_term=2)
p2 = Plane(normal_vector=Vector([1,1,-1]), constant_term=3)
p3 = Plane(normal_vector=Vector([1,0,-2]), constant_term=2)

s = LinearSystem([p0,p1,p2,p3])

s.swap_rows(0,1)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print('test case 1 failed')

s.swap_rows(1,3)
if not (s[0] == p1 and s[1] == p3 and s[2] == p2 and s[3] == p0):
    print('test case 2 failed')

s.swap_rows(3,1)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print('test case 3 failed')

s.multiply_coefficient_and_row(1,0)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print('test case 4 failed')

s.multiply_coefficient_and_row(-1,2)
if not (s[0] == p1 and
        s[1] == p0 and
        s[2] == Plane(normal_vector=Vector([-1,-1,1]), constant_term=-3) and
        s[3] == p3):
    print('test case 5 failed')

s.multiply_coefficient_and_row(10,1)
if not (s[0] == p1 and
        s[1] == Plane(normal_vector=Vector([10,10,10]), constant_term=10) and
        s[2] == Plane(normal_vector=Vector([-1,-1,1]), constant_term=-3) and
        s[3] == p3):
    print('test case 6 failed')

s.add_multiple_times_row_to_row(0,0,1)
if not (s[0] == p1 and
        s[1] == Plane(normal_vector=Vector([10,10,10]), constant_term=10) and
        s[2] == Plane(normal_vector=Vector([-1,-1,1]), constant_term=-3) and
        s[3] == p3):
    print('test case 7 failed')


s.add_multiple_times_row_to_row(1,0,1)
if not (s[0] == p1 and
        s[1] == Plane(normal_vector=Vector([10,11,10]), constant_term=12) and
        s[2] == Plane(normal_vector=Vector([-1,-1,1]), constant_term=-3) and
        s[3] == p3):
    print('test case 8 failed')

s.add_multiple_times_row_to_row(-1,1,0)
if not (s[0] == Plane(normal_vector=Vector([-10,-10,-10]), constant_term=-10) and
        s[1] == Plane(normal_vector=Vector([10,11,10]), constant_term=12) and
        s[2] == Plane(normal_vector=Vector([-1,-1,1]), constant_term=-3) and
        s[3] == p3):
    print('test case 9 failed')


# Triangular form
p1 = Plane(normal_vector=Vector([1,1,1]), constant_term=1)
p2 = Plane(normal_vector=Vector([0,1,1]), constant_term=2)
s = LinearSystem([p1,p2])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == p2):
    print('test case 1 failed')

p1 = Plane(normal_vector=Vector([1,1,1]), constant_term=1)
p2 = Plane(normal_vector=Vector([1,1,1]), constant_term=2)
s = LinearSystem([p1,p2])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == Plane(constant_term=1)):
    print('test case 2 failed')

p1 = Plane(normal_vector=Vector([1,1,1]), constant_term=1)
p2 = Plane(normal_vector=Vector([0,1,0]), constant_term=2)
p3 = Plane(normal_vector=Vector([1,1,-1]), constant_term=3)
p4 = Plane(normal_vector=Vector([1,0,-2]), constant_term=2)
s = LinearSystem([p1,p2,p3,p4])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == p2 and
        t[2] == Plane(normal_vector=Vector([0,0,-2]), constant_term=2) and
        t[3] == Plane()):
    print('test case 3 failed')

p1 = Plane(normal_vector=Vector([0,1,1]), constant_term=1)
p2 = Plane(normal_vector=Vector([1,-1,1]), constant_term=2)
p3 = Plane(normal_vector=Vector([1,2,-5]), constant_term=3)
s = LinearSystem([p1,p2,p3])
t = s.compute_triangular_form()
if not (t[0] == Plane(normal_vector=Vector([1,-1,1]), constant_term=2) and
        t[1] == Plane(normal_vector=Vector([0,1,1]), constant_term=1) and
        t[2] == Plane(normal_vector=Vector([0,0,-9]), constant_term=-2)):
    print('test case 4 failed')


#RREF
p1 = Plane(normal_vector=Vector([1,1,1]), constant_term=1)
p2 = Plane(normal_vector=Vector([0,1,1]), constant_term=2)
s = LinearSystem([p1,p2])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector([1,0,0]), constant_term=-1) and
        r[1] == p2):
    print('test case 1 failed')

p1 = Plane(normal_vector=Vector([1,1,1]), constant_term=1)
p2 = Plane(normal_vector=Vector([1,1,1]), constant_term=2)
s = LinearSystem([p1,p2])
r = s.compute_rref()
if not (r[0] == p1 and
        r[1] == Plane(constant_term=1)):
    print('test case 2 failed')

p1 = Plane(normal_vector=Vector([1,1,1]), constant_term=1)
p2 = Plane(normal_vector=Vector([0,1,0]), constant_term=2)
p3 = Plane(normal_vector=Vector([1,1,-1]), constant_term=3)
p4 = Plane(normal_vector=Vector([1,0,-2]), constant_term=2)
s = LinearSystem([p1,p2,p3,p4])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector([1,0,0]), constant_term=0) and
        r[1] == p2 and
        r[2] == Plane(normal_vector=Vector([0,0,-2]), constant_term=2) and
        r[3] == Plane()):
    print('test case 3 failed')

p1 = Plane(normal_vector=Vector([0,1,1]), constant_term=1)
p2 = Plane(normal_vector=Vector([1,-1,1]), constant_term=2)
p3 = Plane(normal_vector=Vector([1,2,-5]), constant_term=3)
s = LinearSystem([p1,p2,p3])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector([1,0,0]), constant_term=23/9) and
        r[1] == Plane(normal_vector=Vector([0,1,0]), constant_term=7/9) and
        r[2] == Plane(normal_vector=Vector([0,0,1]), constant_term=2/9)):
    print('test case 4 failed')


## Find solutions
# p1 = Plane(normal_vector=Vector([5.862, 1.178, -10.366]), constant_term=-8.15)
# p2 = Plane(normal_vector=Vector([-2.931, -0.589, 5.183]), constant_term=-4.075)
# system = LinearSystem([p1, p2])
# print(system.find_solutions())

# p1 = Plane(Vector([8.631, 5.112, -1.816]), -5.113)
# p2 = Plane(Vector([4.315, 11.132, -5.27]), -6.775)
# p3 = Plane(Vector([-2.158, 3.01, -1.727]), -0.831)
# system = LinearSystem([p1, p2, p3])
# print(system.find_solutions())

# p1 = Plane(Vector([5.262, 2.739, -9.878]), -3.441)
# p2 = Plane(Vector([5.111, 6.358, 7.638]), -2.152)
# p3 = Plane(Vector([2.016, -9.924, -1.367]), -9.278)
# p4 = Plane(Vector([2.167, -13.543, -18.883]), -10.567)
# system = LinearSystem([p1, p2, p3, p4])
# print(system.find_solutions())


# parametrization
p1 = Plane(normal_vector=Vector([0.786, 0.786, 0.588]), constant_term=-0.714)
p2 = Plane(normal_vector=Vector([-0.131, -0.131, 0.244]), constant_term=0.319)
system = LinearSystem([p1, p2])
print(system.find_solutions())

p1 = Plane(Vector([8.631, 5.112, -1.816]), -5.113)
p2 = Plane(Vector([4.315, 11.132, -5.27]), -6.775)
p3 = Plane(Vector([-2.158, 3.01, -1.727]), -0.831)
system = LinearSystem([p1, p2, p3])
print(system.find_solutions())

p1 = Plane(Vector([0.935, 1.76, -9.365]), -9.955)
p2 = Plane(Vector([0.187, 0.352, -1.873]), -1.991)
p3 = Plane(Vector([0.374, 0.704, -3.746]), -3.982)
p4 = Plane(Vector([-0.561, -1.056, 5.619]), 5.973)
system = LinearSystem([p1, p2, p3, p4])
print(system.find_solutions())
