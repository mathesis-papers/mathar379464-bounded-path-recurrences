import json
from fractions import Fraction
from math import comb

import sympy as sp

variable, index = sp.symbols('x n')
parameter = sp.Symbol('z')


def trim(operator):
    result = [sp.cancel(coefficient) for coefficient in operator]
    while result and result[-1] == 0:
        result.pop()
    return result


def add(left, right):
    result = [sp.Integer(0)] * max(len(left), len(right))
    for position, coefficient in enumerate(left):
        result[position] += coefficient
    for position, coefficient in enumerate(right):
        result[position] += coefficient
    return trim(result)


def scale(scalar, operator):
    return trim([scalar * coefficient for coefficient in operator])


def compose(left, right):
    result = [sp.Integer(0)] * (len(left) + len(right))
    for left_order, left_coefficient in enumerate(left):
        for right_order, right_coefficient in enumerate(right):
            for differentiations in range(left_order + 1):
                result[left_order + right_order - differentiations] += (
                    left_coefficient * sp.binomial(left_order, differentiations)
                    * sp.diff(right_coefficient, variable, differentiations))
    return trim(result)


def polynomial(expression, operator, symbol):
    result = []
    for coefficient in sp.Poly(expression, symbol).all_coeffs():
        result = add(compose(result, operator), [coefficient])
    return result


def recurrence_operator(coefficients):
    result = []
    for shift, expression in enumerate(coefficients):
        result = add(result, scale(variable**shift, polynomial(
            expression.subs(index, index + shift), [0, variable], index)))
    return result


def coefficient_rows():
    return {
        1: [3*index*(3*index-4)*(3*index+4),
            -3*(2*index-1)*(18*index**2-18*index-7),
            3*(index-1)*(54*index**2-108*index+47),
            -54*(index-2)*(index-1)*(2*index-3),
            -229*(index-3)*(index-2)*(index-1)],
        2: [3*index*(3*index+4)*(index-3)*(3*index+8),
            3*(-45*index**4+54*index**3+192*index**2-27*index-20),
            9*(index-1)*(30*index**3-72*index**2-7*index+20),
            -3*(index-1)*(index-2)*(90*index**2-234*index+95),
            -(index-1)*(index-2)*(index-3)*(121*index+499),
            229*(index-1)*(index-2)*(index-3)*(index-4)],
        3: [3*(index+4)*(3*index+4)*(3*index+8),
            3*(-63*index**3-297*index**2-349*index-60),
            3*(189*index**3+270*index**2-229*index-140),
            15*(-63*index**3+117*index**2+44*index-64),
            689*index**3-5372*index**2+6946*index-1288,
            (index-4)*(201*index**2+2767*index-3011),
            -(index-5)*(579*index+257)*(index-4),
            229*(index-5)*(index-6)*(index-4)]}


def check_symbolic():
    modified = [-variable, variable*(1-variable)]
    operators = {
        1: (27*parameter**3-48*parameter,
            parameter**3+6*parameter**2+11*parameter+6),
        2: (27*parameter**3+108*parameter**2+96*parameter,
            parameter**3+10*parameter**2+31*parameter+30),
        3: (27*parameter**3+216*parameter**2+528*parameter+384,
            parameter**3+14*parameter**2+63*parameter+90)}
    bases = {
        1: (3*parameter*(3*parameter-1)*(3*parameter+1),
            8*(2*parameter+1)*(4*parameter+1)*(4*parameter+3)),
        2: (3*parameter*(3*parameter+1)*(3*parameter+2),
            8*(2*parameter+1)*(4*parameter+3)*(4*parameter+5)),
        3: (3*(parameter+1)*(3*parameter+1)*(3*parameter+2),
            8*(2*parameter+3)*(4*parameter+3)*(4*parameter+5))}
    rows = coefficient_rows()
    for height in (1, 2, 3):
        first, second = operators[height]
        numerator = sp.prod(4*index+height+offset for offset in range(4))
        denominator = (index+1)*sp.prod(
            3*index+height+offset for offset in (1, 2, 3))
        lower, upper = bases[height]
        assert sp.expand(lower.subs(parameter, index+1)*numerator
                         - upper.subs(parameter, index)*denominator) == 0
        assert sp.expand(64*lower.subs(parameter, parameter/4)-first) == 0
        assert sp.expand(64*upper.subs(parameter, parameter/4)-256*second) == 0
        operator = add(scale((1-variable)**4, polynomial(first, modified, parameter)),
                       scale(-256*variable**4, polynomial(second, modified, parameter)))
        if height == 1:
            operator = scale((1-variable)**-3, operator)
        if height == 2:
            operator = compose(
                [-3/(1-variable)**2, variable/(1-variable)**2], operator)
        assert add(recurrence_operator(rows[height]), scale(-1, operator)) == []
    geometric = 1/(1-variable)
    substitution = (variable/(1-variable))**4
    assert sp.cancel(variable*(1-variable)*sp.diff(geometric, variable)
                     - variable*geometric) == 0
    assert sp.cancel(variable*(1-variable)*geometric*sp.diff(substitution, variable)
                     - 4*geometric*substitution) == 0
    tau = sp.Rational(4, 3)
    critical = sp.Rational(27, 256)
    alpha_squared = 2*tau**4*critical/(12*critical*tau**2)
    assert alpha_squared == sp.Rational(8, 27)
    return ['three factorial-ratio polynomial identities',
            'three operator coefficient identities', 'Euler substitution identity',
            'square-root coefficient 8/27']


def paths(height, maximum):
    states = {0: 1}
    result = [1]
    for length in range(1, maximum+1):
        following = {}
        for level, multiplicity in states.items():
            for increment in (1, 0, -3):
                destination = level+increment
                if destination >= 1-height:
                    following[destination] = following.get(destination, 0)+multiplicity
        states = following
        result.append(states.get(0, 0))
    return result


def closed_form(height, length):
    total = Fraction(0)
    for down_steps in range(length//4+1):
        total += Fraction(height*comb(length, 4*down_steps)
                          * comb(4*down_steps+height, down_steps),
                          4*down_steps+height)
    assert total.denominator == 1
    return total.numerator


def all_floor_formula(floor_depth, length):
    exponent = floor_depth + 1
    total = Fraction(0)
    for down_steps in range(length//4+1):
        zero_free = Fraction(0)
        for correction in range(min(floor_depth//4, down_steps)+1):
            reduced = down_steps-correction
            zero_free += Fraction(
                (-1)**correction * comb(floor_depth-3*correction, correction)
                * exponent * comb(4*reduced+exponent, reduced),
                4*reduced+exponent)
        total += comb(length, 4*down_steps)*zero_free
    assert total.denominator == 1
    return total.numerator


def generalized_binomial(top, bottom):
    result = Fraction(1)
    for offset in range(bottom):
        result *= Fraction(top-offset, offset+1)
    return result


def power_coefficients(coefficients, exponent, maximum):
    result = [Fraction(1)]
    for degree in range(1, maximum+1):
        result.append(sum(
            ((exponent+1)*position-degree)*coefficients[position]
            * result[degree-position]
            for position in range(1, degree+1)) / degree)
    return result


def check_siblings():
    maximum = 120
    for depth in range(9):
        counts = paths(depth+1, maximum)
        assert all(counts[length] == all_floor_formula(depth, length)
                   for length in range(maximum+1))
    assert paths(4, maximum) == [closed_form(4, length)
                                 for length in range(maximum+1)]
    assert paths(5, 4)[4] == 5
    assert closed_form(5, 4) == 6
    coefficient_limit = 16
    coefficients = [Fraction(2*comb(4*degree+2, degree), 4*degree+2)
                    for degree in range(coefficient_limit+1)]
    assert coefficients[:3] == [1, 2, 9]
    comparisons = 0
    for multiplier in range(-5, 6):
        for degree in range(1, coefficient_limit+1):
            exponent = multiplier*degree
            actual = power_coefficients(coefficients, exponent, degree)[degree]
            expected = 2*multiplier*generalized_binomial(
                (2*multiplier+4)*degree-1, degree-1)
            assert actual == expected
            comparisons += 1
    assert power_coefficients(coefficients, -4, 2)[2] == 4
    assert closed_form(1, 4) == 2
    assert paths(1, 3)[3] == 1
    matrix_values = [1]
    for degree in range(31):
        assert matrix_values[0] == Fraction(comb(4*degree, degree), 3*degree+1)
        expected_row_sum = Fraction(2*comb(4*degree+2, degree), 4*degree+2)
        assert sum(matrix_values) == expected_row_sum
        following = [0]*(len(matrix_values)+1)
        for row, multiplicity in enumerate(matrix_values):
            following[0] += multiplicity*comb(row+2, 2)
            for column in range(1, row+2):
                following[column] += multiplicity*comb(row-column+3, 2)
        matrix_values = following
    return {
        'all_floors': {'floor_depth': [0, 8], 'length': [0, maximum],
                       'exact_comparisons': 9*(maximum+1)},
        'floor_minus_three_Raney_sum': [0, maximum],
        'unasserted_Raney_extension': {'floor_depth': 4, 'length': 4,
                                      'actual': 5, 'naive_sum': 6},
        'power_coefficients': {'multipliers': [-5, 5], 'degrees': [1, 16],
                               'exact_comparisons': comparisons},
        'A069271_inverse_power_correction': {'degree': 2, 'value': 4,
                                             'printed_value': -4},
        'A069271_PARI_factor_correction': {'index': 0, 'value': 1,
                                          'printed_value': 2},
        'A127902_offset_correction': {'index': 4, 'value': 1,
                                      'unshifted_formula': 2},
        'production_matrix': {'powers': [0, 30],
                              'upper_left_and_row_sum_comparisons': 62},
    }


def main():
    maximum = 240
    symbolic = check_symbolic()
    rows = coefficient_rows()
    ranges = {}
    expected_prefixes = {1: [1, 1, 1, 1, 2, 6, 16, 36],
                         2: [1, 1, 1, 1, 3, 11, 31, 71],
                         3: [1, 1, 1, 1, 4, 16, 46, 106]}
    for height in (1, 2, 3):
        counts = paths(height, maximum)
        assert counts[:8] == expected_prefixes[height]
        assert all(counts[length] == closed_form(height, length)
                   for length in range(maximum+1))
        functions = [sp.lambdify(index, expression, 'math') for expression in rows[height]]
        minimum = len(functions)-1
        for length in range(minimum, maximum+1):
            assert sum(function(length)*counts[length-shift]
                       for shift, function in enumerate(functions)) == 0
        ranges[str(height)] = {'path_formula': [0, maximum],
                               'recurrence': [minimum, maximum],
                               'initial_terms': counts[:12]}
    vertex_limit = 20
    trees = [0] * (vertex_limit+1)
    for vertices in range(1, vertex_limit+1):
        trees[vertices] = (1 if vertices == 1 else 0) + trees[vertices-1]
        trees[vertices] += sum(
            trees[first]*trees[second]*trees[third]*trees[vertices-1-first-second-third]
            for first in range(1, vertices)
            for second in range(1, vertices-first)
            for third in range(1, vertices-first-second))
    assert trees[1:] == paths(1, vertex_limit-1)
    assert trees[2] == 1
    siblings = check_siblings()
    print(json.dumps({'status': 'pass', 'arithmetic': 'exact integers and rationals',
                      'symbolic_identities': symbolic, 'finite_checks': ranges,
                      'sibling_checks': siblings,
                      'ordered_tree_vertices': [1, vertex_limit],
                      'labeled_tree_counterexample': {'vertices': 2,
                          'distinctly_labeled_rooted_trees': 2, 'A127902': trees[2]}},
                     indent=2))


if __name__ == '__main__':
    main()
