"""Independent closed-form benchmarks; run: python -m unittest discover -s tests -v."""

import math
import unittest

from beam_balance import Load, Support, polynomial_roots, solve_beam


class BeamBenchmarks(unittest.TestCase):
    def assert_close(self, actual, expected, tolerance=1e-8):
        self.assertAlmostEqual(actual, expected, delta=tolerance * max(1, abs(expected)))

    def simple(self, length, loads, ei=1600):
        return solve_beam(length, ei, [Support("Pinned", 0), Support("Roller", length)], loads)

    def test_reference_center_point(self):
        result = self.simple(6, [Load("Point", 3, 30)], 30000)
        self.assert_close(result.reactions[0].force, 15)
        self.assert_close(result.reactions[1].force, 15)
        self.assert_close(result.extrema("M")[1][0], 45)
        self.assert_close(result.extrema("M")[1][1], 3)
        self.assert_close(result.deflection(3)[0], -30*6**3/(48*30000))
        self.assert_close(result.deflection(3)[1], 0)

    def test_simple_uniform_load(self):
        result = self.simple(4, [Load("UDL", 0, 6, 4)])
        self.assert_close(result.reactions[0].force, 12)
        self.assert_close(result.shear_moment(2)[1], 12)
        self.assert_close(result.deflection(2)[0], -5*6*4**4/(384*1600))
        self.assert_close(result.deflection(0)[1], -6*4**3/(24*1600))
        self.assert_close(result.deflection(4)[1], 6*4**3/(24*1600))

    def test_partial_udl_exact_moment_extremum(self):
        result = self.simple(5, [Load("UDL", 1, 4, 3)])
        self.assert_close(result.reactions[0].force, 4.8)
        self.assert_close(result.reactions[1].force, 3.2)
        maximum, x = result.extrema("M")[1]
        self.assert_close(maximum, 7.68)
        self.assert_close(x, 2.2)

    def test_left_cantilever_point(self):
        result = solve_beam(4, 1600, [Support("Fixed", 0)], [Load("Point", 4, 10)])
        self.assert_close(result.reactions[0].force, 10)
        self.assert_close(result.reactions[0].moment, 40)
        self.assert_close(result.shear_moment(0)[1], -40)
        self.assert_close(result.deflection(4)[0], -10*4**3/(3*1600))
        self.assert_close(result.deflection(4)[1], -10*4**2/(2*1600))

    def test_right_cantilever_point(self):
        result = solve_beam(4, 1600, [Support("Fixed", 4)], [Load("Point", 0, 10)])
        self.assert_close(result.reactions[0].force, 10)
        self.assert_close(result.reactions[0].moment, -40)
        self.assert_close(result.shear_moment(4, "left")[1], -40)
        self.assert_close(result.deflection(0)[0], -10*4**3/(3*1600))
        self.assert_close(result.deflection(0)[1], 10*4**2/(2*1600))
        self.assert_close(result.deflection(4)[0], 0)
        self.assert_close(result.deflection(4)[1], 0)

    def test_cantilever_udl(self):
        for fixed_x, free_x, couple in ((0, 4, 48), (4, 0, -48)):
            with self.subTest(fixed_x=fixed_x):
                result = solve_beam(4, 1600, [Support("Fixed", fixed_x)], [Load("UDL", 0, 6, 4)])
                self.assert_close(result.reactions[0].force, 24)
                self.assert_close(result.reactions[0].moment, couple)
                self.assert_close(result.deflection(free_x)[0], -6*4**4/(8*1600))

    def test_fixed_fixed_uniform_load(self):
        result = solve_beam(4, 1600, [Support("Fixed", 0), Support("Fixed", 4)], [Load("UDL", 0, 6, 4)])
        self.assert_close(result.reactions[0].force, 12)
        self.assert_close(result.reactions[1].force, 12)
        self.assert_close(result.reactions[0].moment, 8)
        self.assert_close(result.reactions[1].moment, -8)
        self.assert_close(result.shear_moment(0)[1], -8)
        self.assert_close(result.shear_moment(4, "left")[1], -8)
        self.assert_close(result.shear_moment(2)[1], 4)
        self.assert_close(result.deflection(2)[0], -6*4**4/(384*1600))
        self.assert_close(result.extrema("y")[0][1], 2)

    def test_fixed_fixed_triangle_zero_nodal_dofs(self):
        result = solve_beam(8, 12000, [Support("Fixed", 0), Support("Fixed", 8)], [Load("Linear", 0, 0, 8, 3)])
        self.assert_close(result.reactions[0].force, 3.6)
        self.assert_close(result.reactions[1].force, 8.4)
        self.assert_close(result.reactions[0].moment, 6.4)
        self.assert_close(result.reactions[1].moment, -9.6)
        self.assert_close(result.deflection(8)[0], 0)
        self.assert_close(result.deflection(8)[1], 0)
        self.assertLess(result.deflection(4)[0], 0)

    def test_propped_cantilever(self):
        length, w = 4, 6
        result = solve_beam(length, 1600, [Support("Fixed", 0), Support("Roller", length)], [Load("UDL", 0, w, length)])
        self.assert_close(result.reactions[0].force, 5*w*length/8)
        self.assert_close(result.reactions[1].force, 3*w*length/8)
        self.assert_close(result.reactions[0].moment, w*length**2/8)
        maximum, x = result.extrema("M")[1]
        self.assert_close(maximum, 9*w*length**2/128)
        self.assert_close(x, 5*length/8)
        self.assert_close(result.deflection(length)[0], 0)

    def test_two_equal_continuous_spans(self):
        span, w = 4, 6
        result = solve_beam(2*span, 1600, [Support("Pinned", 0), Support("Roller", span), Support("Roller", 2*span)], [Load("UDL", 0, w, 2*span)])
        self.assert_close(result.reactions[0].force, 3*w*span/8)
        self.assert_close(result.reactions[1].force, 5*w*span/4)
        self.assert_close(result.reactions[2].force, 3*w*span/8)
        self.assert_close(result.shear_moment(span)[1], -w*span**2/8)
        self.assert_close(result.deflection(span)[0], 0)

    def test_triangular_load(self):
        length, w, ei = 6, 9, 2400
        result = self.simple(length, [Load("Linear", 0, 0, length, w)], ei)
        self.assert_close(result.reactions[0].force, w*length/6)
        self.assert_close(result.reactions[1].force, w*length/3)
        maximum, x = result.extrema("M")[1]
        self.assert_close(maximum, w*length**2/(9*math.sqrt(3)))
        self.assert_close(x, length/math.sqrt(3))
        self.assert_close(result.deflection(0)[1], -7*w*length**3/(360*ei))
        self.assert_close(result.deflection(length)[1], w*length**3/(45*ei))

    def test_trapezoid_superposition(self):
        linear = self.simple(6, [Load("Linear", 0, 3, 6, 9)])
        separate = self.simple(6, [Load("UDL", 0, 3, 6), Load("Linear", 0, 0, 6, 6)])
        for x in (0, 1, 3, 5, 6):
            for actual, expected in zip(linear.shear_moment(x), separate.shear_moment(x)):
                self.assert_close(actual, expected)
            self.assert_close(linear.deflection(x)[0], separate.deflection(x)[0])

    def test_overlapping_uniform_loads(self):
        result = self.simple(10, [Load("UDL", 0, 2, 10), Load("UDL", 2, 4, 8)])
        self.assert_close(result.reactions[0].force, 22)
        self.assert_close(result.reactions[1].force, 22)
        self.assert_close(result.shear_moment(2)[1], 40)
        self.assert_close(result.shear_moment(5)[1], 67)
        self.assert_close(result.shear_moment(8)[1], 40)
        self.assert_close(result.deflection(5)[0], -8207/(12*1600))

    def test_applied_moment_jump(self):
        length, a, couple = 6, 2, 12
        result = self.simple(length, [Load("Moment", a, couple)])
        self.assert_close(result.reactions[0].force, couple/length)
        self.assert_close(result.reactions[1].force, -couple/length)
        before = result.shear_moment(a, "left")
        after = result.shear_moment(a)
        self.assert_close(before[1], couple*a/length)
        self.assert_close(after[1], -couple*(length-a)/length)
        self.assert_close(after[1]-before[1], -couple)
        self.assert_close(after[0]-before[0], 0)
        self.assert_close(result.extrema("M")[0][0], -8)
        self.assert_close(result.extrema("M")[1][0], 4)

    def test_point_force_jump_and_boundaries(self):
        result = self.simple(4, [Load("Point", 1, 10)])
        before = result.shear_moment(1, "left")
        after = result.shear_moment(1)
        self.assert_close(after[0]-before[0], -10)
        self.assert_close(after[1], before[1])
        self.assert_close(result.shear_moment(0, "left")[0], 0)
        self.assert_close(result.shear_moment(4)[0], 0)
        self.assert_close(result.shear_moment(4)[1], 0)

    def test_overhang_negative_reaction(self):
        result = solve_beam(6, 1600, [Support("Pinned", 2), Support("Roller", 5)], [Load("Point", 6, 10)])
        self.assert_close(result.reactions[0].force, -10/3)
        self.assert_close(result.reactions[1].force, 40/3)
        self.assert_close(result.shear_moment(5)[1], -10)
        self.assert_close(result.deflection(2)[0], 0)
        self.assert_close(result.deflection(5)[0], 0)

    def test_couple_at_free_cantilever_end(self):
        result = solve_beam(4, 1600, [Support("Fixed", 0)], [Load("Moment", 4, 10)])
        self.assert_close(result.reactions[0].force, 0)
        self.assert_close(result.reactions[0].moment, -10)
        self.assert_close(result.shear_moment(2)[1], 10)
        self.assert_close(result.deflection(4)[0], 10*4**2/(2*1600))

    def test_load_applied_directly_to_support(self):
        result = self.simple(4, [Load("Point", 0, 10), Load("Point", 4, -3)])
        self.assert_close(result.reactions[0].force, 10)
        self.assert_close(result.reactions[1].force, -3)
        self.assert_close(result.deflection(2)[0], 0)
        self.assert_close(result.shear_moment(2)[1], 0)

    def test_sign_reversal(self):
        down = self.simple(4, [Load("Point", 2, 10)])
        up = self.simple(4, [Load("Point", 2, -10)])
        self.assert_close(up.deflection(2)[0], -down.deflection(2)[0])
        self.assert_close(up.shear_moment(2)[1], -down.shear_moment(2)[1])

    def test_zero_load_model(self):
        result = self.simple(4, [])
        self.assertEqual(result.deflection(2), (0, 0))
        self.assertEqual(result.shear_moment(2), (0, 0))

    def test_ei_scales_deflection_only(self):
        first = self.simple(4, [Load("UDL", 0, 6, 4)], 800)
        second = self.simple(4, [Load("UDL", 0, 6, 4)], 1600)
        self.assert_close(first.deflection(2)[0], 2*second.deflection(2)[0])
        self.assert_close(first.reactions[0].force, second.reactions[0].force)

    def test_shear_extremum_inside_sign_changing_load(self):
        result = self.simple(4, [Load("Linear", 0, -4, 4, 4)])
        maximum, x = result.extrema("V")[1]
        self.assert_close(x, 2)
        self.assert_close(maximum, 4/3)

    def test_balanced_signed_linear_loads(self):
        result = self.simple(8, [Load("Linear", 0, -3, 4, 3), Load("Linear", 4, 3, 8, -3)], 12000)
        self.assert_close(result.reactions[0].force, 0)
        self.assert_close(result.reactions[1].force, 0)
        self.assert_close(result.equilibrium_error[0], 0)
        self.assert_close(result.equilibrium_error[1], 0)
        self.assert_close(result.shear_moment(2)[0], 3)
        self.assert_close(result.shear_moment(6)[0], -3)
        self.assert_close(result.shear_moment(4)[1], 8)
        self.assert_close(result.deflection(4)[1], 0)
        self.assert_close(result.deflection(4)[0], -224/(5*12000))
        self.assert_close(result.deflection(0)[1], -16/12000)

    def test_solution_contains_formula_and_current_model(self):
        result = self.simple(4, [Load("UDL", 0, 6, 4)])
        text = result.solution_text()
        self.assertIn("Kff df = ff", text)
        self.assertIn("L = 4 m", text)
        self.assertIn("R = 12", text)
        self.assertIn("EI y″ = M", text)

    def test_polynomial_root_isolation(self):
        roots = polynomial_roots([0.04, -0.4, 1])  # (x-.2)^2: repeated root
        self.assertEqual(len(roots), 1)
        self.assert_close(roots[0], .2)
        roots = polynomial_roots([-.08, .66, -1.5, 1])  # roots .2,.5,.8
        self.assertEqual(len(roots), 3)
        for actual, expected in zip(roots, (.2, .5, .8)):
            self.assert_close(actual, expected)


class InvalidModels(unittest.TestCase):
    def test_invalid_cases(self):
        simple = [Support("Pinned", 0), Support("Roller", 4)]
        cases = [
            (0, 1600, simple, []), (-4, 1600, simple, []),
            (4, 0, simple, []), (4, -1, simple, []),
            (float("nan"), 1600, simple, []), (4, float("inf"), simple, []),
            (4, 1600, [], []), (4, 1600, [Support("Pinned", 0)], []),
            (4, 1600, [Support("Pinned", 0), Support("Roller", 0)], []),
            (4, 1600, [Support("Fixed", 2)], []),
            (4, 1600, [Support("Fixed", float("nan"))], []),
            (4, 1600, [Support("Fixed", 5)], []),
            (4, 1600, [Support("Spring", 0)], []),
            (4, 1600, simple, [Load("Point", -1, 10)]),
            (4, 1600, simple, [Load("Point", 5, 10)]),
            (4, 1600, simple, [Load("Point", 2, float("nan"))]),
            (4, 1600, simple, [Load("Moment", float("inf"), 10)]),
            (4, 1600, simple, [Load("UDL", 2, 4, 2)]),
            (4, 1600, simple, [Load("UDL", 2, 4, 1)]),
            (4, 1600, simple, [Load("UDL", 0, 4, 5)]),
            (4, 1600, simple, [Load("Linear", 0, 4, 4, float("inf"))]),
            (4, 1600, simple, [Load("Unknown", 2, 4)]),
        ]
        for model in cases:
            with self.subTest(model=model), self.assertRaises(ValueError):
                solve_beam(*model)

    def test_invalid_read_position(self):
        result = solve_beam(4, 1600, [Support("Fixed", 0)], [])
        for x in (-1, 5, float("nan"), float("inf")):
            with self.subTest(x=x), self.assertRaises(ValueError):
                result.deflection(x)

    def test_extreme_numerical_model_reports_error(self):
        with self.assertRaises(ValueError):
            solve_beam(1e-200, 1e200, [Support("Fixed", 0)], [Load("Point", 1e-200, 1)])


if __name__ == "__main__":
    unittest.main()
