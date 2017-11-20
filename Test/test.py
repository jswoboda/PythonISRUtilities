import numpy.testing as npt
import isrutilities.mathutils as isrmath

def test_diric():
    x = isrmath.diric([5.,8.],2)
    npt.assert_allclose(x,[-0.80114362, -0.65364362])


if __name__ == '__main__':
 #   test_diric()
    npt.run_module_suite()
