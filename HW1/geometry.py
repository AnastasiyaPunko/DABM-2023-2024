import numpy as _np

def project_point_on_line(point, line_point_and_vector):
    r0 = _np.array(point)
    r1 = _np.array(line_point_and_vector[0])
    s1 = _np.array(line_point_and_vector[1])

    return r1 - ((r1 - r0).dot(s1.T)) / s1.dot(s1.T) * s1


def get_axis_rotation_matrix(vector, angle):
    if _np.linalg.norm(vector) != 1:
        vector = vector / _np.linalg.norm(vector)

    cf = _np.cos(angle)
    sf = _np.sin(angle)
    vf = 1 - cf
    x, y, z = vector

    return _np.array([[x * x * vf + cf, x * y * vf - z * sf, x * z * vf + y * sf],
                      [x * y * vf + z * sf, y * y * vf + cf, y * z * vf - x * sf],
                      [x * z * vf - y * sf, y * z * vf + x * sf, z * z * vf + cf]])

def place_next_dihedral_atom(prev_atoms, dih_ang_val, planar_ang_val, dist):
        A1 = prev_atoms[-2] + get_axis_rotation_matrix(prev_atoms[-1] - prev_atoms[-2],
                                                       dih_ang_val).dot(prev_atoms[-3] - prev_atoms[-2])
        v0 = (prev_atoms[-1] - prev_atoms[-2]) / _np.linalg.norm(prev_atoms[-1] - prev_atoms[-2])
        E = project_point_on_line(A1,
                                  (prev_atoms[-1],
                                   prev_atoms[-1] - prev_atoms[-2]))
        v1 = (A1 - E) / _np.linalg.norm(A1 - E)
        return prev_atoms[-1] - v0 * dist * _np.cos(planar_ang_val) + v1 * dist * _np.sin(planar_ang_val)