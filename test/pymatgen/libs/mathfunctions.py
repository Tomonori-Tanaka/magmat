import numpy as np
import sys

def rotvec(vec_in, mat, mode='N'):
    """
    Perform matrix x vector multiplication.

    vec_out = mat      * vec_in   (mode = 'N')
             (mat)^{t} * vec_in   (mode = 'T')

    vec_in: input vector. ndarray([3])
    mat: rotation matrix. ndarray(3x3)
    :return: vec_out: output vector.  ndarray([3])
    """
    vec_in = np.array(vec_in)
    vec_out = np.empty_like(vec_in)
    mat = np.array(mat)
    if mode=='N':
        # vec_out = mat@vec_in
        for i in range(3):
            vec_out[i] = mat[i][0] * vec_in[0] + mat[i][1] * vec_in[1] + mat[i][2] * vec_in[2]

    elif mode=='T':
        # vec_out = mat.T@vec_in
        for i in range(3):
            vec_out[i] = mat[0][i] * vec_in[0] + mat[1][i] * vec_in[1] + mat[2][i] * vec_in[2]
    else:
        print("Invalid mode ", mode)
        sys.exit(1)

    return vec_out



if __name__ == '__main__':
    vec_in = np.array([1, 2])
    mat = np.array([[1, 1], [0, 1]])
    vec_out = rotvec(vec_in, mat, 'T')
    print(vec_out)