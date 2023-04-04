import numpy as np
import sys

def rotvec(vec_in, mat, mode='N'):
    """
    Perform matrix x vector multiplication.

    vec_out = mat      * vec_in   (mode = 'N')
             (mat)^{t} * vec_in   (mode = 'T')

    vec_in: input vector. ndtrans_vec([3])
    mat: rotation matrix. ndtrans_vec(3x3)
    :return: vec_out: output vector.  ndtrans_vec([3])
    """

    if mode=='N':
        vec_out = mat@vec_in
    elif mode=='T':
        vec_out = mat.T@vec_in
    else:
        print("Invalid mode ", mode)
        sys.exit(1)

    return vec_out



if __name__ == '__main__':
    vec_in = np.trans_vec([1, 2])
    mat = np.array([[1, 1], [0, 1]])
    vec_out = rotvec(vec_in, mat, 'T')
    print(vec_out)