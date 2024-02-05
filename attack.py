import time
from fpylll import BKZ as BKZ_FPYLLL, GSO, IntegerMatrix, FPLLL
from fpylll.tools.quality import basis_quality
from fpylll.algorithms.bkz2 import BKZReduction
from random import randint
import json
import keygen
from keygen import NTRUKeyGenerator
from utils import (get_key_norm, get_norm, is_it_zero, is_it_ternary, is_it_pm_2, add_vectors_with_centerlifting,
                   substract_vectors_with_centerlifting, divide_by_2, run_all, parse_args, dump_blocksize_for_group,dump_seed)

FPLLL.set_precision(120)
"""
Implementation of the lattice reduction for both cyclic group(NTRU) and dihedral group(DiTRU)


How to:
   

Examples:

for dihderal group: 
    python attack.py 14  -q=128 --verbose=True --dump=True --group="dihedral" --h="[115, 42, 117, 108, 73, 3, 53, 29, 108, 34, 72, 5, 36, 101]"
    python attack.py 89 --verbose=True --group="cyclic" --dump=True --h="[403, 317, 342, 288, 441, 171, 102, 79, 33, 92, 295, 146, 160, 287, 480, 264, 167, 164, 396, 216, 493, 219, 351, 163, 140, 384, 62, 290, 218, 410, 151, 101, 369, 8, 398, 118, 231, 152, 428, 233, 172, 55, 307, 480, 58, 86, 469, 45, 10, 317, 121, 399, 234, 26, 108, 498, 325, 234, 37, 228, 456, 6, 371, 475, 310, 3, 22, 378, 419, 196, 93, 158, 124, 409, 286, 187, 216, 490, 71, 69, 122, 261, 388, 405, 71, 397, 343, 470, 337]" --filename="test" --bkz_betas="3:15"
"""

class Attack:
    """
    Attack on NTRU over cyclic group and dihedral group.
    We find the first block size that find the key for equivalent instances
    It includes the block size needed to find a {ternary key, non_ternary}

    Examples:

    for dihedral group:
            python attack.py 14  -q=128 --verbose=True --group="dihedral" --h="[115, 42, 117, 108, 73, 3, 53, 29, 108, 34, 72, 5, 36, 101]"
    """

    def __init__(self, params):
        self.n = params['n']  # The order of the group
        self.q = params['q']  # The used modulo
        self.seed = params['seed']  # NTRU key seed
        self.group = params['group']  # cyclic or dihedral
        #self.nsamples = params['nsamples']  # number of NTRU samples
        self.blocksizes = params['blocksizes']  # range [a,b] with starting blocksize a, last blocksize b-1
        self.ntours = params['tours']  # number of bkz tours
        self.nthreads = params['threads']  # number of threads
        self.verbose = params['verbose']  # verbose mode
        self.filename = params['filename'] ##file name
        self.dump     = params['dump']
        keyGenSuccess = False
        if params['h'] != None:
            self.h = json.loads(params['h'])  # for user-input h, we don't know f,g.
            self.f = None
            self.g = None
        else:
            self.h = None
        self.generator = NTRUKeyGenerator(self.n, self.q, self.seed, self.h, dihedral=(self.group == "dihedral"))

        while not keyGenSuccess:
            try:
                self.keyseed = self.generator.newSeed()
                #print("seed: ", self.seed)
                #print(self.generator.get_key(self.keyseed))
                self.lattice = self.generator.get_lattice(self.keyseed)
                # self.lattice contains a tuple: as (lattice, plus_lattice, minus_lattice)
                # in the case of the cyclic group both of plus_lattice and minus_lattice are None

                keyGenSuccess = True
            except:
                self.seed = self.seed+1
                self.generator.update_seed(self.seed)

        ### At this point we have the lattice of the element either in Z_qC_n or Z_qD_n
        ## Also the object generator is saving f,g corresponding to the lattice.
        self.f = self.generator.get_f()
        self.g = self.generator.get_g()
        self.h = self.generator.get_h()
        if self.group == "cyclic":
            self.dim = self.n * 2  # For cyclic group the lattice dimension for the reduction is two times the order
            self.threshold = 4 * get_key_norm(self.n)
        else:
            self.dim = self.n  # For dihedral group the lattice dimension where the reduction to applied equals the
            # order of the group
            self.threshold = 2 * get_key_norm(self.n)  # threshold is two times of the key norm, (check in two smaller
            # lattices) and then pull back to the original lattice, will
            # get almost 4*||key||
        if self.dim <= 178:
            self.float_type = "long double"
        else:
            self.float_type = "mpfr"

        self.basis = IntegerMatrix.from_matrix(self.lattice[0], int_type="long")
        self.M = GSO.Mat(self.basis, float_type=self.float_type,
                         U=IntegerMatrix.identity(self.basis.nrows, int_type=self.basis.int_type),
                         UinvT=IntegerMatrix.identity(self.basis.nrows, int_type=self.basis.int_type))

        # For dihedral group, we build save the two lattices {plus, minus}
        if self.group == "dihedral":
            self.plus_basis = IntegerMatrix.from_matrix(self.lattice[1], int_type="long")
            self.plus_M = GSO.Mat(self.plus_basis, float_type=self.float_type,
                                  U=IntegerMatrix.identity(self.plus_basis.nrows, int_type=self.plus_basis.int_type),
                                  UinvT=IntegerMatrix.identity(self.plus_basis.nrows,
                                                               int_type=self.plus_basis.int_type))

            self.minus_basis = IntegerMatrix.from_matrix(self.lattice[2], int_type="long")
            self.minus_M = GSO.Mat(self.minus_basis, float_type=self.float_type,
                                   U=IntegerMatrix.identity(self.minus_basis.nrows, int_type=self.minus_basis.int_type),
                                   UinvT=IntegerMatrix.identity(self.minus_basis.nrows,
                                                                int_type=self.minus_basis.int_type))


        #Open file and write the seed
           #print(self.M.B)
           # print()
           # print(self.plus_M.B)
           # print()
           # print(self.minus_M.B)
    def __call__(self):
        self.progressive_search()  # call the function that retrieves the key



    def check_for_cyclic(self, keys_found_tuple, keys):
        """
        Upon a reduced basis of a lattice for GR-NTRU based on the cyclic group,
        check for the ternary/non-ternary key.
        keys_found_tuple: a tuple refers if the (non-ternary-found, ternary-found).
        keys: it is passed as (None, None) initially, then will save the results of the found keys over rounds
        Output: a tuple (k1, k2) where ki itself is a tuple as ([f,g], norm).
        if no key is returned, returns "failure".
        """

        key1_found = keys_found_tuple[0]
        key2_found = keys_found_tuple[1]
        k1 = keys[0]  ##(key, norm)
        k2 = keys[1]  ##(key, norm)
        norms = {}
        # print(self.M.B)
        B = self.M.B
        for i in range(self.dim):
            norms[i] = B[i].norm()
        sorted_norms = sorted(norms.items(), key=lambda x: x[1])
        for i in range(self.dim):
            if sorted_norms[i][1] > self.threshold:
                if key1_found or key2_found:
                    return (k1, k2)
                return "failure"

            fg = list(B[sorted_norms[i][0]])
            f = fg[self.n:]
            g = fg[:self.n]
            if not is_it_zero(g):
                if not key1_found and self.generator.is_invertible_R_p(f):
                    k1 = (fg, sorted_norms[i][1])  # (key, its norm)
                    key1_found = True

            if not key2_found and is_it_ternary(fg):
                if self.generator.is_invertible_R_p(f):
                    k2 = (fg, sorted_norms[i][1])  # (key, its norm)
                    key2_found = True

            if key1_found and key2_found:
                return (k1, k2)

        return "failure"

    def check_for_dihedral(self, keys_found_tuple, keys):
        """
                Upon a reduced basis of a lattice for GR-NTRU based on the dihedral group,
                check for the ternary/non-ternary key.
                The reduced basis are two matrices of dimension n= 2N, where n is the order
                of the dihedral group.
                keys_found_tuple: a tuple refers of the (non_ternary found, ternary found).
                keys: it is passed as (None, None) initially, then will save the results of the found keys over rounds.
                Output: a tuple (k1, k2) where ki itself is a tuple as ([f,g], norm).
                if no key is found, returns "failure".
        """
        key1_found = keys_found_tuple[0]
        key2_found = keys_found_tuple[1]
        k1 = keys[0]
        k2 = keys[1]
        N = int(self.dim / 2)  # half the order of the dihedral (should be a prime).
        norms_plus = {}
        norms_minus = {}
        B_plus = self.plus_M.B  # reduced basis for the plus mat
        B_minus = self.minus_M.B  # reduced basis for the minus mat
        for i in range(self.dim):
            norms_plus[i] = B_plus[i].norm()
        sorted_norms_plus = sorted(norms_plus.items(), key=lambda x: x[1])

        for i in range(self.dim):
            norms_minus[i] = B_minus[i].norm()
        sorted_norms_minus = sorted(norms_minus.items(), key=lambda x: x[1])

        for i in range(self.dim):
            if sorted_norms_plus[i][1] > self.threshold:
                if key1_found or key2_found:
                    return (k1, k2)
                return "failure"
            t = list(B_plus[sorted_norms_plus[i][0]])
            #print("N: ", N)
            g0 = t[0:N]
            f0 = t[ N:2*N]

            if not is_it_zero(g0):
                for j in range(self.dim):
                    if sorted_norms_minus[j][1] > self.threshold:
                        break
                    t = list(B_minus[sorted_norms_minus[j][0]])
                    g1 = t[0:N]
                    f1 = t[N:2* N]

                    if not is_it_zero(g1):

                        fp0 = add_vectors_with_centerlifting(f0, f1, N, self.q)
                        fp1 = substract_vectors_with_centerlifting(f0, f1, N, self.q)

                        gp0 = add_vectors_with_centerlifting(g0, g1, N, self.q)
                        gp1 = substract_vectors_with_centerlifting(g0, g1, N, self.q)
                        #print("fp0: ", fp0)
                        #print("gp0: ", gp0)
                        F = fp0 + fp1  # concatenating (fp0, fp1)
                        G = gp0 + gp1  # concatenating  (gp0, gp1)

                        if not key1_found and self.generator.is_invertible_R_p(F):
                            #print("(f0,g0): ", f0+g0)
                            #print("first vecor norm: ", get_norm(g0 + f0))
                            #print("(f1,g1): ", f1+g1)
                            #print("second vecor norm: ", get_norm(g1 + f1))
                            k1 = (F + G, get_norm(F + G))
                            key1_found = True

                        if not key2_found and is_it_pm_2(F + G):
                            F = divide_by_2(F)
                            G = divide_by_2(G)

                            if self.generator.is_invertible_R_p(F):
                                k2 = (F + G, get_norm((F + G)))
                                key2_found = True

                            if key1_found and key2_found:
                                return (k1, k2)
        #print("reached here")
        return "failure"

    def progressive_search(self):
        """
         Apply reduction algorithm with increased block sizes and return the minimum blocksize
         that retrieves both a non-ternary and ternary keys
        """
        key1_found = False  # The non ternary key found?.
        key2_found = False # The ternary key found?.
        key1 = (None, None) ## (The non ternary key, its norm)
        key2 = (None, None) ## (The ternary key, its norm)
        key_tuple = [key1, key2]
        beta = [0]*2 ## block size needed to retrieve the (non-ternary key, ternary key)

        T0_global = time.time()
        if self.group == "cyclic":
            # print("before: ", self.M.B)
            self.bkz = BKZReduction(self.M)
            self.bkz.lll_obj()
            key_tuple = self.check_for_cyclic((key1_found, key2_found), key_tuple)
        else:
            # Apply LLL to both of the basis of the dihedral group
            self.bkz1 = BKZReduction(self.plus_M)
            self.bkz1.lll_obj()

            self.bkz2 = BKZReduction(self.minus_M)
            self.bkz2.lll_obj()
            key_tuple = self.check_for_dihedral((key1_found, key2_found), key_tuple)

        if self.verbose:
        # print("group:", self.group)
            fmt = "{'initial LLL applied: underlying group':'%8s', 'total walltime': %.3f}"
            print(fmt % (self.group, time.time() - T0_global))
        if key_tuple == "failure":
            key_tuple = [key1, key2]
            if self.verbose:
                print("failure")
        else:
            if key_tuple[0][0]!=None:
                if self.verbose:
                    print("(Non ternary key, its norm)", key_tuple[0])
                key1_found = True
                beta[0] = 2
            if key_tuple[1][0]!=None:
                if self.verbose:
                    print(("(Ternary key ,its norm)", key_tuple[1]))
                key2_found = True
                beta[1] = 2

        if not (key1_found and key2_found):

            for blocksize in self.blocksizes: ##apply bkz with increasing block size

                T0_local = time.time()
                if self.verbose:
                    print("New round with block size: ", blocksize)

                if self.group =="cyclic":
                    for t in range(self.ntours):  # runs BKZ tours
                         par = BKZ_FPYLLL.Param(blocksize,
                                                strategies=BKZ_FPYLLL.DEFAULT_STRATEGY,
                                                max_loops=8)
                         self.bkz(par)  ##apply bkz with the required parameters the check again for key
                         #print("new tour")
                         key_tuple = self.check_for_cyclic((key1_found, key2_found), key_tuple)
                         #print("key tuple after calling the function: ", key_tuple[0], "   ", key_tuple[1])
                         if key_tuple =="failure":
                             if self.verbose:
                                 print("failure")
                             key_tuple = [key1, key2]
                         else:

                             if key_tuple[0][0] != None and not  key1_found:
                                 if self.verbose:
                                     print("(Non ternary key, its norm)", key_tuple[0])
                                 key1_found = True
                                 beta[0] = blocksize
                             if key_tuple[1][0] != None:
                                 if self.verbose:
                                     print(("(Ternary key ,its norm)", key_tuple[1]))
                                 key2_found = True
                                 beta[1] = blocksize
                         if key1_found and key2_found:
                             break

                else:  # for dihedral group
                    for t in range(self.ntours):  # runs BKZ tours
                         par = BKZ_FPYLLL.Param(blocksize,
                                                strategies=BKZ_FPYLLL.DEFAULT_STRATEGY,
                                                max_loops=8)
                         self.bkz1(par)
                         self.bkz2(par)
                         key_tuple = self.check_for_dihedral((key1_found, key2_found), key_tuple)
                         if key_tuple == "failure":
                             key_tuple = [key1, key2]
                             if self.verbose:
                                 print("failure")
                         else:
                             if key_tuple[0][0] != None and not key1_found:
                                 if self.verbose:
                                     print("(Non ternary key, its norm)", key_tuple[0])
                                 key1_found = True
                                 beta[0] = blocksize
                             if key_tuple[1][0] != None:
                                 if self.verbose:
                                     print(("(Ternary key ,its norm)", key_tuple[1]))
                                 key2_found = True
                                 beta[1] = blocksize
                         if key1_found and key2_found:
                             break
                if key1_found and key2_found:
                    break
                fmt = "{'BKZ 'beta': %2d, underlying group':'%8s', 'total walltime': %.3f}"
                if self.verbose:
                    print(fmt % (blocksize,self.group, time.time() - T0_local))
        if self.verbose:
            print("Block size need to find (non-ternary key, ternary key) is ({},{})".format(beta[0], beta[1]))
        if self.dump:
            #print("non ternary key: ", key_tuple[0])
            #print("ternary key: ", key_tuple[1])
            dump_seed(self.seed,self.group,self.filename)
            dump_blocksize_for_group(self.f, self.g, self.h, key_tuple, beta, self.filename, self.group, self.seed, time.time()-T0_global)

def main_call(params):
    attack_inst = Attack(params)
    return attack_inst()


if __name__ == '__main__':
    print("main")
    all_parsed_params = parse_args()
    run_all(main_call, all_parsed_params)
