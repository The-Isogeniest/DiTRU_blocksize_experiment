from random import Random
import math
#from sage.matrix.matrix_space import MatrixSpace
#from  sage.rings.finite_rings.integer_mod_ring import  IntegerModRing
from sage.all import *
from utils import is_it_ternary
"""
Implementation of NTRU key generation for initial variant and NTREncrypt variant


How to:
    Create a NTRUKeyGenerator object with suitable parameters,
    generate seed via newSeed() and then call getKey() or get Lattice

Examples:

for dihderal group: 
    python attack.py 14  -q=128 --verbose=True --group="dihedral" --h="[115, 42, 117, 108, 73, 3, 53, 29, 108, 34, 72, 5, 36, 101]"

"""


class NTRUKeyGenerator:

    """
    Create NTRUKeyGenerator object with parameters
        -n: the order of the group
        -q: the required modulo
        -seed: the seed to generate a random bitstring

    """
    def __init__(self, n, q=0, seed=None, h=None, dihedral=False):

        self.n   = n
        self.q   = q
        self.p   = 3
        if h==None:
            self.h = None
        else:
            self.h = h
        self.cache = {}
        self.seed = seed
        self.sample_fixed_type = 30*n
        self.sample_key_bits = 2*self.sample_fixed_type
        self.d = math.floor(n/3)
        self.FFp = IntegerModRing(self.p)
        self.FFq = IntegerModRing(self.q)
        self.dihedral = False
        if dihedral:
            self.dihedral = True
            if self.n%2!=0:
                raise ValueError("The order of dihedral group is even!!!")
        if seed == None:
            self.seed = randint(0, 2 ** 64)

        self.f = None
        self.g = None # we are going to access f,g in the attack class just to pass them and print them in a file

    """
    Input: Integer s.
    Output: Random bit array of length s.
    """
    def randomBitArray(self,s):

        random = Random()
        random.seed(self.seed)
        return [random.randrange(2) for i in range(s)]
    """
    set the seed as a random bit array of length sample_key_bits.
    To be used as seed in getKey() and getLattice().
    """
    def newSeed(self):
        self.seed = self.randomBitArray(self.sample_key_bits)
        return self.seed

    """
    
    """

    def update_seed(self,new_seed):
        self.seed =new_seed
    """
    Input: two elements representing two polynomials from the ring Z_mod[x]/(X^n-1)
    Output: their multiplication
    """

    def ZCn_multiply(self, element1, element2, mod):
        multi_result = [0] * self.n
        # ai*aj*ri*rj = ai*aj(r_(i+j)%n)
        for i in range(self.n):
            for j in range(self.n):
                multi_result[(i + j) % self.n] = (multi_result[(i + j) % self.n] + element1[i] * element2[
                    j]) % mod  # ai*aj*ri*rj = ai*aj(r_(i+j))

        return multi_result
    """
    Input: two elements representing two polynomials from the ring  Z_qDN =~Z_q[x,y]/(x^N-1, y^2-1, yx-x^(N-1)y)
    Output: the result of multiplication 
    """
    def ZDn_multiply(self, element1, element2, mod):
        #print("inside")
        multi_result = [0] * self.n
        n = int(self.n/2) ##half the order of the group


        for i in range(n):
            for j in range(n):
                multi_result[(i + j) % n] = (multi_result[(i + j) % n] + element1[i] * element2[
                    j]) % mod

        for i in range(n):
            for j in range(n):

                multi_result[n + (j - i) % n] = (multi_result[n + (j - i) % n] + element1[i] * element2[
                    j + n]) % mod
        for i in range(n):
            for j in range(n):

                multi_result[n + (i + j) % n] = (multi_result[n + (i + j) % n] + element1[i + n] * element2[
                    j]) % mod
        for i in range(n):
            for j in range(n):
                multi_result[(-i + j) % n] = (multi_result[(-i + j) % n] + element1[i + n] * element2[
                    j + n]) % mod
        return multi_result

    """
           From https://ntru.org/f/ntru-20190330.pdf, Section 1.10.5.
           Input: A bit array b of length sample_fixed_type_bits.
           Output: A ternary polynomial with exactly d1 coefficients equal to 1 and d2  coefficients equal to −1.
    """
    def fixed_type(self, b, d1, d2):
        A = [0] * (self.n )
        v = [0] * (self.n )
        i = 0
        while i < d1:
            A[i] = 1
            for j in range(30):
                A[i] += 2 ** (2 + j) * b[30 * i + j]
            i += 1
        while i < d1+d2:
            A[i] = 2
            for j in range(30):
                A[i] += 2 ** (2 + j) * b[30 * i + j]
            i += 1

        while i < self.n:
            for j in range(30):
                A[i] += 2 ** (2 + j) * b[30 * i + j]
            i += 1

        A.sort()

        for i in range(self.n ):
            v[i] = A[i] % 4
            if v[i] ==2:
                v[i] =-1

        return v


    """
           Input: arr: an array , n: an integer 
           Output: shifting to left the array by n positions
    """

    def shiftLbyn(self, arr, n=0):
        return arr[n::] + arr[:n:]

    """
    
    Matrix representation for an element of Z_qC_n (right circulant matrix)
    It's also an auxiliary matrix for matrix representation for an element in Z_qD_n
    Input: the first row that represents an element f,g, or h.
    FF: the space over it, the matrix to be constructed either IntegerModRing(3)
    or IntgerModRing(q)
    

    Output: the matrix that represents the group ring element for a cyclic group
    """

    def get_A(self, first_row, FF):

        n = len(first_row)
        a = first_row
        m = []
        for i in range(n):
            m.append(a)
            a = self.shiftLbyn(a, -1)

        MS2 = MatrixSpace(FF, n, n)
        B = MS2.matrix(m)

        return B


    """
    An auxiliary matrix used to represent an element in Z_qD_n
    Input: the first row that represents an element f, g or h.
    FF: the space over it, the matrix to be constructed either IntegerModRing(3)
    or IntgerModRing(q)

    Output: the matrix that represents the group ring element for a cyclic group
            the output is a left circulant matrix.
    
    """
    def get_B(self,first_row,FF):
        n = len(first_row)
        a = first_row
        m = []
        for i in range(n):
            m.append(a)
            a = self.shiftLbyn(a, 1)

        MS2 = MatrixSpace(FF, n, n)
        B = MS2.matrix(m)

        return B

    """
    Input: an element of Z_qD_n.
    FF: the space over it, the matrix to be constructed either IntegerModRing(3)
    or IntgerModRing(q).
    Output: the matrix representation of the element.
    """

    def Zdn_matrix(self, element,FF):

        n = int(self.n/2) # half the order for dihedral group.
        A = self.get_A(element[:n],FF)
        #print("A: ", A)
        B = self.get_B(element[n:],FF)
        #print("B: ", B)
        M = block_matrix(2,2,[A,B,B,A])
        ApB = A+B
        AmB = A-B
        return (M, ApB, AmB)

    """
       Returns the matrix corresponding to the element
       For cyclic group is the matrix A.
       For dihedral group is the matrix [A,B,B,A].
    """
    def element_to_matrix(self,element,FF):
        if self.dihedral:
            return self.Zdn_matrix(element,FF)
        else:
            return (self.get_A(element,FF), None,None)


    """
    Sample f, g corresponding to the initial variant of NTRU defined in  Hoffstein  book 

    Input: A bit array fg_bits of length sample_key_bits.
    Output: two vectors represent f,g.
    """

    def sample_fg(self, fg_bits):

        f_bits = fg_bits[0:self.sample_fixed_type]
        g_bits = fg_bits[self.sample_fixed_type:]

        f = self.fixed_type(f_bits, self.d+1, self.d)
        Fp_mat, _, _ = self.element_to_matrix(f,self.FFp) ##eithe Z_qC_n or Z_qD_n matrix depending on the group
        if Fp_mat.is_invertible():
            Fq_mat, _,_ = self.element_to_matrix(f, self.FFq)
            if Fq_mat.is_invertible():
                Fp = Fp_mat.inverse()[0]  #inverse of f mod (p, X^n-1)
                Fq = Fq_mat.inverse()[0]  #inverse of f mod (q, X^n-1)
                g = self.fixed_type(g_bits, self.d, self.d)
                return (f,g, Fp, Fq)
        ## We reach here if the inverse doesn't exist for the previous seed
        ## Therefor, we generate a new seed
        #print("generating a new seed")
        #self.seed = randint(0, 2 ** 64)
        raise ValueError("Not invertible key")

    """
    Input: an element of the underlying group-ring.
    The function builds the matrix of the group ring element and 
    return True if it's invertible, otherwise, it returns False.
    """
    def is_invertible_R_p(self,element):
        Fp_mat, _, _ = self.element_to_matrix(element, self.FFp)  ##eithe Z_qC_n or Z_qD_n matrix depending on the group
        if Fp_mat.is_invertible():
            return True
        return False


    """
        Input: seed
        Output: a key (f,g,h) where h = g*f^-1 mod(q, X^n-1)
    """
    def get_key(self, seed):
        if self.h!=None:  ##h has been initialized in the constructor
            return (None, None,self.h)

        seedT = tuple(seed)
        if seedT in self.cache:
            return self.cache[seedT]

        else:
            f,g,Fp, Fq = self.sample_fg(seed)
            if self.dihedral:
                #print("here")
                h = self.ZDn_multiply(Fq,g, self.q)
                #print("h",h)
            else:
                h = self.ZCn_multiply(Fq,g,self.q)
        self.cache[tuple(self.seed)] = (f,g,h)
        return (f,g,h)

    """
    Input: seed
    Output: Coppersmith-Shamir basis that define for a cyclic group/Dihedral group
    For a dihedral group the output contains three lattices : The original lattice 
    and the other two lattices for reduction {plus, minus} lattices.
    """

    def get_lattice(self, seed):
        """
        Generate Coppersmith-Shamir lattice
        """
        n = self.n
        q = self.q
        plus_basis = None
        minus_basis = None
        f, g, h = self.get_key(seed)

        # we will save the values of f, g, h to access them from the attack file
        self.f = f
        self.g = g
        self.h = h
        B = [[0] * 2 * n for i in range(2 * n)]
        for i in range(n):
            B[i][i] = q
        for i in range(n):
            B[n + i][n + i] = 1

        element_mat, plus_mat, minus_mat = self.element_to_matrix(h,self.FFq)
        for i in range(n):
            for j in range(n):
                B[n + i][j] = int(element_mat[i][j])

        if self.dihedral:
            # Construct plus basis
            plus_basis= [[0] * n for i in range(n)]
            half_n = int(n/2)
            for i in range(half_n):
                plus_basis[i][i] = q
            for i in range(half_n):
                plus_basis[half_n + i][half_n + i] = 1

            for i in range(half_n):
                for j in range(half_n):
                    plus_basis[half_n + i][j] = int(plus_mat[i][j])


            # Construct minus basis
            minus_basis = [[0] * n for i in range(n)]
            for i in range(half_n):
                minus_basis[i][i] = q
            for i in range(half_n):
                minus_basis[half_n + i][half_n + i] = 1

            for i in range(half_n):
                for j in range(half_n):
                    minus_basis[half_n + i][j] = int(minus_mat[i][j])

        return (B, plus_basis, minus_basis)

    """
    The function returns the value of f
    """
    def get_f(self):
        return self.f

    """
    The function returns the value of g
    """
    def get_g(self):
        return  self.g


    """
    The function returns the public key h
    """
    def get_h(self):
        return self.h

    """
    Input: seed
    Output: Reduction lattices
            - in the case of cyclic group, it's one lattice
            - For dihedral group, two lattices of dimension 2*n
    """

    def get_reduction_lattices(self, seed):
        n = self.n
        q = self.q
        f, g, h = self.get_key(seed)
def main():

    n = 7
    q = 64
    print(is_it_ternary([1, -1, 0, -1, 1, 0, 0, -1, 0, 0, 1, 1, 0, 0]))
    for i in range(100):

        seed = randint(0,2**64)
        keygen = NTRUKeyGenerator(n, q,seed=seed)

        keyseed = keygen.newSeed()
        f,g,h= keygen.get_key(keyseed)
        print(h)
        print(g)



if __name__ == "__main__":
    main()