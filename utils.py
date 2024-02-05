import argparse
import six
import copy
import sys, os
import re
from math import sqrt, floor, ceil, log, exp, log2
from datetime import datetime
import random
from random import randint
from multiprocessing import Pool
import matplotlib.pyplot as plt
from random import randint
from fpylll.util import gaussian_heuristic
from fpylll import IntegerMatrix
from mpmath import mp, erfc as mperfc
from scipy.special import  erfcinv
import csv
mp.dps = 1200
error_rate = 2**-64
import warnings



max_beta = 150 # maximal dimension we run SVP in
cross_over_beta = 70

def remove_zeros(B):
    """
    removes zero rows from matrix B
    """
    cr = 0
    for i in range(B.nrows):
      if not B[i].is_zero(): cr+=1

    B_ = [0]*cr
    cr = 0
    for i in range(B.nrows):
      if not B[i].is_zero():
        B_[cr] = B[i]
        cr+=1

    return IntegerMatrix.from_matrix(B_, int_type="long")

def testAndMakeDir(path):
  if not os.path.isdir(path):
      os.makedirs(path)


def get_norm(vector):
    """
    Calculates and returns the norm of the vector.
    """
    # print("inside the function: \n")
    # print("iside vector: ", vector)

    s = 0
    for v in vector:
        s += v * v
    # print("inside the norm: ", np.sqrt(s))
    return sqrt(s)

def get_key_norm(n):
    """
    Input: n the order of the group
    Output: the norm of the key of the form (f,g)
    where f in T(d+1,d) and g in T(d,d).
    """

    d = int(n/3)
    return sqrt(4*d+1)

def is_it_ternary(l):
    """
    Input: a list l.
    Check if the list is ternary and returns True otherwise returns True
    """
    for i in l:
        if i!=1 and i!=-1 and i!=0 :
            return False
    return True
def is_it_zero(l):
    """
    Input: list l
    Output: True if the list entries are all zeros, False otherwise
    """
    for i in l:
        if i!=0:
            return False
    return True

def get_q_no_error(d,p):
    """
    The function returns the value of q that gives no decryption failure for variant of NTRU
    that has: h = gf^-1
    Input: d = int(order of the group/3)
           p usually 3
    """
    value= p*(6*d+1)
    q= 2**(len(bin(value))-2)
    return q

"""
Calculating the decryption failures in the following functions is happening according to 
calculations in the paper https://link.springer.com/chapter/10.1007/978-3-642-01957-9_27 .
The calculations are returned for two variants of NTRU
            - The initial variant of NTRU where h = g*f^-1 (we call it old)
            - Another variant of NTRU where h = (3f+1)^-1*g
            Also the estimation is being calculated for both cyclic group and dihedral
"""


def sigma_square(n, dihedral=False, old=True, p=3):
    """
    input:  - n the order of the group
            - dihedral: refers if the dihedral group is used
            - old: 1 means the old NTRU design mentioned in Hoffstein book f^-1*g otherwise (3f+1)^-1*g design
            - p: the modulo.
    output: returns sigma^2 for Xi
    where Xi is a random variable refers to coefficients of the a = f*c where c is the ciphertext
    For the old implementation of NTRU a = prg+fm where r,g in T(d,d) and f in T(d+1, d). [any coefficient must be smaller
    than q/2 for correct decryption]
    For the second implementation a = prg+(3f+1)*m for f, r, g as mentioned above.
    Therefore any coefficient in  rg+fm should be smaller than (q-1)/2p
    """
    if dihedral:
        n = int(n/2) # r which is sampled from T(dr_dihedral, dr_dihedral) can be thought as r0+y*r1 where r0 sampled
                     # from T(dr_cyclic, dr_cyclic) that's why we divide by 2.
    dr = int(n/3)
    df = int(n/3)
    m = 1
    if dihedral:
        m = 2
    t = 1
    #print("old: ", old)
    if old:
        t = p**2
    ss = m * (4 * (t * dr + df) / 3)
    #print("ss: ",ss)
    return ss

def decryption_failure(n, c, sigma):
    """
    Input: n: the order of the group.
           c: the threshold value (greater than this value a decryption failure occurs)
           sigma of the coefficients' distribution.
   Output: decryption failure
    """
    return n*mperfc(c/(sigma*sqrt(2)))


def get_q(n, sigma, error_rate, old=True, dihedral=False, p=3,powoftwo=1):
     """
     Input  - n: the order of the group.
            - sigma of the coefficients' distribution.
            - error_rate: the accepted error rate (any error rate smaller is accepted).
            - old:  - old: 1 means the old NTRU design mentioned in Hoffstein book f^-1*g otherwise (3f+1)^-1*g design.
            - dihedral: if 1 means dihedral group otherwise cyclic group.
            - p: the modulus
            - powoftwo: means return the first power of two that statisfies the required error rate
              if 0 returns the first q not the power of two.
     Output: the first power of two that gives decryption failure smaller than the input error-rate.

     """
     #print("n: ", n)
     #print("sigma: ", sigma)
     #print(error_rate)
     #print("dihedral: ", dihedral)
     #if dihedral:
     #    n = int(n/2) ## We don't need to divide by 2 for dihedral
     c = erfcinv(error_rate/n)*(sigma*sqrt(2))
     #print("c = ",c)
     if old:
         q_orginal = 2*c # q should be greater than this value
     else:
         q_orginal = 2*p*c+2

     q_orginal = ceil(q_orginal)
     q = q_orginal
     if powoftwo:
         q = 2 ** (len(bin(int(q_orginal))) - 2)

     return q
def is_it_pm_2(l):
    """
    Input: a list.
    Return True if all entries are two, minus two, or zeros.
    Otherwise: False.
    """
    for i in l:
        if i!=2 and i!=-2 and i!=0 :
            return False
    return True


def divide_by_2(l):
    """
    Input: a list of {2,-2,0}
        divide the coefficients by 2 and return the resultant list.

    """
    for i in range(len(l)):
        if l[i]>0:
            l[i] =1
        elif l[i]<0:
            l[i] = -1
    return l
def add_vectors_with_centerlifting(l1, l2, n, q):
    """
    Coefficients-wise adding of coefficients of the correspondence vectors
    with centrelifting
    Input: l1, l2 two lists representing two vectors in the lattice.
    n: the vector length.
    Output: the resultant vector after adding them.
    """
    res = [0]*n
    for i in range(n):
        res[i] = (l1[i]+l2[i])%q
        if res[i]>int(q/2):
            res[i] = res[i]-q
    return res


def substract_vectors_with_centerlifting(l1, l2, n, q):
    """
    Coefficients-wise adding of coefficients of the correspondence vectors
    with centrelifting
    Input: l1, l2 two lists representing two vectors in the lattice.
        n: the vector length
    Output: the resultant vector after adding them.
    """
    res = [0] * n
    for i in range(n):
        res[i] = (l1[i] - l2[i]) % q
        if res[i] > int(q / 2):
            res[i] = res[i] - q
    return res





def dump_seed(seed, group,filename):
    """
    Input: the seed and the file name
    Output: write the seed to the file to later add the trails and the betas.
    """
    org_seed = seed
    seed = seed - (seed % 10 ** 8)
    path = "keys_dumps/" + group + "/" +"seeds/"
    testAndMakeDir(path)

    filename += "_" + str(seed) + ".txt"
    with open(path + filename, "a+") as f:
        print("seed: ",org_seed, file=f)



def process_sample(sample):
    """

    Input: a sample as  a list of info
    a sample can be a list holding infor like: ['f', 'g', 'key norm', 'h', 'k1 (non-ternary)', 'k1-norm', 'k2 (ternary)', 'k2-norm', 'beta1', 'beta2', 'total time (seconds)']

    The function processes the sample and returns a list of strings values corresponding to the values in the list.

    """
    l = []

    for i in range(len(sample)):

        s = str(sample[i])
        string_to_write = ""
        for i in range(len(s)):
            if (s[i] != "]" and s[i] != "["):
                string_to_write += s[i]
        l.append(string_to_write)
    return l


def create_file(seed, group,filename):
    """
    Input: the seed, the group and the file name
    The function creates a file with the specified path and write the header to the file.
    the header = [f, g, norm, h, f1_prime, f1_norm, f_prime2, f2_norm, beta1, beta2, total_time]
    The function writes the header into a csv file, if not already existed.
    """
    #print("file created: ")

    seed = seed - (seed % 10 ** 8)
    path = "keys_dumps/" +group + "/records/"
    testAndMakeDir(path)
    header = ['f', 'g', 'key norm', 'h', 'k1 (non-ternary)', 'k1-norm', 'k2 (ternary)', 'k2-norm', 'beta1', 'beta2', 'total time (seconds)']

    filename += "_" + str(seed) + ".csv"
    isExisting = os.path.exists(path+filename)
    if not isExisting:
        with open(path + filename, "w", newline='') as wfl:
           csvwriter = csv.writer(wfl, delimiter=',')
           csvwriter.writerow([val for val in header])

def dump_blocksize_for_group(f,g,h,key_tuple, beta, filename, group, seed, total_time):
  """
  Input:
         f,g,h: the original key (f,g): is the private key and h is the corresponding public key.
         key_tuple: it is an array key_tuple[0] = (non-ternary-key, its norm) and key_tuple[1] = (ternary key, its norm)
         key_tuple can be a string "failure" if we couldn't find the key.
         beta: an array: beta[0]: the blocksize needed to find the non-ternary key and beta[1]: the blocksize needed to
         find the ternary key.
         filename: the file name (recommended to be "n_q")
         group: cyclic or dihedral.
         seed: the seed.
         total)time: the total time the attack took to find the keys
  """
  sample = []
  sample.append(f)
  sample.append(g)
  if f== None:
      sample.append("NA")
  else:
      sample.append(get_norm(f+g))
  sample.append(h)

  if key_tuple =="failure":
      sample.append("failure") ## not able to find the non-ternary key
      sample.append("NA")  ## no norm
      sample.append("failure") ## not able to find the ternary key
      sample.append("NA") ## no norm
      sample.append("failure for tried betas") ## beta for the non ternary key
      sample.append("failure for tried betas") ## beta for the ternary key
  else:
      sample.append(key_tuple[0][0]) ## non ternary key
      sample.append(key_tuple[0][1]) ## non-ternary norm
      sample.append(key_tuple[1][0]) ## ternary key
      sample.append(key_tuple[1][1]) ## ternary norm
      sample.append(beta[0]) ## beta1 needed to find the non-ternary key
      sample.append(beta[1]) ## beta2 needed to find the ternary key

  sample.append(total_time)

  seed = seed - (seed % 10 ** 8)
  path = "keys_dumps/" + group + "/records/"
  testAndMakeDir(path)

  to_write = process_sample(sample)
  filename += "_" + str(seed) + ".csv"
  with open(path + filename, "a+", newline='') as wfl:
      csvwriter = csv.writer(wfl, delimiter=',')
      csvwriter.writerow([val for val in to_write])
      # print( str(beta1) + "\t" + str(beta2) + "\t" + str(total_time) + "\t"+ datetime.now().strftime('%Y-%m-%d %H:%M:%S'), file = f )


def rough_estimate_on_betas(n,q, dihedral=False):
    """
    use the output of this function in case the use did not provide us with blocksizes
    """
    if dihedral:
        n = int(n/2)
    if n<150: beta_low = 10
    else: beta_low = floor( 0.28*4.*n/( (log(q)/log(n))**2 + 1))
    return list(range(beta_low, max_beta))




def parse_args():
    parser = argparse.ArgumentParser(description='Parse NTRU attack params.')

    #main parameters
    parser.add_argument('n', type=int, help="ring dimension")
    parser.add_argument('--group', type=str, help="cyclic or dihedral", default="cyclic")
    parser.add_argument('-q', type=int, dest="q", default=None, help="NTRU modulus")
    parser.add_argument('--nsamples', type=int, default=None, dest="nsamples", help="Number of samples/rows of rot(h) used")
    parser.add_argument('--seed',  type=int, dest="seed", default=None, help="randomness seed")
    parser.add_argument('--h', dest="h", default=None, help="Uses given input as h, instead of creating a random instance.")
    parser.add_argument('--dump', dest='dump', default=False, help="flag to dump intermediate bases")
    # number of runs, number of threads
    parser.add_argument('-t', '--trials', type=int, dest="trials", default=1,
                        help="number of experiments to run per dimension")
    parser.add_argument('-w', '--workers', type=int, dest="workers", default=1,
                        help="number of parallel experiments to run")
    parser.add_argument('--threads', type=int, dest="threads", default = 1, help="number of threads used by 1 worker")


    parser.add_argument('--bkz_betas', type=str, dest="blocksizes", default=None, help="bkz block sizes as string of the form: min_beta:max_beta:step")
    parser.add_argument('--bkz_tours', type=int, dest="tours", default=8, help="number of tours of bkz reduction")


    parser.add_argument('--verbose', dest="verbose", default=False, help="verbosity")
    parser.add_argument('--dry-run', dest="dry_run", default=False,
                        help="Show parameters that would be used but don't run any actual experiments.")
    parser.add_argument('--show-defaults', dest="show_defaults", action='store_true',
                        help="Show default parameters and exit.")

    parser.add_argument('--filename', dest='filename', default=None, help="prefix of the dump filenames")



    args, unknown = parser.parse_known_args()


    fmt = "{key:%ds}: {value}"%20

    if len(unknown)>0:
        print('Parameters', unknown, 'are not recognized and will be ignored')

    all_defaults = {key: parser.get_default(key) for key in vars(args)}

    if args.show_defaults:
        for k, v in six.iteritems(all_defaults):
            print(fmt.format(key=k, value=v))
        exit(0)

    all_params = check_parsed_params(vars(args))

    if args.dry_run:
        for k, v in six.iteritems(all_params):
            print(fmt.format(key=k, value=v))
        exit(0)


    return all_params

def check_parsed_params(params):

    if not params['group'] in ["cyclic", "dihedral"]:
        raise ValueError("group =%s not recognized." % params['group'])


    if params['nsamples']==None: params['nsamples'] = params['n']
    else: assert(params['nsamples'] > 0 and params['nsamples']<=params['n'])



    if params['blocksizes'] ==None:
        params['blocksizes'] = rough_estimate_on_betas(params['n'], params['q'], dihedral= params['group']=="dihedral")
    else: params['blocksizes'] = eval("range(%s)" % re.sub(":", ",", params['blocksizes']))


    assert(len(params['blocksizes'])>0)



    if params['seed']==None:
        params['seed'] = randint(0, 2**64)

    if params['q']==None:
        ## for dihedral group, calculate for cyclic the power of two that gives
        ## error less than 2**-100, then q' = sqrt(2)*q for dihedral.

        n = params['n']
        if params['group'] == "dihedral":
            n = int(n/2)

        ss = sigma_square(n) ##cyclic
        sigma = sqrt(ss)   #cyclic
        q =get_q(n, sigma, error_rate) ### for cyclic always here
        if params['group'] == "dihedral":
            new_rate = decryption_failure(n, int(q/2), sigma) ##get the exact error rate for the cyclic for the selected q
            ssd = sigma_square(params['n'], dihedral=True)  ##dihedral
            sigma = sqrt(ssd)  # cyclic
            q = get_q(params['n'], sigma, float(new_rate), powoftwo=False)
            print()
        params['q'] = q

    if params['filename']==None:
        params['filename'] = str(params['n'])+'_'+str(params['q'])
    return params

def run_all(f, params):
    jobs = []

    original_seed = params['seed']
    if params['dump']:
        create_file(original_seed, params['group'], params['filename']) ##Create excel file to start saving the records
        #dump_seed(original_seed,params['group'],params['filename'])
    for t in range(params['trials']):
        params_  = copy.deepcopy(params)
        params_['seed'] = original_seed+random.randint(1,1000)
        jobs.append(params_)
    if params['workers'] == 1:
        for job in jobs:
            res = f(copy.deepcopy(job))
    else:
        pool = Pool(params['workers'])
        pool.map(f, jobs)
