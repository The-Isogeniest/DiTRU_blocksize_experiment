# DiTRU: A Resurrection of NTRU over Dihedral Group
This repository contains the scripts accompanying the article

## Providing experimental results to identify the hardness of solving the SVP against NTRU lattice over cyclic group vs. DiTRU lattice after one layer of Gentry's attack



# Requirements

* [fpylll](https://github.com/fplll/fpylll)

* [SageMath 9.5+](https://www.sagemath.org/) 


# Description of files
Short description of the content:
* `attack.py` main file to run lattice reduction attack on NTRU over cyclic group and DiTRU (i.e., GR-NTRU over dihedral group)
* `keygen.py` contains the functions for key generation for both NTRU and DiTRU
* `utils.py` contains helping functions
* folder `keys_dumps` contains two folders, `cyclic` and `dihedral`: each of them contains two subfolders, `records` and `seeds.`
  * subfolder `records` contains records corresponding to the original keys and the retrieved keys according to the attack
    ; for a cyclic group, the records save `f`, `g`, `key norm,` and `h` for the original key 
    `k1 (non-ternary)`: the non-ternary key retrieved by the attack, `k1-norm`: its norm,`k2 (ternary)`: the ternary key 
     retrieved by the attack, `k2-norm`: its norm, `beta1`: blocksize needed to retrieve k1, `beta2`: blocksize needed to
     retrieve the ternary key, `total time (seconds)`: total time of the attack.
  
    For the dihedral group: the folder contains the same data, but the collected after applying one layer of Gentry's attack
  * subfolder `seeds` contains the seeds for which the records have been generated; running the attack using the same 
    seed reproduces the same results.
# How to use

Run `attack.py` with the following parameters.

* `n`: defines the order of the underlying group. The order 'n' is prime for the cyclic group, and `n = 2N` for `N` is prime 
for the dihedral group. Integer. Obligatory parameter. For NTRU over the cyclic group, the SVP is solved over a
lattice of dimension `2n`, while for the DiTRU over the dihedral group, we perform one layer of Gentry's attack before 
solving the SVP over a lattice of dimension `2N`.
* `-q`: NTRU modulus. If not entered, for the cyclic group, it is calculated as the first power of two, which guarantees
that the decryption failure does not exceed the specified threshold. For DiTRU: it is calculated as the exact `q`
that achieves the same decryption failure compared to NTRU over the same dimension.
* `--seed`: randomness seed to generate a key and build the corresponding lattices.
* `--group`: `cyclic` for standard NTRU over the cyclic group, `dihedral` for DiTRU over the dihedral group.
* `--bkz_betas`: a string as 'b1:b2', where b1 indicates the first blocksize and b2 indicates the last blocksize
to try upon running progressive bkz.
* `--dump`: True to save the results into files, False otherwise.
* `--verbose`: True to print detailed output in the console while running the program, False otherwise.
* `--filename`: the file's name where to save the results.

Parameters related to running experiments in parallel: 
* `-t` (or `--trials`): number of experiments to run per GR-NTRU dimension 
* `-w` (or `--workers`): number of parallel experiments to run
* `--threads`: number of threads used by one worker


# Experiments

You can solve the SVP using progressive BKZ for a lattice of dimension 188 by running the command.
```
python attack.py 89   --verbose=True --dump=True --group="cyclic"   --bkz_betas=3:50

```

It takes approximately less than 40 seconds on a laptop.
It generates a random seed and an instance corresponding to the seed and runs the attack.

You can specify the seed to generate a specific example or reproduce one of our experiments.


```
python attack.py 89   --verbose=True --dump=True --group="cyclic"   --bkz_betas=3:50 --seed=3301211018778602037
```


You can also specify the number of trials to run.


```
python attack.py 89   --verbose=True --dump=True --group="cyclic"   --bkz_betas=3:50 --trials=100
```

For a dihedral group of order `2N`, first, we apply one layer of Gentry's attack, then the SVP is solved in two lattices 
of dimensions `2N`, and the solution is pulled back to the original lattice of dimension `2n=4N`.


###  For instance, to find a short vector in the lattice associated with DiTRU of order `n=254`, one can run the command


Running the command 
```
python attack.py 254   --verbose=True --dump=True --group="dihedral"  --trials=20 --bkz_betas=50:65 --seed=11230911219425525815
 
```
The previous command finds the following non-ternary vector in the DiTRU lattice.
```
0, 0, 3, 1, 1, 2, 0, 1, 0, -2, -1, 1, -3, 1, 1, 1, 1, 0, -1, -1, -1, 0, 2, 0, 0, -2, -2, 1, 2, 0, -2, 1, 0, 3, -2, -1, 
-2, -2, 3, 1, -1, 2, 1, 1, 1, -3, -2, -1, -1, 1, 1, -2, -2, -2, 0, -1, -1, 1, -1, -3, 1, 1, 0, 1, -2, -2, -1, 2, -2, -1, 
-1, -2, -1, 0, 0, 3, 0, 0, -2, 1, -1, 1, 0, -2, 0, -2, -3, -1, 0, -1, 1, -1, 0, 1, 0, -1, -1, 2, 0, 0, 1, 3, -2, 3, 1, 0
, 0, 0, 1, 0, 3, 1, 2, 0, -1, -2, 0, 2, 1, 0, 1, -1, 0, -2, 0, 2, 2, 0, 0, -1, -1, 1, 0, 0, -1, -4, 0, 1, 1, -1, -1, -1,
-1, 1, 0, -1, 3, -1, 0, 2, 0, -2, -2, 0, 1, -2, -2, 0, 1, 2, 1, 2, -1, 0, 0, -1, 1, -1, 2, -1, -3, 1, 1, 0, -1, 1, -1, 
-3, 2, 0, 0, -2, -1, 1, 1, 1, 1, -1, -1, -2, -1, -2, -2, 1, 2, 2, -1, 3, 0, -1, 0, 2, 1, 0, 0, 2, 1, 1, 3, 0, -2, -2, 0, 
1, 3, 4, -1, 3, -3, 0, 1, 0, 3, 1, -2, 0, 0, -1, -1, 2, -1, 3, 0, 0, 2, 1, 0, 1, -1, 0, 0, 1, -2, 2, 0, 1, -2, 1, -1, 0, 
-2, 0, -2, -2, 0, -1, 0, -1, 0, -1, 0, -3, -2, 1, -1, -3, 3, 0, 1, 0, -1, 0, -1, 0, 1, 1, 2, -2, -2, 1, 0, 2, 0, -1, 0, 
-2, -2, -1, -2, -2, 2, -2, -1, 1, -1, -1, -1, 1, -3, 2, 1, 2, -2, 1, 1, -2, -2, -1, -1, -2, 0, 0, 0, 1, -1, -3, -2, 2, 2
, -1, -2, -2, -1, -1, -2, 1, -1, 1, -2, 3, -1, -3, 1, 0, 0, 0, 3, -1, 1, -1, -1, 1, -2, 0, 1, 1, 2, 0, -1, 1, -2, -1, -1, 
1, -3, 2, 1, 1, -1, 2, 0, 1, -1, 2, -1, 0, 0, 3, 2, -2, 0, 0, 1, -1, -1, 0, -2, -4, -1, -2, 0, 0, -1, -2, 3, 0, 1, -2, -1,
0, -1, 1, -1, -1, -2, 3, 0, 1, 2, 1, 0, 1, -1, -2, 2, 2, 3, -4, 0, 0, -1, 0, 0, 2, 1, 0, 0, 0, 0, -1, -3, -1, 1, -1, 1, 
-1, 0, -1, 0, -2, -1, 1, 2, 0, 3, 3, 0, 2, 0, 2, -1, 3, 1, 2, 2, 0, -1, 0, 0, -1, 1, 2, 1, -1, 1, 0, 1, 1, 1, -1, 2, 0, 
-2, 1, 1, 1, -1, -1, 3, 0, 0, -1, 3, 2, -2, 1, 1, 2, -1, 3, 3, -1, -2, 1, -1, -1, 2, 0, 3, 1, -2, -1, 0, -4, -1, 2, 2, 
0, 2, 1, 3, 1, -2, 0, 0, -1, -2, 4
```
The norm of the previous key is `34.5253530032641`; the attack further searches for the ternary key and finds
the following vector in the associated lattice.

```
-1, 0, -1, 1, 1, -1, 0, 0, 0, -1, -1, 0, 1, 0, 0, -1, -1, -1, -1, -1, 0, 1, -1, 1, -1, -1, 0, 1, -1, 0, 0, 0, -1, 0, -1, 
1, 1, -1, 1, -1, 0, -1, 0, 0, -1, -1, 0, 1, 1, 1, 0, 1, 0, 0, 1, -1, -1, 1, -1, 0, -1, -1, 1, 0, 1, 0, -1, -1, 1, 0, 0, 
-1, -1, 1, 0, 0, 1, -1, 0, 0, -1, 1, 0, 0, 1, 1, 1, 1, 1, 0, -1, 1, 0, 1, 0, 1, -1, 1, 1, -1, 1, 0, 1, 1, 1, 1, -1, 1, 
0, 0, 1, 0, 1, 0, -1, 0, 1, 0, -1, 1, -1, 0, 1, 0, 1, 0, -1, -1, 1, -1, -1, 1, -1, 0, 0, 0, -1, -1, 1, 1, 0, 1, 0, 1, -1
, 0, -1, 0, 1, -1, 1, -1, 1, 1, -1, -1, -1, 1, -1, 0, 1, -1, -1, -1, 0, -1, 0, 0, -1, -1, 0, 1, 1, 0, 1, 0, -1, 1, -1, 1,
-1, 0, 0, -1, 0, 1, 1, 0, -1, -1, 0, 0, -1, 1, 1, 0, -1, -1, -1, -1, 0, 1, 0, 1, 0, -1, -1, 1, 0, 1, 1, 1, -1, 1, -1, 1,
0, 0, 0, -1, 0, 1, 0, 1, 0, -1, -1, 1, 0, 0, -1, 0, 1, 0, 0, 1, 1, -1, -1, 1, 1, 1, 0, 1, 0, 1, -1, 0, -1, -1, 1, 1, 0, 
0, 1, 1, -1, -1, -1, 1, 1, 0, 0, 1, 0, -1, 1, 1, 0, -1, 1, -1, 1, -1, 1, 1, 1, -1, 1, 0, -1, 0, 0, -1, 0, 0, 1, 1, 1, 1,
0, 1, -1, -1, 1, 1, 0, -1, -1, 0, 1, 0, 1, 1, 1, 0, 1, -1, 1, 0, -1, 1, 0, 0, 1, 1, -1, 0, -1, 1, -1, 0, 1, 1, 1, 1, 1, 0, 
0, 0, 0, 1, 1, 1, 0, 0, -1, 1, 0, 1, 0, 1, 0, -1, -1, 1, -1, 0, 1, 0, 0, -1, 0, -1, 0, -1, -1, -1, 0, 1, 1, 0, 0, -1, -1, 
1, 0, 1, -1, 1, 1, 1, -1, 1, 0, -1, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0, 1, -1, -1, -1, 1, 0, 0, 1, 0, 0, -1, 0, -1, 0, 0, 0, 0,
1, 1, -1, 0, 1, 0, 0, 1, 0, 1, 0, 0, -1, 0, -1, -1, 0, -1, 0, -1, 1, -1, -1, -1, 0, 1, -1, 0, -1, 1, -1, -1, 0, -1, -1, 0, 
0, -1, -1, 1, -1, -1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1, -1, -1, 0, 0, 0, -1, -1, -1, 0, 0, 1, 1, 1, -1, 1, 0, 0, -1, 1, 1,
1, 1, 1, -1, -1, 1, -1, 0, -1, -1, 0, 0, -1, 0, -1, -1, -1, 1, 0, -1, 1, 0, 0, 0, -1, 0, -1, -1, 1, 1, -1, 1, 0
```
The norm of the found ternary key is `18.3575597506858`, and both of the found keys have been found at blocksize=63 and 
took approximately`30.7` core days on a system with Linux (Ubuntu 22.04.2
LTS) on Intel(R) Xeon(R) CPU E3-1246 v3 @ 3.50GHz and 32 GB installed RAM.

The following table summarizes the experimental results for the tested dimensions; the dimension of the lattice where the SVP is solved is `2N`.
The blocksize value is averaged over `100` trials where the lattice's dimension is smaller than `226` and averaged over at least `20` trials for larger dimensions.

  
| N                   | 71 | 73 | 79 | 83 | 89 | 97 | 101 | 107 | 109 | 113 | 127 | 131|
| ------------------  |----|----|----|----|----|----|-----|-----|-----|-----|-----|----|
| $\beta (cyclic)$    |2.28|2.48|3.02|3.64|5.22|8.94|11.05|16.22|18.56|26.68|52.8|57.35|
| $\beta (dihedral)$  |2.87 |3  |3.88|5.06|7.18|11.68|15.63|26.57|34.06|43.95|63.6|_   | 

![summary](https://github.com/The-Isogeniest/DiTRU_blocksize_experiment/assets/106514085/efe13a33-1f77-45c6-ac9a-8353009a9c46)
