This program is intended for charge density wave analysis in bulk and low dimensional materials.

Electronic susceptibility (a.k.a. Lindhard function) Chi0(q) from band energies on a full k-mesh is computed and written in Fortran and parallelised over q with MPI.
The real part and the imaginary part (nesting function), both are written at every q-point to one file, `chi_<E_F>_<T>.dat`, along with the contributions of single k-points if requested.

A second program, `chi_bz.x`, puts the result onto the Brillouin zone for overlay plots and comparison with experiments. 
By default it reads Quantum ESPRESSO (QE) output directly, so no manual pre-processing is needed, and writes the intermediate files so the data it used can be checked.

!!We are using constant matrix approximation wherein the el-ph matrix elements are assumed to be unitary!!

## Layout and build

```
src/lindhard_chi.f90     the calculation (MPI; modules and program) ->  bin/lindhard_chi.x
src/chi_bz.f90           chi0 on the Brillouin zone (serial)        ->  bin/chi_bz.x
examples/                chi.in and bz.in examples, plot_bz.py (2D maps and 1D cross sections as PDF)
```

```bash
make                         # mpif90 -O3  ->  bin/lindhard_chi.x, bin/chi_bz.x  (module files in build/)
make clean                   # remove bin/ and build/
```

Requires a Fortran 2008 compiler and MPI (`use mpi`); no other dependencies. `chi_bz.x` makes no MPI calls and runs without `mpirun`.

## Workflow with QE

1. Run a `scf` calculation as usual
2. Then run a, `nscf` calculation on the full unshifted k-mesh, with the k-point list and energies printed in output with verbosity set to `high`, `nosym = .true.` and `noinv = .true.`
3. For more than 40000 k-points raise `npk` parameter in `Modules/parameters.f90` of QE before compiling pw.x
4. Use bands.x on the 'nscf' run with `lsym = .false.` to generate `bandst.dat.gnu` file
5. Change `max_lines` in `PP/src/bands.f90` to match the number of k-points `npk`
6. Run `lindhard_chi.x` with a input file `chi.in` (shown in examples folder) as follows:

# mpirun -np 16 bin/lindhard_chi.x chi.in       
or
# mpirun -np 16 bin/lindhard_chi.x < chi.in

7. The result is written to `chi_<E_F>_<T>.dat` (e.g. `chi_6.4521_300.dat`): `q1 q2 q3  Re_Chi  Im_Chi`
8. Runs with another `temp` or `ef` write new files instead of overwriting
9. Use `chi_bz.x` for the Brillouin-zone overlay then use the supplied python script to plot and generate pdf file as follows:

# chi_bz.x bz.in
# python3 plot_bz.py chi_6.4521_300_ab
or
# python3 plot_bz.py (looks for relevant files in the current directory)

# use `python3 plot_bz.py --help` to understand plotting options


## Convergence

Chi(q) on a finite mesh is a sum of sharply peaked terms smeared only by `k_BT`. Below the temperature at which `k_BT` is comparable to the energy change between neighbouring k-points `E_F`, the result is dominated by sampling noise; the same holds for `sigma` in the imaginary part. 

Check that the features you interpret are stable with respect to both `T` and the mesh density, and whether bands just below or above `E_F` (printed at the start of the run) should be included.
