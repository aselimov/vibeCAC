# min_style

```sh
min_style minimizer
```

## Inputs

`minimizer` - Either `fire` or `cg`.

## Examples

```sh
min_style fire
min_style cg
```

## Description

This code sets the minimizer for future `minimize` commands.
Currently CAC has implemented both conjugate gradient minimization, `cg`, and the fast inertial relaxation engine [1], `fire`. 

## References

[1] Guénolé, Julien et al. "Assessment and optimization of the fast inertial relaxation engine (fire) for energy minimization in atomistic simulations and its implementation in lammps" Computational Materials Science (2020)
