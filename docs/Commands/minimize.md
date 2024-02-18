# minimize 

```sh
minimize ftol etol max_iter
```

## Inputs

`ftol` - force tolerance

`etol` - energy tolerance

`max_iter` - maximum number of iterations to calculate before exiting.

## Examples

```sh
minimize 10d-10 10d-10 10000
```

## Description

The minimize command is used to perform energy minimization on the CAC model. 
The minimization method is chosen by the [min_style](./min_style.md) command. 
Minimization is considered to be complete under the following conditions:

1. \( \frac{(E_i - E_{i-1})}{E_i} < etol \), where \( i \) is the current iteration.
2. No atom has a force component greater than `ftol`
3. The number of iterations is greater than `max_iter`


