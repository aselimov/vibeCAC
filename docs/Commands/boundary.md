# boundary

```sh
boundary xyz
```

## Inputs

`xyz` - string describing the boundary conditions in which `x`, `y`, and `z` can have values of either `p` for periodic or `s` for shrink-wrapped.

## Examples

```sh
boundary ppp
boundary sps
```

## Description 

Boundary sets the boundary conditions of the model being either periodic or shrink-wrapped. 
Because boundary data is included within the restart file format, this command is not strictly required.
If included this command will override the boundary conditions defined within the data file.
