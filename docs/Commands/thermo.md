# thermo

```sh
thermo step
```

## Inputs

`step` - Output every `step` steps.

## Examples

```
thermo 50
```

## Descriptions

Currently the thermo output is fixed. 
For energy minimization, the thermo output contains only the energy and second norm of the global force vector.
For dynamics, the thermo output also contains the kinetic energy and the temperature.

