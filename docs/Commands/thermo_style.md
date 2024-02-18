# thermo_style

```sh
thermo_style args
```

## Inputs

`args` - List of thermo_style options. Acceptable options are listed below:

## Examples

```
thermo 50
```

## Descriptions

The `thermo_style` command sets the thermo output for CAC. 
By default the thermostyle command only outputs the potential energy and the global force norm.
Acceptable options for thermo_style are listed below:

- `pe`: potential energy
- `fnorm`: global force norm
- `temp`: temperature
- `temp_c`: temperature in CG 
- `temp_a`: temperature in atoms
- `ke`: kinetic energy
- `ke_c`: kinetic energy in CG 
- `ke_a`: kinetic energy in atoms
- `lx`, `ly`, `lz`: Box dimensions in x, y, z respectively
- `px`, `py`, `pz`: Box pressure in x, y, z respectively
