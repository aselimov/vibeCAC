# neighbor

```sh
neighbor skin
```

## Inputs

`skin` - The skin distance in angstrom added to the cutoff radius for construction of verlet lists.

## Examples

```sh
neighbor 1.5
```

## Description

This command specifies the additional skin distance which is used when constructing the verlet neighbor lists.
When an atom has moved greater than the half of the skin distance, then the neighbor lists are reconstructed.

