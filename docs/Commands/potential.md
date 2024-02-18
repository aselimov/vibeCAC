# potential

```sh
potential pot_type file elements ...
```

## Inputs

`pot_type` - Potential type, currently accepts either "eam" for setfl formatted eam/alloy potentials, or "fs" for eam finnis sinclair style potentials.
`file` - Potential file, currently only accepts .eam.alloy style potential files

`elements ...` - List of atomic elements mapping the atom type number to atomic species. Can be excluded resulting in mapping atomic types based on ordering within the potential file.

## Examples

```sh
potential CuNi.eam.alloy Ni Cu
potential Cu_mishin1.eam.alloy 
```

## Description

This command specifies the potential information.
**Currently only .eam.alloy style potentials are accepted.** 
The `elements ...` portion of the code is used to map the numeric atom type to the potential file.
For the first example, `Ni Cu` maps atoms with type 1 to Ni and atoms with type 2 to Cu.
Inputting `Cu Ni` would instead reverse the order and map atoms with type 1 to Cu and vice-versa.
This can be excluded, as in the second example.
If the CuNi.eam.alloy file is ordered with Ni first and Cu second, then `potential CuNi.eam.alloy Ni Cu` is equivalent to `potential CuNi.eam.alloy`.
