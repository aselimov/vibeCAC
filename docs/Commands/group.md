# group

```sh
group id selection shape ...
```

## Inputs 

`id` - assigned id for the group.

`selection` - either `atoms`, `elements`, or `all`. 

`shape` - either `block`, `sphere`, or `type`.

**Additional arguments depend on the group shape passed and are listed below:** 

```sh
group id selection block xlo xhi ylo yhi zlo zhi
```
where `{x,y,z}lo` are the lower bounds along `x,y,z` and `{x,y,z}hi` are the upper bounds along `x,y,z`.

```sh
group id selection sphere x y z r
```
where `x y z` define the centroid of the sphere and `r` defines the radius of the sphere.

```sh
group id selection type n
```
where `n` is a number which selects the `n`th atom type.

## Examples
```sh
group tophalf all block -inf inf -inf inf inf*0.5 inf
group precip atoms sphere 100 78 78 10
group Ni atoms type 2
```

## Description

The group option is used in conjunction with other commands which require a group input. 
**A default `all` group is defined which contains all atoms and elements.**
A subsection of all atoms or elements can be selected using the group command.
`id` is used to refer to the group in other commands.
`selection` determines which element types are selected within the group. 
A group can be defined containing only `atoms`, `elements`, or `both`. 
If element types are within the bounds of the group but are not specified with `selection` they are not defined as members of the group.

All positions for the group command can be specified using position specifications described [here](../Misc/position.md).



