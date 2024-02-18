# read_data

```
read_data filename
```

## Inputs 

`filename` - Input file to be read in restart format.

## Examples 

```sh
read_data model_1.restart
```

## Description

Model building capabilities are not included into the CAC simulation tool. 
Instead the [CAC model builder (CACmb)](https://gitlab.com/aselimov/cacmb) must be used to create models in the .restart format.
The read_data command reads in one .restart formatted data file and distributes the model across all processors.



