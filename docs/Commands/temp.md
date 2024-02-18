# temp

```sh
temp command group_id args
```

## Inputs

`command`  - Command for temp option, either `create` or `control`.

`group_id` - Group id to apply temperature control to 

**Additional arguments depend on the command passed and are listed below** 

```sh
temp create group_id T_target
```

- `T_target` - Target temperature

```sh
temp control group_id T_target time_constant
```

- `T_target` - Target temperature 
- `time_constant` - time constant 
