using Arrhenius
using PyCall

```
1. config Conda.jl
2. add PyCall
3. build PyCall
```
ct = pyimport("cantera")

gas = ct.Solution("gri30")

gas.X
