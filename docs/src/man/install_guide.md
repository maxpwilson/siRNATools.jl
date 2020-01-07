# Installation Guide

First download and install Julia from [here](https://julialang.org/downloads/) if not already installed

Next start Julia and type ] into the command line.  The screen should look like below.

```julia
(v1.3) Pkg> 

```

Next add the dependencies as shown below.
```julia
(v1.3) Pkg> add CSV, DataFrames, StatsBase, StringDistances, GZip, ProgressMeter, BSON, Flux, Query, PyCall
```

Finally install the package.

```julia
(v1.3) Pkg> add https://github.com/maxpwilson/siRNATools.jl.git

```