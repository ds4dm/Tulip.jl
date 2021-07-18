# TulipCL

App for building a command-line executable of the [Tulip.jl](https://github.com/ds4dm/Tulip.jl) interior-point solver.
All commands shown below are executed in the `Tulip.jl/app/` directory.

## Installation instructions

1. Download and install Julia (version 1.3.1 or newer).

1. Install `PackageCompiler`
    ```bash
    $ julia -e 'using Pkg; Pkg.add("PackageCompiler")'
    ```

1. Instantiate the current environment
    ```bash
    $ julia --project=. -e 'using Pkg; Pkg.instantiate()'

1. Build the command-line executable
    ```julia
    $ julia -q --project=.
    julia> using PackageCompiler
    julia> create_app(".", "tulip_cl", force=true, precompile_execution_file="precompile_app.jl", app_name="tulip_cl");
    julia> exit()
    ```
    The executable will be located at `tulip_cl/bin/tulip_cl`.

    While building the command-line executable, a precompilation script is executed wherein a small number of problems will be solved, thereby producing a number of logs. This process is normal.
    For more information on how to build the app, take a look at the [PackageCompiler documentation](https://julialang.github.io/PackageCompiler.jl/dev/apps/).

### Using a different version of Tulip

Following the completion of step 3. above, this directory will contain a `Manifest.toml` that specifies the version of all packages, including that of `Tulip`.
By default, this will be the most recently tagged version.

To build the executable with a different version/branch of Tulip follow these instructions, tag the particular version/branch before creating the app.

For instance, to try out a local development branch, running
```julia
julia> ]
pkg> dev ..
```
will use the current state of the Tulip.jl repository.

Finally, follow the last step in the above installation guide to generate the command-line executable.

## Running the command-line executable

Once the build step is performed, the executable can be called from the command line as follows:
```bash
tulip_cl [options] finst
```
where `finst` is the problem file. For instance,
```bash
tulip_cl --Threads 1 --TimeLimit 3600 afiro.mps
```
will solve the problem `afiro.mps` using one thread and up to 1 hour of computing time.

Currently, possible user options are

| Option name | Type | Default | Description |
|-------------|------|---------|-------------|
| `Presolve`  | `Int`     | `1`   | Set to `0` to disable presolve, `1` to activate it |
| `Threads`   | `Int`     | `1`   | Maximum number of threads |
| `TimeLimit` | `Float64` | `Inf` | Time limit, in seconds |
| `IterationsLimit` | `Int` | `500` | Maximum number of barrier iterations |
| `Method` | `String` | `HSD` | Interior-point method |

For more information, run `tulip_cl --help`, or look at Tulip's [documentation](https://ds4dm.github.io/Tulip.jl/stable/) for more details on parameters.