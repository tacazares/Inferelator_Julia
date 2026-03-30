# Developer Notes — InferelatorJL

Personal reference for working on this codebase.

---

## Starting Julia

```bash
julia          # from any terminal
```

Julia has three REPL modes. You switch between them with single keystrokes:

| Mode | Prompt | Enter with | Exit with |
|---|---|---|---|
| Normal | `julia>` | default | — |
| Package (Pkg) | `(@v1.7) pkg>` | `]` | Backspace |
| Shell | `shell>` | `;` | Backspace |
| Help | `help?>` | `?` | Backspace |

---

## One-time setup (do once per machine)

```julia
julia> ]
(@v1.7) pkg> dev /Users/owop7y/Desktop/InferelatorJL
(@v1.7) pkg> instantiate
```

- `dev` registers the package by path — Julia reads directly from your folder, no copying.
- `instantiate` downloads and installs all dependencies listed in `Project.toml`.

You only need to run these **once**. After that the package is permanently available in
your default environment.

---

## Every session

```julia
using Revise            # must come BEFORE the package
using InferelatorJL     # loads the package
```

**Always load Revise first.** It monitors `src/` for changes and patches the running
session when you save a file — no restart needed.

---

## Revise: what it can and cannot do

| Change | Revise handles it? |
|---|---|
| Edit a function body | ✅ Live — takes effect on next call |
| Add a new function | ✅ Live |
| Add a new `include(...)` to `InferelatorJL.jl` | ✅ Live |
| Add or rename a field in a struct (`Types.jl`) | ❌ Must restart Julia |
| Change a struct's field type | ❌ Must restart Julia |
| Rename or delete a struct | ❌ Must restart Julia |

Struct changes are the only thing that forces a restart. Everything else is live.

---

## Running examples

From the REPL after loading the package:
```julia
include("examples/interactive_pipeline.jl")
include("examples/utilityExamples.jl")
```

From the terminal (no REPL):
```julia
julia examples/run_pipeline.jl
```

---

## Checking what changed / what's loaded

```julia
# Which version is active?
using Pkg; Pkg.status("InferelatorJL")

# Where is it loaded from?
pathof(InferelatorJL)

# What does the package export?
names(InferelatorJL)

# What fields does a struct have?
fieldnames(GeneExpressionData)
fieldnames(GrnData)

# Inspect a loaded struct at runtime
data = GeneExpressionData()
propertynames(data)
```

---

## Switching between dev and release versions

```julia
# Currently using dev (your local folder):
] status InferelatorJL     # shows   InferelatorJL [path] /Users/owop7y/Desktop/InferelatorJL

# Switch to the released version (once it is published):
] free InferelatorJL       # removes the dev pin
] add InferelatorJL        # installs from registry

# Switch back to dev:
] dev /Users/owop7y/Desktop/InferelatorJL
```

---

## Project environments (keeping work isolated)

Every directory can have its own `Project.toml`. To work inside a specific project
environment (e.g., a collaborator's analysis folder):

```julia
] activate /path/to/project    # switch to that environment
] status                       # see what is installed there
] activate                     # return to your default environment (@v1.7)
```

When you `] dev .` from inside a project folder, the package is only registered
in that project, not globally.

---

## Running tests

```julia
] test InferelatorJL
```

Or from inside the package environment:
```julia
] activate .
] test
```

---

## Common errors and what they mean

| Error | Cause | Fix |
|---|---|---|
| `UndefVarError: InferelatorJL not defined` | Package not loaded | `using InferelatorJL` |
| `UndefVarError: Revise not defined` | Revise not installed | `] add Revise` |
| `Cannot redefine struct` | Changed `Types.jl` | Restart Julia |
| `MethodError: no method matching ...` | Wrong argument types or order | Check function signature with `?functionname` |
| `KeyError` on a dict field | Field name wrong | `fieldnames(StructType)` to check |
| `Warning: Package X does not have Y in its dependencies` | Missing dep in `Project.toml` | Add it with `] add Y` then add UUID to Project.toml |

---

## Package structure quick reference

```
src/
  InferelatorJL.jl      ← module entry point, all includes and using statements
  Types.jl              ← ALL struct definitions (edit here for new fields)
  API.jl                ← public API (loadData, buildNetwork, etc.)
  data/                 ← data loading functions
  prior/                ← TF merging
  utils/                ← DataUtils, NetworkIO, PartialCorrelation
  grn/                  ← pipeline core (PrepareGRN, BuildGRN, AggregateNetworks, RefineTFA)
  metrics/              ← PR/ROC evaluation and plotting

examples/
  interactive_pipeline.jl      ← step-by-step, public API (your main reference)
  interactive_pipeline_dev.jl  ← step-by-step, internal calls (your dev reference)
  run_pipeline.jl              ← function-wrapped, public API
  run_pipeline_dev.jl          ← function-wrapped, internal calls
  utilityExamples.jl           ← utility function demos, no real data needed
  plotPR.jl                    ← PR curve evaluation
```

---

## Adding a new exported function

1. Write the function in the appropriate `src/` file.
2. Add `export myFunction` to the export block in `src/InferelatorJL.jl`.
3. Add a docstring above the function definition (Julia uses `"""..."""`).
4. If it is a utility function without real-data requirements, add an example
   to `examples/utilityExamples.jl`.

## Adding a new dependency

1. `] add PackageName` (this updates `Project.toml` and `Manifest.toml` automatically).
2. Add `using PackageName` in the appropriate `src/` file.
3. Add a `[compat]` bound in `Project.toml` for the new package.
