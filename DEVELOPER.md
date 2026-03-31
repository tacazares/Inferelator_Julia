# Developer Notes — InferelatorJL

Personal reference for working on this codebase.

---

## Starting Julia

```bash
julia          # from any terminal
```

Julia has four REPL modes. You switch between them with single keystrokes:

| Mode | Prompt | Enter with | Exit with |
|---|---|---|---|
| Normal | `julia>` | default | — |
| Package (Pkg) | `(@v1.X) pkg>` | `]` | Backspace |
| Shell | `shell>` | `;` | Backspace |
| Help | `help?>` | `?` | Backspace |

The `@v1.X` in the Pkg prompt shows your Julia version (e.g. `@v1.12` locally, `@v1.7` on the cluster).

---

## One-time setup — local Mac (Julia 1.12)

```julia
julia> ]
(@v1.12) pkg> dev /Users/owop7y/Desktop/InferelatorJL
(@v1.12) pkg> instantiate
```

Then rebuild PyCall so PyPlot works (required once per Julia version):

```julia
(@v1.12) pkg> build PyCall
```

- `dev` registers the package by path — Julia reads directly from your folder, no copying.
- `instantiate` resolves and downloads all dependencies, generating a local `Manifest.toml`.
- `build PyCall` compiles the Python–Julia bridge for your current Python and Julia version.

---

## One-time setup — cluster (Julia 1.7.3)

The cluster does not use `dev` mode — you work directly inside the package folder.
Run this once after cloning or pulling the repo:

```bash
# On the cluster, inside the InferelatorJL directory:
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.build("PyCall")'
```

This generates a fresh `Manifest.toml` on the cluster, resolved specifically for Julia 1.7.3.
The cluster Manifest is gitignored — it never conflicts with your local one.

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
```bash
julia --project=. examples/run_pipeline.jl
```

---

## Running tests

Run the test suite directly — works on all Julia versions:

```bash
# From the terminal, inside the InferelatorJL directory:
julia --project=. test/runtests.jl
```

> **Why not `] test`?**
> In Julia ≥ 1.10, `Pkg.test()` creates a sandbox environment and calls
> `check_registered()` on all test dependencies. Stdlib packages like `Test`
> are not in the General registry, so this fails with
> `"expected package Test to be registered"`.
> Running the test file directly with `julia --project=.` bypasses the sandbox
> and works correctly on all Julia versions.

---

## Working across machines and Julia versions

`Manifest.toml` is **gitignored** — it is machine and Julia-version specific.
Only `Project.toml` is committed. Each machine generates its own Manifest.

| Machine | What to do once | What gets committed |
|---|---|---|
| Local Mac (1.12) | `] instantiate` + `] build PyCall` | nothing (Manifest is gitignored) |
| Cluster (1.7.3) | `Pkg.instantiate()` + `Pkg.build("PyCall")` | nothing |
| Any new machine | same as above | nothing |

**After `git pull` that changes `Project.toml`** (e.g. new dependency added), re-run:
```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

**After upgrading Julia** to a new version, re-run both:
```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.build("PyCall")'
```

**If `Pkg.instantiate()` fails** with `"empty intersection"` or `"unsatisfiable requirements"`:
a specific package's compat bound in `Project.toml` is too narrow for the Julia version
being used. The error names the package. Widen its bound in `Project.toml` and re-run.
Example: `GLMNet = "0.4, 0.5, 0.6, 0.7"` → `"0.4, 0.5, 0.6, 0.7, 0.8"`.

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
] activate                     # return to your default environment (@v1.X)
```

When you `] dev .` from inside a project folder, the package is only registered
in that project, not globally.

---

## Common errors and what they mean

| Error | Cause | Fix |
|---|---|---|
| `UndefVarError: InferelatorJL not defined` | Package not loaded | `using InferelatorJL` |
| `UndefVarError: Revise not defined` | Revise not installed | `] add Revise` |
| `Cannot redefine struct` | Changed `Types.jl` | Restart Julia |
| `MethodError: no method matching ...` | Wrong argument types or order | Check function signature with `?functionname` |
| `KeyError` on a dict field | Field name wrong | `fieldnames(StructType)` to check |
| `PyCall not properly installed` | PyCall not built for this Julia version | `] build PyCall` |
| `expected package Test to be registered` | `] test` sandbox bug in Julia ≥ 1.10 | Use `julia --project=. test/runtests.jl` instead |
| `empty intersection between X@Y.Z and project compatibility` | Compat bound too narrow for installed version | Widen bound for that package in `Project.toml`, then `] instantiate` |
| `Warning: Package X does not have Y in its dependencies` | Missing dep in `Project.toml` | Add it with `] add Y` |

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

test/
  runtests.jl           ← unit tests (run with: julia --project=. test/runtests.jl)
```

---

## Adding a new exported function

1. Write the function in the appropriate `src/` file.
2. Add `export myFunction` to the export block in `src/InferelatorJL.jl`.
3. Add a docstring above the function definition (Julia uses `"""..."""`).
4. If it is a utility function without real-data requirements, add an example
   to `examples/utilityExamples.jl`.

---

## Adding a new dependency

1. `] add PackageName` — this updates both `Project.toml` and `Manifest.toml`.
2. Add `using PackageName` in the appropriate `src/` file.
3. Add a `[compat]` bound in `Project.toml` for the new package.
4. **Commit only `Project.toml`** — `Manifest.toml` is gitignored and should not be committed.
   Other machines will regenerate their own Manifest via `Pkg.instantiate()`.
