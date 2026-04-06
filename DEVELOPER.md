# Contributing to InferelatorJL

Guidelines for developers working on this codebase. Covers environment setup,
development workflow, adding new features, and common pitfalls.

---

## Getting started with Julia

```bash
julia          # launch from any terminal
```

Julia has four REPL modes. Switch between them with single keystrokes:

| Mode | Prompt | Enter with | Exit with |
|---|---|---|---|
| Normal | `julia>` | default | — |
| Package (Pkg) | `(@v1.X) pkg>` | `]` | Backspace |
| Shell | `shell>` | `;` | Backspace |
| Help | `help?>` | `?` | Backspace |

The `@v1.X` in the Pkg prompt reflects your Julia version (e.g. `@v1.12` locally, `@v1.7` on the cluster).

---

## One-time setup (per machine / Julia version)

**Interactive session (REPL available):**

```
julia> ]
(@v1.X) pkg> activate /path/to/InferelatorJL   ← step 1: enter package environment
(InferelatorJL) pkg> instantiate                ← step 2: install all dependencies
(InferelatorJL) pkg> build PyCall               ← step 3: build Python bridge
(InferelatorJL) pkg> activate                   ← step 4: return to global (no argument)
(@v1.X) pkg> dev /path/to/InferelatorJL         ← step 5: register package globally
```

> **Why this order matters:**
> `dev` must be run from the global environment (`@v1.X`), not from inside the
> package environment. Running `dev .` while InferelatorJL is active will error with
> `"has the same name or UUID as the active project"`.
> Always return to global with `activate` (no argument) before running `dev`.

**Non-interactive (no REPL, e.g. a batch job):**

```bash
julia --project=/path/to/InferelatorJL -e 'using Pkg; Pkg.instantiate()'
julia --project=/path/to/InferelatorJL -e 'using Pkg; Pkg.build("PyCall")'
```

> `dev` is not needed in non-interactive mode — scripts are run with
> `julia --project=/path/to/InferelatorJL` which activates the environment automatically.

- `instantiate` resolves all dependencies from `Project.toml` and generates a local `Manifest.toml`.
- `build PyCall` compiles the Python–Julia bridge — required once per Julia version.
- `dev` registers the package globally so `using InferelatorJL` works from any session without activation.
- `Manifest.toml` is gitignored — each machine generates its own.

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

**Any new machine — do once:**

```
# Interactive
] activate /path/to/InferelatorJL
] instantiate
] build PyCall

# Non-interactive
julia --project=/path/to/InferelatorJL -e 'using Pkg; Pkg.instantiate(); Pkg.build("PyCall")'
```

**After `git pull` that changes `Project.toml`** (new dependency added):

```
# Interactive              |  Non-interactive
] activate /path/...       |  julia --project=. -e 'using Pkg; Pkg.instantiate()'
] instantiate              |
```

**After upgrading Julia** to a new version — re-run both instantiate and build PyCall (same commands as above).

**If `Pkg.instantiate()` fails** with `"empty intersection"` or `"unsatisfiable requirements"`:
a specific package's compat bound in `Project.toml` is too narrow for the Julia version
being used. The error names the package — widen its bound in `Project.toml` and re-run.
Example: `GLMNet = "0.4, 0.5, 0.6, 0.7"` → `"0.4, 0.5, 0.6, 0.7, 0.8"`.

---

## Inspecting the loaded package

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
# Currently using dev (local source folder):
] status InferelatorJL     # shows path to local folder

# Switch to the released version (once published):
] free InferelatorJL       # removes the dev pin
] add InferelatorJL        # installs from registry

# Switch back to dev:
] dev /path/to/InferelatorJL
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

## Common errors and fixes

| Error | Cause | Fix |
|---|---|---|
| `UndefVarError: InferelatorJL not defined` | Package not loaded | `using InferelatorJL` |
| `UndefVarError: Revise not defined` | Revise not installed | `] add Revise` |
| `Cannot redefine struct` | Changed `Types.jl` | Restart Julia |
| `MethodError: no method matching ...` | Wrong argument types or order | Check signature with `?functionname` |
| `KeyError` on a dict field | Field name wrong | `fieldnames(StructType)` to check |
| `PyCall not properly installed` | PyCall not built for this Julia version | `] build PyCall` |
| `has the same name or UUID as the active project` | Ran `] dev .` inside the package environment | `] activate` first, then `] dev /path/to/InferelatorJL` |
| `expected package Test to be registered` | `] test` sandbox bug in Julia ≥ 1.10 | Use `julia --project=. test/runtests.jl` instead |
| `empty intersection between X@Y.Z and project compatibility` | Compat bound too narrow | Widen bound in `Project.toml`, then `] instantiate` |
| `Warning: Package X does not have Y in its dependencies` | Missing dep in `Project.toml` | Add it with `] add Y` |

---

## Package structure

```
src/
  InferelatorJL.jl      ← module entry point, all includes and using statements
  Types.jl              ← ALL struct definitions (edit here for new fields)
  API.jl                ← public API (loadData, buildNetwork, evaluateNetwork, etc.)
  data/                 ← data loading functions
  prior/                ← TF merging
  utils/                ← DataUtils, NetworkIO, PartialCorrelation
  grn/                  ← pipeline core (PrepareGRN, BuildGRN, AggregateNetworks,
                           RefineTFA, PlotInstability)
  metrics/              ← PR/ROC evaluation and plotting

examples/
  interactive_pipeline.jl      ← step-by-step, public API
  interactive_pipeline_dev.jl  ← step-by-step, internal calls
  run_pipeline.jl              ← function-wrapped, public API
  run_pipeline_dev.jl          ← function-wrapped, internal calls
  utilityExamples.jl           ← utility function demos, no real data needed
  plotPR.jl                    ← standalone PR curve evaluation and plotting

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
5. Add a test to `test/runtests.jl`.

---

## Adding a new internal (non-exported) function

1. Write the function in the appropriate `src/` file.
2. Add its name to the internal comment block at the bottom of `src/InferelatorJL.jl`
   so it is discoverable.
3. It remains accessible via `import InferelatorJL: myFunction` for power users.

---

## Adding a new dependency

1. `] add PackageName` — updates both `Project.toml` and `Manifest.toml`.
2. Add `using PackageName` in the appropriate `src/` file.
3. Add a `[compat]` bound in `Project.toml` for the new package.
4. **Commit only `Project.toml`** — `Manifest.toml` is gitignored.
   Other machines regenerate their own Manifest via `Pkg.instantiate()`.
