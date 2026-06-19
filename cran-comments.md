# cran-comments for GLOWr 0.1.0

## Test environment
- R 4.3.3 on x86_64-conda-linux-gnu
- check command: `R CMD check --as-cran --no-manual`
- NOTE: no LaTeX (pdflatex) was available in this environment, so the PDF manual was not built (`--no-manual`). The Rd files themselves pass all Rd checks; the PDF manual builds on a LaTeX-equipped machine (e.g. CRAN).

## R CMD check results
- 0 ERROR | 0 WARNING | 5 NOTE

- No ERRORs.

- No WARNINGs.

### NOTEs (verbatim from the check log)

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Zheyang Wu <zheyangwu@wpi.edu>’

Unknown, possibly misspelled, fields in DESCRIPTION:
  ‘Remotes’
```

```
* checking for future file timestamps ... NOTE
unable to verify current time
```

```
* checking compilation flags used ... NOTE
Compilation used the following non-portable flag(s):
  ‘-march=nocona’
```

```
* checking compiled code ... NOTE
Warning in read_symbols_from_object_file(so) :
  this requires 'nm' to be on the PATH
Warning in read_symbols_from_object_file(so) :
  this requires 'nm' to be on the PATH
Warning in read_symbols_from_object_file(so) :
  this requires 'nm' to be on the PATH
File ‘GLOWr/libs/GLOWr.so’:
  Found no calls to: ‘R_registerRoutines’, ‘R_useDynamicSymbols’

It is good practice to register native routines and to disable symbol
search.

See ‘Writing portable packages’ in the ‘Writing R Extensions’ manual.
```
