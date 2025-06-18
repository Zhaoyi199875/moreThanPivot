# More Than Pivot for Maximal Clique Enumeration

This project provides implementations of four maximal clique enumeration algorithms, each integrated with our proposed moreThanPivot pruning algorithm:

- BK_degenMTP – moreThanPivot-enhanced variant of BK_degen.
- RMCE_degenMTP – moreThanPivot-enhanced variant of RMCE_degen.
- BK_rcdMTP – moreThanPivot-enhanced variant of BK_rcd.
- HBBMCMTP – moreThanPivot-enhanced variant of HBBMC.

## Compilation

Use the `make` command to compile the project. This will generate four executables:

- `degenmtp` – corresponds to BK_degenMTP
- `rmcemtp` – corresponds to RMCE_degenMTP
- `rcdmtp` – corresponds to BK_rcdMTP
- `hbbmcmtp` – corresponds to HBBMCMTP

```bash
make
```

## Usage

To run an algorithm on a dataset, use the following command format:

```bash
./degenmtp ./dataset
./rmcemtp ./dataset
./rcdmtp ./dataset
./hbbmcmtp ./dataset
```

## Dataset Format

The input graph file is a plain text file with the following format:
```bash
number of vertices
number of directed_edges
v1,v2
v2,v1
...
