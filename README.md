# Maximal Clique Enumeration with moreThanPivot

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
