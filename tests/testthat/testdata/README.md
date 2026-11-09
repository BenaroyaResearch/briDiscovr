# Example data for briDiscovr tests

Small flow cytometry dataset used by the end-to-end tests (see `../helper-exampleData.R`).
Copied from `data/dataTest/` in the BenaroyaResearch/briDiscovr-2 repository
(created 2025-11-14 by the `briDiscovr_create_test_dataset` script).

## Contents

- `fcs/`: 20 downsampled, compensated FCS 3.0 files from a B-cell panel:
  10 donors (`donor1`–`donor10`) × 2 cell subsets.
  - `Tet-` (tetramer-negative B cells): the **parent population**; 2,000 events per donor.
  - `Tet+` (tetramer-positive B cells): the rare subset; 784–2,000 events per donor.
- `fcsInfo.csv`: one row per FCS file (`subject`, `cellSubset`, `filename`).
  Filenames are bare file names. The tests turn them into absolute paths
  at run time, so this file works on any machine.
- `markerInfo.csv`: marker table (`fixed` = common marker name, `desc` = FCS
  description, `useToCluster`, `normalizationMethod`). 10 clustering markers;
  normalization methods include `zScore`, `none`, `warpSet`, and `warpSet3`.

The original `fcsInfo` file in briDiscovr-2 held absolute paths from the
creator's computer, so it was rewritten here; the FCS and markerInfo contents
are unchanged.
