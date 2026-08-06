# Method and cleanup notes

- Promoter-associated states E1–E4 use midpoint overlap; E5–E9 use a configurable fractional overlap with historical default 0.5.
- sfCRE means overlap with the same sheep state; sdCRE means overlap with another E1–E9 sheep state; soCRE has no qualifying sheep-state overlap; ssCRE is represented by unmapped human intervals.
- Human TSR assignment is base-pair weighted. Ties are resolved in the order `all_common > broad > intermediate > ts_*`.
- Reuse of the sheep projection is guarded by bytewise equality of all paired E5 sfCRE files unless explicitly disabled.
- The final LDSC BED set contains 68 files under the historical tissue/class exclusions encoded in `08_collect_gwas_beds.sh`.
- LDSC is an external dependency. Version and reference checksums must be recorded by the executing user; they were not recoverable as a complete public environment lockfile.
- A public release still needs the authors' chosen software license and citation metadata.
