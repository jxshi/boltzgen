# GLP-1R 6X18 reference audit

This is an analysis of observed reference geometry, not a Boltz API result.

The input is the user-provided 6X18.cif. Receptor author chain R corresponds to label chain F; native GLP1 author P corresponds to label E. The design input contains receptor R only, with no native peptide, auxiliary proteins, or waters.

Canonical target: human GLP1R P43220 residues 29–423, 395 sequence positions and 384 residues with coordinates. Missing 130–135 (KRGERS) and 339–343 (LMCKT) remain in the sequence and are not concatenated away. F260 in 6X18 is restored to canonical L260 by retaining N/CA/C/O/CB only; the Leu side chain is not modeled. Original retained coordinates are unchanged. Native disulfides are preserved.

Numbering: original label_seq_id = author + 9; prepared label_seq_id = author − 28; API zero-based index = author − 29. These API indices are valid only with the complete prepared sequence and no crop/reindexing.

The prior canonical CIF passed an independent, atom-by-atom check against the original mmCIF. The complete mapping includes all 395 positions, including unobserved residues. The prior residue_map.csv used a coordinate-only sequence_index that skipped missing residues; that column must not be used as an API index.

| Site | Original label F | Prepared label R | API index | Nearest native GLP1 | Heavy-atom distance Å | Atom pair |
|---|---:|---:|---:|---|---:|---|
| R190 | 199 | 162 | 161 | E9 | 2.7785 | NH2–OE2 |
| Q234 | 243 | 206 | 205 | H7 | 3.1159 | OE1–NE2 |
| T298 | 307 | 270 | 269 | S17 | 2.4623 | OG1–OG |
| R299 | 308 | 271 | 270 | S17 | 2.6687 | NH1–OG |
| N300 | 309 | 272 | 271 | S14 | 2.7246 | N–OG |
| W306 | 315 | 278 | 277 | H7 | 3.6124 | CH2–NE2 |
| R310 | 319 | 282 | 281 | H7 | 3.3274 | CG–N |
| R380 | 389 | 352 | 351 | D15 | 2.6209 | NH2–OD2 |
| L384 | 393 | 356 | 355 | T11 | 3.7870 | CD2–CG2 |
| L388 | 397 | 360 | 359 | F12 | 3.6925 | CD1–CG |
| L32 | 41 | 4 | 3 | E21 | 2.9608 | N–OE2 |
| W39 | 48 | 11 | 10 | R36 | 3.4207 | CH2–NE |
| Y69 | 78 | 41 | 40 | V33 | 3.5496 | CD1–CG2 |
| W91 | 100 | 63 | 62 | I29 | 3.1345 | CZ2–CG2 |

A sparse first design hypothesis can guide Q234/W306/R380 (API 205/277/351) and assess the other sites afterward. This does not impose GLP1-like peptide orientation or N-terminal activation geometry. ECD anchoring and N-terminal pocket engagement require explicit pose inspection after prediction.

Distances are minimum Euclidean distances between observed heavy atoms in model 1. Full residue-pair contacts at ≤4.5 Å are in native_contacts_4p5A.tsv. Missing atoms, missing loops, and absent glycans limit the interpretation. Structure confidence or contacts alone cannot establish agonism.

The reference GLP1 sequence has 30 residues (author 7–36), HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR. No terminal chemical modification is added or inferred from the plain sequence.

Reproduce with Python 3.9+ (standard library only):

```sh
python3 scripts/prepare_reference.py
```

Sources: [6X18](https://www.rcsb.org/structure/6X18), [coordinate archive](https://files.rcsb.org/download/6X18.cif), [UniProt P43220](https://www.uniprot.org/uniprotkb/P43220/entry). Cached source hashes and exact preparation choices are recorded in structure_audit.json.
