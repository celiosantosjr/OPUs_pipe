# summarized_OPUs_annotation.tsv.xz

It is a tab separated xzipped table containing the annotation
for the generated OPUs (Operational Protein Units). This annotation
follows the KO terms and is given as the KO term and the percent of
proteins in the OPU assigned to it. Proteins without annotation
in the emapper2 output given will be assigned as UNKNOWN and their
percent in the group will be given always when present.

The columns are:

| **name** | **description** |
| :---: | :---: |
| OPU | Operational protein unit created to accommodate the protein |
| annotation | KO term assigned to this protein or 'UNKNOWN' (percent of proteins in the group%) | 

**Example:**

```
OPU	annotation
OPU0	UNKNOWN (100.00%)
OPU1	ko:K03043 (100.00%)
OPU10	ko:K01601,ko:K01963 (100.00%)
OPU100	ko:K00097,ko:K22024 (100.00%)
OPU1000	ko:K00128 (100.00%)
OPU10000	UNKNOWN (100.00%)
OPU100000	ko:K06997 (100.00%)
OPU1000000	ko:K03723,ko:K05365 (100.00%)
OPU1000001	ko:K01256 (100.00%)
```

---

[Click here to go back](../README.md)
