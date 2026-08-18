                    GEO
                     │
                     ▼
               GSE300475
                     │
          ┌──────────┴──────────┐
          ▼                     ▼
   Supplementary files      Metadata
          │                     │
          ▼                     ▼
       10X data            mapping table
          │
          ▼
      Read10X()
          │
          ▼
   ┌───────────────┐
   │ Gene Expression│
   │ Antibody Capture│
   └───────┬───────┘
           │
           ▼
   CreateSeuratObject
           │
           ▼
     RNA + HTO match
           │
           ▼
       HTO Assay
           │
           ▼
     CLR Normalize
           │
           ▼
       HTODemux
           │
      ┌────┼─────┐
      ▼    ▼     ▼
  Singlet Doublet Negative
      │
      ▼
    QC plots
      │
      ▼
 Keep Singlets
      │
      ▼
 Process all samples
      │
      ▼
      Merge
      │
      ▼
 merged_singlets.rds