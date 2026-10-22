# Install: genomic_context

No install needed. Pure Python standard library, no external tool.

The module needs, for each input CDS FASTA, the genome FASTA and the GFF3
annotation the CDS were taken from.

## Files

Put them in one directory (not the input directory), named by the input
file stem:

```
genomes/
  Athaliana.gff3          # or .gff, optionally .gz
  Athaliana.fa            # or .fasta / .fna, optionally .gz
  Athaliana.seqid_map.tsv # optional, see below
```

```
python run_features.py --modules genomic_context \
    --batch cds/ --out-dir outputs/ --module-ref genomic_context=genomes/
```

The script can also be run on its own:

```
python modules/genomic_context/compute_genomic_context.py cds/Athaliana.fna \
    --gff genomes/Athaliana.gff3 --genome genomes/Athaliana.fa --out Athaliana.tsv
```

## Matching CDS to the annotation

Each CDS FASTA ID (first word of the header), or the `[protein_id=...]`
tag of NCBI `cds_from_genomic` headers, is looked up in the GFF as:
- a transcript (`mRNA`/`transcript`) `ID` or `Name`;
- a `CDS` `protein_id`, `Protein_Accession` or `Name`.

Each is also tried with or without a `.p` suffix and an `rna-` or
`transcript:` prefix. If none of the IDs match, the run stops with an
error.

## Sequence names

The GFF and the genome FASTA must name sequences the same way. NCBI
GenBank (GCA) and RefSeq (GCF) copies of one assembly do not: `CM052939.1`
in one is `NC_079511.1` in the other. The run stops with an error when
most matched CDS lie on sequences missing from the genome.

Fix it with `<stem>.seqid_map.tsv` (or `--seqid-map` when running the
script directly): tab-separated, first column the GFF name, second the
genome name; further columns are ignored. NCBI's assembly report lists
both names for every sequence.

## Output columns

| Column | Meaning |
|---|---|
| `cds_exons` | number of CDS segments |
| `cds_span` | genomic span from the first to the last coding base |
| `cds_intron_total` | total intron length between CDS segments |
| `kozak_score` | sum of log10 base frequencies over the 12 nt from −6 to +6 around the start codon; empty if the CDS does not start with ATG |
| `overlaps_gene` | 1 if the CDS span overlaps another annotated gene (either strand) |
| `nearest_gene_dist` | bp to the nearest other gene on the same sequence (0 if overlapping) |
| `ssr_count`, `ssr_bp`, `ssr_frac` | perfect tandem repeats in the CDS (motif 1–6 nt, MISA minimum copy numbers 10/6/5/5/5/5) |

Notes:
- Only CDS coordinates are used, because UTR annotation differs widely
  between genomes.
- The Kozak base frequencies are estimated from all sequences in the
  input file. Scores are therefore comparable within one input, and depend
  on which sequences it contains (for example all isoforms vs one per gene).
- Isoforms of one gene can start at different codons, with different
  start contexts. One isoform per gene (for example the longest) keeps
  genes with many isoforms from dominating the frequencies, but the chosen
  isoform is not necessarily the main translated one. Use the same
  isoform rule for every genome you want to compare.
- `overlaps_gene` and `nearest_gene_dist` are empty when the GFF has no
  `gene` features on that sequence.
- Sequences that cannot be placed keep a row with empty feature columns.
  The stderr summary reports how many were found in the GFF, placed on the
  genome, and given a Kozak window.
