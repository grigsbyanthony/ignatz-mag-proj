<p align="center">
  <img src="meta/ignatz_bacteria.png" alt="Bacteria" width="100"/>
</p>

<p align="center">
  <img src="meta/ignatzschineria-banner.gif" alt="Ignatzschineria banner" width="600"/>
</p>

# File directory
```

```

# Ignatzschineria MAG Analysis Pipeline

## PHASE 1: QUALITY ASSESSMENT
```
checkm2 predict --threads 8 --input *.fna --output-directory results/quality_assessment/checkm2_results

for file in *.fna; do
    base=$(basename "$file" .fna)
    echo "Processing $base..."
    seqkit stats -a "$file" > "results/basic_stats/${base}_detailed_stats.txt"
done

for file in *.fna; do
    seqkit stats -T "$file" | tail -n 1 >> results/basic_stats/genome_summary.csv
done
```

## PHASE 2: GENE PREDICTION
```
for file in *.fna; do
    base=$(basename "$file" .fna)
    echo "Predicting genes for $base..."
    
    prodigal -i "$file" \
             -o "results/gene_prediction/${base}_genes.gff" \
             -a "results/gene_prediction/${base}_proteins.faa" \
             -d "results/gene_prediction/${base}_nucleotides.fna" \
             -f gff -p meta
done
```

## PHASE 3: FUNCTIONAL ANNOTATION
```
download_eggnog_data.py -y --data_dir databases/eggnog

for proteins in results/gene_prediction/*_proteins.faa; do
        base=$(basename "$proteins" _proteins.faa)
        echo "Annotating $base..."
        emapper.py -i "$proteins" \
                   --output "${base}_eggnog" \
                   --output_dir results/functional_annotation/eggnog \
                   --data_dir databases/eggnog \
                   --cpu 8 --override
done
```

## PHASE 4: DECOMPOSITION ANALYSIS
```

```

## PHASE 5: SECONDARY METABOLITE ANALYSIS
```

```
