import pandas as pd
import sys

# sys.argv[1]: results from nominal eQTL pass
# sys.argv[2]: nominal p-value thresholds for each eGene
# sys.argv[3]: outfile

nom_thresholds = pd.read_csv(sys.argv[2])

with open(sys.argv[3], "a+") as output:
  # Write header
    output.write('gene_id gene_chr gene_start gene_end gene_strand num_cis_variants dist variantID variant_chr variant_start variant_end nom_pval r_squared beta best_hit pval_nominal_threshold\n')
    with open(sys.argv[1], "r") as f:
        for line in f:
            qtldata = line.rstrip().split(" ")
            gene_id = qtldata[0]

            # Find gene_id in nom_thresholds and get corresponding nom pvalue threshold
            pval_threshold = nom_thresholds[nom_thresholds["gene_id"] == gene_id]["pval_nominal_threshold"].values[0]

            # If satisfies threshold, write to file
            nom_pval = float(qtldata[11])
            if nom_pval <= pval_threshold:
              new_line = line.rstrip() + " " + str(pval_threshold)
                output.write(new_line)
