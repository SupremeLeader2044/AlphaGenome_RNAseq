from alphagenome import colab_utils
from alphagenome.data import gene_annotation, genome, transcript as transcript_utils
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np
import os

#==========AlphaGenome Function==========
def AlphaGenomePredictATAC(api_key, chromosome, position, interval_size, ref, alt, ontology_terms, length_keys):
    dna_model = dna_client.create(api_key)
    interval = genome.Interval(
        chromosome=chromosome,
        start=max(position - interval_size // 2, 1),
        end=position + interval_size // 2
    ).resize(length_keys)

    variant = genome.Variant(chromosome, position, ref, alt)

    ATAC_output = dna_model.predict_variant(
        interval=interval,
        variant=variant,
        ontology_terms=ontology_terms,
        requested_outputs=[dna_client.OutputType.ATAC]
    )

    return ATAC_output

#==========ATAC Visualization==========
def AlphaGenomeATACVisualisation(ATAC_output, interval_size):
    ref_atac = ATAC_output.reference.atac
    alt_atac = ATAC_output.alternate.atac

    # Get start position and total number of values
    interval_start = ref_atac.interval.start
    resolution = ref_atac.resolution
    num_points = ref_atac.values.shape[0]

    positions = range(interval_start, interval_start + num_points * resolution, resolution)
    ref_values = ref_atac.values.flatten()
    alt_values = alt_atac.values.flatten()

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(positions, ref_values, label='Reference', color='steelblue', linewidth=1)
    ax.plot(positions, alt_values, label='Alternate', color='tomato', linewidth=1)

    # Optional: mark variant at center if applicable
    variant_pos = interval_start + (num_points // 2)
    ax.axvline(x=variant_pos, color='black', linestyle='--', label='Variant')

    biosample = ref_atac.metadata.iloc[0]["biosample_name"]
    ontology_id = ref_atac.metadata.iloc[0]["ontology_curie"]
    ax.set_title(f"ATAC-seq: {biosample} ({ontology_id})")
    ax.set_xlabel("bp")
    ax.set_ylabel("Chromatin Accessibility")
    ax.legend()
    return fig


#==========Main Driver Function==========
def ExportATACToPDF():
    api_key = "Your API Key"
    interval_size = 1048576
    ontology_terms = ["UBERON:0002078"]
    length_keys = dna_client.SEQUENCE_LENGTH_1MB

    input_path = input("Enter Input Excel file path: ").replace("\\", "/").replace("\"", "")
    variant_dataset = pd.read_excel(input_path)
    input_folder = os.path.dirname(input_path)
    output_pdf_path = os.path.join(input_folder, "ATAC_Report.pdf")

    with PdfPages(output_pdf_path) as pdf:
        for i in range(len(variant_dataset)):
            unique_varID = variant_dataset.iloc[i, 3]
            chromosome_number = variant_dataset.iloc[i, 11]
            variant_position = int(variant_dataset.iloc[i, 12])
            ref_base = variant_dataset.iloc[i, 13]
            alt_base = variant_dataset.iloc[i, 14]
            MPRA_score = variant_dataset.iloc[i, 22]

            atac_output = AlphaGenomePredictATAC(
                api_key=api_key,
                chromosome=chromosome_number,
                position=variant_position,
                interval_size=interval_size,
                ref=ref_base,
                alt=alt_base,
                ontology_terms=ontology_terms,
                length_keys=length_keys
            )

            fig = AlphaGenomeATACVisualisation(atac_output, interval_size)

            fig.suptitle(
                f"Variant ID: {unique_varID} | Chr: {chromosome_number} | Pos: {variant_position} | Ref: {ref_base} → Alt: {alt_base} | MPRA Score: {MPRA_score}",
                fontsize=6
            )

            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            print(f"Plot {unique_varID} embedded successfully")

    print(f"\n ATAC Signal Report saved to: {output_pdf_path}")

#==========Run Pipeline==========
ExportATACToPDF()
