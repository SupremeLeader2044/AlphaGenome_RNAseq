# Github @SupremeLeader2044 @vanted7580 
#!pip install alphagenome
#imports
from alphagenome import colab_utils
from alphagenome.data import gene_annotation
from alphagenome.data import genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.interpretation import ism
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
from alphagenome.visualization import plot_components

import matplotlib.pyplot as plt
import pandas as pd
import os

#----------Section 1: Functions----------

#==========AlphaGenome Functiion==========
def AlphaGenomePredictVariant(api_key,
                              chromosome,
                              position,
                              interval_size,
                              ref,
                              alt,
                              ontology_terms,
                              length_keys):
    #The input arguments are: API key, chromosome number (e.g.: chr11), ontology term,  
    dna_model = dna_client.create(api_key)  # an API key is requried, it can be found on AlphaGenome's github page
    interval = genome.Interval(chromosome=chromosome,
                               start=max(position - interval_size // 2, 1),
                               end=position + interval_size // 2
                              ).resize(length_keys)  # using information above, find the interval of interest
    variant = genome.Variant(chromosome, position, ref, alt)  # variant (mutation) information

    RNAseq_predicted_output = dna_model.predict_variant(
        interval=interval,
        variant=variant,
        ontology_terms=ontology_terms,
        requested_outputs=[dna_client.OutputType.RNA_SEQ],
    )

    variant_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']  # Variant Score
    RNAseq_variant_scores = dna_model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=[variant_scorer]
    )

    return RNAseq_predicted_output, RNAseq_variant_scores  # return the output
#==========End of Function==========


#==========Variant Score Refining Function==========
def VariantScoreRefining(RNAseq_variant_scores,
                         key_biosamples,
                         raw_score_max,
                         raw_score_min,
                         quantile_score_max,
                         quantile_score_min):

    matrix_raw = pd.DataFrame(
        variant_scorers.tidy_scores([RNAseq_variant_scores], match_gene_strand=True)
    )
    print("matrix_raw created with shape:", matrix_raw.shape)
    print(matrix_raw.head(1))  # Show preview of first few rows

    matrix_refined = (
        matrix_raw
        .loc[matrix_raw['biosample_name'].isin(key_biosamples)]
        .loc[lambda df: df['raw_score'].abs().between(raw_score_min, raw_score_max)]
        .loc[lambda df: df['quantile_score'].abs().between(quantile_score_min, quantile_score_max)]
        .copy()
    )

    # compute per-gene averages
    matrix_refined['avg_raw_score'] = (
        matrix_refined.groupby('gene_name')['raw_score']
                      .transform('mean')
    )
    matrix_refined['avg_quantile_score'] = (
        matrix_refined.groupby('gene_name')['quantile_score']
                      .transform('mean')
    )

    # overwrite the originals
    matrix_refined['raw_score']      = matrix_refined['avg_raw_score']
    matrix_refined['quantile_score'] = matrix_refined['avg_quantile_score']

    # drop helper columns
    matrix_refined.drop(
        ['avg_raw_score', 'avg_quantile_score'],
        axis=1,
        inplace=True
    )

    return matrix_refined
#==========End of Function==========


#==========Excel Processing Function==========
def ProcessExcelVariantData():
    #These variables are customisable
    api_key = "Your API Key"  # You can request a valid API key from AlphaGenome's Github Page
    interval_size = 16384  # placeholder, not really important
    ontology_terms = ["UBERON:0002078"] #You can find the ontology term you need on ontology_term.csv
    length_keys = dna_client.SEQUENCE_LENGTH_1MB  # The size of sequence you wish to investigate

    key_biosamples = [
        "right cardiac atrium",
        "left cardiac atrium",
        "heart right ventricle",
        "heart left ventricle",
        "cardiac muscle cell",
        "cardiac septum",
        "regular cardiac myocyte",
        "Right ventricle myocardium inferior",
        "Right ventricle myocardium superior",
        "left ventricle myocardium inferior",
        "left ventricle myocardium superior",
    ]  # filter: what type of biosample you wish to investigate

    raw_score_max = 100    # filter: maximum variant raw score
    raw_score_min = 0  # filter: minimum variant raw score
    quantile_score_max = 100   # filter: maximum variant quantile score
    quantile_score_min = 0 # filter: minimum variant quantile score

    input_path = input("Enter Input Excel file path: ").replace("\\", "/").replace("\"", "")  # take input and format it
    variant_dataset = pd.read_excel(input_path)  # read the input file
    input_folder = os.path.dirname(input_path)  # identify the input folder
    output_path = os.path.join(input_folder, "Variant_Score_Results.xlsx")  # create an output file in the same input folder

    print((variant_dataset.head(3)))

    combined_results = pd.DataFrame()  # a placeholder for the output matrix

    for i in range(0, len(variant_dataset)):  # this loops reads the data in the input file
        unique_varID = variant_dataset.iloc[i, 3]      # extract variant ID from the file. In this case, it is the 4th column of the input matrix
        chromosome_number = variant_dataset.iloc[i, 11]  # extract chromosome ID from the file. In this case, it is the 12th column of the input matrix
        variant_position = int(variant_dataset.iloc[i, 12])  # extract location of the variant. In this case, it is the 13th column of the input matrix
        ref_base = variant_dataset.iloc[i, 13]  # extract original base. In this case, it is the 14th column of the input matrix
        alt_base = variant_dataset.iloc[i, 14]  # extract altered base. In this case, it is the 15th column of the input matrix
        MPRA_score = variant_dataset.iloc[i, 22]  # extract MPRA score. In this case, it is the 23th column of the input matrix

        i_predicted_output, ith_variant_score = AlphaGenomePredictVariant(
            api_key=api_key,
            chromosome=chromosome_number,
            position=variant_position,
            interval_size=interval_size,
            ref=ref_base,
            alt=alt_base,
            ontology_terms=ontology_terms,
            length_keys=length_keys
        )

        i_matrix_refined = VariantScoreRefining(
            RNAseq_variant_scores=ith_variant_score,
            key_biosamples=key_biosamples,
            raw_score_max=raw_score_max,
            raw_score_min=raw_score_min,
            quantile_score_max=quantile_score_max,
            quantile_score_min=quantile_score_min
        )

        varID_dataframe = pd.DataFrame({'unique_varID': [unique_varID] * len(i_matrix_refined)}) #extract variant ID from your input excel sheet
        mpra_dataframe = pd.DataFrame({'MPRA_score':    [MPRA_score]   * len(i_matrix_refined)}) #You may or may not have this score, change this line accordingly

        i_matrix_final = pd.concat(
            [varID_dataframe.reset_index(drop=True),
             i_matrix_refined.reset_index(drop=True),
             mpra_dataframe.reset_index(drop=True)],
            axis=1
        )

        # ====== SORT BY gene_name HERE ======
        i_matrix_final = (
            i_matrix_final
            .sort_values(by='gene_name')
            .reset_index(drop=True)
        )

        combined_results = pd.concat([
            combined_results,
            i_matrix_final,
            pd.DataFrame([{}])  # blank row
        ], ignore_index=True)

    print(combined_results)
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        combined_results.to_excel(writer, sheet_name="AllVariants", index=False)
#==========End of Function==========

#----------Section 2: Running the Code----------
ProcessExcelVariantData()
#----------End of Section 2----------