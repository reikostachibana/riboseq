import ribopy
from ribopy import Ribo
from ribopy_functions import get_cds_range_lookup, get_psite_offset
import numpy as np
import multiprocessing
import time
import pickle
import logging
import math
import gzip

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Function to process a single transcript and return coverage
def process_transcript(transcript, exp, min_len, max_len, alias, cds_range, offset, ribo_path):
    try:
        # Initialize a new Ribo object within the worker process
        if alias == True:
            ribo_object = Ribo(ribo_path, alias=ribopy.api.alias.apris_human_alias)
        else:
            ribo_object = Ribo(ribo_path)
        
        start, stop = cds_range[transcript]

        coverages = []
        for i in range(min_len, max_len + 1): 
            if offset[i] <= start:
                coverage = ribo_object.get_coverage(experiment=exp, range_lower=i, range_upper=i, alias=alias)\
                           [transcript][start - offset[i] : stop - offset[i]]
                coverages.append(coverage)
            else:
                coverage = ribo_object.get_coverage(experiment=exp, range_lower=i, range_upper=i, alias=alias)\
                           [transcript][: stop - offset[i]]
                coverage = np.concatenate((np.zeros(offset[i] - start), coverage))
                coverages.append(coverage)
                
        coverage = sum(coverages, np.zeros_like(coverages[0]))
        
        return transcript, coverage
    except Exception as e:
        logging.error(f"Error processing transcript {transcript}: {e}")
        return transcript, None

# Wrapper function to pass to multiprocessing.Pool.imap_unordered
def process_wrapper(args):
    return process_transcript(*args)

def apris_human_alias(x):
    return x.split("|")[4]

if __name__ == '__main__':
    try:
        # Initialize variables
        output_file = 'HSC_vehicle_vs_GMP_vehicle_negForwarded.pkl.gz'
        experiments = ["GMP_vehicle_A", "GMP_vehicle_B", "HSC_vehicle_A", "HSC_vehicle_B"]
        file_path = "output/genes_list/HSC_vehicle_vs_GMP_vehicle_negForwarded.txt"
        min_len = 20
        max_len = 39
        alias_int = 1 # Mouse
        ribo_path = "all.ribo"
        if alias_int == 1:
            alias = True
            ribo_object = Ribo(ribo_path, alias=ribopy.api.alias.apris_human_alias)
        else: 
            alias = False
            ribo_object = Ribo(ribo_path)

        with open(file_path, 'r') as file:
            gene_list = [line.strip() for line in file]

        transcripts = []
        for gene in gene_list:
            transcript = apris_human_alias(next(str(name) for name in ribo_object.transcript_names if gene in name))
            transcripts.append(transcript)
            
        cds_range = get_cds_range_lookup(ribo_object)

        all_coverage_dict = {}
        for exp in experiments:
            logging.info(f"Starting {exp}...")
            coverage_dict = {}
            offset = get_psite_offset(ribo_object, exp, min_len, max_len)

            # Parallelize transcript processing
            with multiprocessing.Pool() as pool:
                for transcript, coverage in pool.imap_unordered(
                        process_wrapper,
                        [(t, exp, min_len, max_len, alias, cds_range, offset, ribo_path) for t in transcripts]
                    ):
                    # Accumulate the coverage for each transcript in the dictionary
                    coverage_dict[transcript] = coverage
            all_coverage_dict[exp] = coverage_dict

        with gzip.open(output_file, 'wb') as f:
            pickle.dump(all_coverage_dict, f)

        logging.info(f"Saved as {output_file}.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")