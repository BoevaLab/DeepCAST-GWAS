from pathlib import Path
import pandas as pd
import re
import concurrent.futures

from config import GFF_FILE_PATH
from logger import setup_logger

logger = setup_logger(__name__)

def generate_coding_regions_chr(chr: int, gff_df: pd.DataFrame):
    start_col = []
    end_col = []

    for _, row in gff_df.iterrows():
        if row['type'] == 'exon':
            start_col.append(row['start'])
            end_col.append(row['end'])

    chr_col = [chr for _ in range(len(start_col))]

    return pd.DataFrame({'chr':chr_col, 'start':start_col, 'end':end_col}, dtype=int)

def call_generate_coding_regions_chr(tuple):
    return generate_coding_regions_chr(*tuple)

def generate_coding_regions(columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]):

    if Path(GFF_FILE_PATH).exists():
        gff_df = pd.read_csv(GFF_FILE_PATH, sep="\t", comment='#', header=None, names=columns)

    else:
        raise ValueError(f'Specified gff file {GFF_FILE_PATH} does not exist. Refer to src/config.py')

    pattern = re.compile(r"chromosome=(\d+)")
    chr_dfs: list[tuple[int, pd.DataFrame]] = []
    non_std_chr_seqs = []
    non_ref_chr_seqs = []

    for seqid, sub_df in gff_df.groupby('seqid'):
        if str(seqid).split('_')[0] == 'NW' or str(seqid).split('_')[0] == 'NT':
            # Those are not the reference genome but scaffolds and so on
            non_ref_chr_seqs.append(seqid)
            continue
        match = pattern.search(sub_df.iloc[0]['attributes'])
        if match:
            chr = int(match.group(1))
            chr_dfs.append((chr, sub_df))
        else:
            non_std_chr_seqs.append(seqid)
            # print(f'''No match when parsing chromosome in gff file for '''
            #       f'''first row of sequence {seqid}. First row attributes: '''
            #       f'''{sub_df.iloc[0]['attributes']}''')
    
    logger.debug(f"Skipped {len(non_ref_chr_seqs)} sequences since they're not reference genome: \n {non_ref_chr_seqs}")
    logger.debug(f"Skipped {len(non_std_chr_seqs)} sequences that don't match chromosomes we are looking at: \n {non_std_chr_seqs}")

    with concurrent.futures.ProcessPoolExecutor() as executor:
        chr_coding_regions: list[pd.DataFrame] = list(executor.map(call_generate_coding_regions_chr, chr_dfs))

    if chr_coding_regions:
        return pd.concat(chr_coding_regions).reset_index(drop=True)
    else:
        return pd.DataFrame(columns=['chr', 'start', 'end'])