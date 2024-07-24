import logging
import sys

import numpy as np
import pysam
from pysam.utils import SamtoolsError
from tqdm import tqdm


def bam_index(alignment, filetype):
    """Index bam file
    Args:
        alignment(str): Path to alignment file
        filetype(str): Type of alignment file
    Returns:
        str: path to index file
    """
    basename = ".".join(alignment.split(".")[:-1])
    try:
        if filetype == "bam":
            pysam.index(alignment)
            return f"{basename}.bam.bai"
        elif filetype == "cram":
            pysam.index(alignment)
            return f"{basename}.bam.crai"
    except (ValueError, SamtoolsError) as e:
        logging.error(f"[ERROR]: {e}")
        logging.error(
            f"[ERROR]: An error occured while attempting to index {alignment}. Is it sorted ?"
        )
        sys.exit(1)


def filter_bam(bam, method, output, splits):
    """Filter bam file to remove unaligned references

    Args:
        bam (str): Path to alignement file
        method (str): Method to infer references present in the input alignment
                      file
        output (str): Path to output alignment file
        splits (int): Number of output alignment files
    """

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    mode = {"bam": "rb", "cram": "rc"}
    if bam == "-":
        filetype = "bam"
    else:
        filetype = bam.split(".")[-1]
    try:
        alignment = pysam.AlignmentFile(bam, mode[filetype])
        alignment.check_index()
    except ValueError:
        logging.warning(f"[WARNING]: Indexing {bam}")
        index = bam_index(bam, filetype)
        alignment = pysam.AlignmentFile(bam, mode[filetype], index_filename=index)
    except KeyError:
        logging.error(f"[ERROR]: .{filetype} is not a valid filetype")
        sys.exit(1)

    total_refs = alignment.nreferences
    present_refs = set()
    logging.info("Step 1/4: reading alignment file")
    logging.info(f"* {alignment.mapped} aligned reads")
    logging.info(f"* {total_refs} reference sequences")
    if method.lower() == "index_stat":
        for ref_stat in tqdm(
            alignment.get_index_statistics(), total=total_refs, unit="references"
        ):
            refname = ref_stat[0]
            nb_mapped_reads = ref_stat[1]
            if nb_mapped_reads > 0:
                present_refs.add(refname)
        refs = tuple(present_refs)
    elif method.lower() == "parse":
        for read in tqdm(alignment.fetch(), total=alignment.mapped, unit="reads"):
            if not read.is_unmapped:
                present_refs.add(read.reference_name)
        refs = tuple(present_refs)
    reflens = list()

    logging.info("Step 2/4: getting references length")
    for ref in tqdm(refs, unit="references"):
        reflens.append(alignment.get_reference_length(ref))

    logging.info("Step 3/4: recreating header")
    if splits == 1:
        header = alignment.header.to_dict()
        header["SQ"] = [{"SN": r, "LN": rl} for (r, rl) in zip(refs, reflens)]
        header = pysam.AlignmentHeader.from_dict(header)
    else:  # split into multiple headers
        header = alignment.header.to_dict()
        splitted_sqs = split_contigs(header["SQ"], splits)
        splits_ref_map = {
            ref['SN']: i
            for i, h in enumerate(splitted_sqs)
            for ref in h
        }
        headers = [header.copy() for i in range(splits)]
        for i in range(splits):
            headers[i]['SQ'] = splitted_sqs[i]
        headers = [pysam.AlignmentHeader.from_dict(h) for h in headers]

    logging.info("Step 4/4: Writing alignment file")
    if splits == 1:
        outbam = pysam.AlignmentFile(output, "wb", header=header)
        for ref in tqdm(refs, unit="references"):
            for read in alignment.fetch(ref):
                read_dict = read.to_dict()
                read_out = pysam.AlignedSegment.from_dict(read_dict, header)
                outbam.write(read_out)
        alignment.close()
        outbam.close()
        if output == "-":
            output_message = "STDOUT"
        else:
            output_message = output
    else:  # Split in multiple BAM files
        if output.split(".")[-1] in ['bam', 'cram']:
            outprefix = ".".join(output.split(".")[:-1])
        else:
            outprefix = bam
        outbam = [pysam.AlignmentFile(f"{outprefix}.{i:02d}.{filetype}",
                                      "wb", header=headers[i])
                  for i in range(splits)]
        for ref in tqdm(refs, unit="references"):
            for read in alignment.fetch(ref):
                read_dict = read.to_dict()
                read_out = pysam.AlignedSegment.from_dict(read_dict,
                                                          headers[splits_ref_map[ref]])
                outbam[splits_ref_map[ref]].write(read_out)
        alignment.close()
        for i in range(splits):
            outbam[i].close()
        output_message = output
    logging.info(
        f"* {total_refs - len(present_refs)} references with unaligned reads were removed from index"
    )
    logging.info(f"* Output bam file written to: {output_message}")


def split_contigs(sq, splits):
    """Split header of contigs into multiple dictionary so that every one has
       a similar total length of contigs
    Args:
        sq(dict): Dictionary of reference contigs
    Returns:
        list: List of dictionary of reference contigs
    """
    reflens = np.asarray([ref['LN'] for ref in sq])
    cumsum_rl = np.cumsum(reflens)
    split_sum = cumsum_rl[-1] // splits
    cumsum_splits = np.array(range(1, splits+1)) * split_sum
    indices = np.searchsorted(cumsum_rl, cumsum_splits)
    sq_splits = []
    for i in range(splits):
        sq_splits.append([])
    for j, ref in enumerate(sq):
        sq_splits[np.argmax(j < indices)].append(ref)
    return sq_splits
