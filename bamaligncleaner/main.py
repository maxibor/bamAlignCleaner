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


def filter_bam(bam, method, output, splits, splitmode):
    """Filter bam file to remove unaligned references

    Args:
        bam (str): Path to alignement file
        method (str): Method to infer references present in the input alignment
                      file
        output (str): Path to output alignment file
        splits (int): Number of output alignment files
        splitmode (str): Method to split the contigs into multiple output
                         alignment files
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
    logging.info("Step 1/4: reading alignment file")
    logging.info(f"* {alignment.mapped} aligned reads")
    logging.info(f"* {total_refs} reference sequences")
    if method.lower() == "index_stat":
        present_refs = set()
        n_reads = {}
        for ref_stat in tqdm(
            alignment.get_index_statistics(), total=total_refs, unit="references"
        ):
            refname = ref_stat[0]
            nb_mapped_reads = ref_stat[1]
            if nb_mapped_reads > 0:
                present_refs.add(refname)
                n_reads[refname] = nb_mapped_reads
        refs = tuple(present_refs)
    elif method.lower() == "parse":
        observed_refs = {}
        n_reads = {}
        for read in tqdm(alignment.fetch(), total=alignment.mapped, unit="reads"):
            if not read.is_unmapped:
                if read.reference_name not in observed_refs:
                    observed_refs[read.reference_name] = []
                if (read.is_paired and
                    read.next_reference_name not in observed_refs[read.reference_name]):
                    observed_refs[read.reference_name].append(read.next_reference_name)
                if read.reference_name not in n_reads:
                    n_reads[read.reference_name] = 1
                else:
                    n_reads[read.reference_name] += 1
        present_refs = set(observed_refs.keys())
        refs = tuple(present_refs)
    reflens = list()

    logging.info("Step 2/4: getting references length")
    for ref in tqdm(refs, unit="references"):
        reflens.append(alignment.get_reference_length(ref))

    logging.info("Step 3/4: recreating header")
    if splits == 1:
        header = alignment.header.to_dict()
        sq = [{"SN": r, "LN": rl} for (r, rl) in zip(refs, reflens)]
        if method.lower() == "parse": 
            sq = extend_sqs(sq, observed_refs, header['SQ'])
        header['SQ'] = sq
        header = pysam.AlignmentHeader.from_dict(header)
    else:  # split into multiple headers
        header = alignment.header.to_dict()
        if splitmode == "contigs":
            splitted_sqs = split_contigs(header["SQ"], splits)
        else:
            splitted_sqs = split_reads(header['SQ'], n_reads, splits)
        splits_ref_map = {
            ref['SN']: i
            for i, h in enumerate(splitted_sqs)
            for ref in h
        }
        headers = [header.copy() for i in range(splits)]
        for i in range(splits):
            headers[i]['SQ'] = extend_sqs(splitted_sqs[i],
                                          observed_refs,
                                          header['SQ'])
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
            outprefix = output
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
        splits(int): Number of splits
    Returns:
        list: List of dictionary of reference contigs
    """
    reflens = np.asarray([ref['LN'] for ref in sq])
    cumsum_rl = np.cumsum(reflens)
    split_sum = cumsum_rl[-1] // splits
    cumsum_splits = np.array(range(1, splits+1)) * split_sum
    indices = np.searchsorted(cumsum_rl, cumsum_splits)
    sq_splits = []
    for _ in range(splits):
        sq_splits.append([])
    for j, ref in enumerate(sq):
        sq_splits[np.argmax(j < indices)].append(ref)
    return sq_splits


def split_reads(sq, ref_stats, splits):
    """Split header of contigs into multiple dictionary so that every one has
       a similar number of aligned reads
    Args:
        sq(dict): Dictionary of reference contigs
        ref_stats(dict): Number of aligned reads per contig
        splits(int): Number of splits
    Returns:
        list: List of dictionary of reference contigs
    """
    nreads = np.asarray([ref_stats[ref['SN']] for ref in sq])
    cumsum_nr = np.cumsum(nreads)
    split_sum = cumsum_nr[-1] // splits
    cumsum_splits = np.array(range(1, splits+1)) * split_sum
    indices = np.searchsorted(cumsum_nr, cumsum_splits)
    sq_splits = []
    for i in range(splits):
        sq_splits.append([])
    for j, ref in enumerate(sq):
        sq_splits[np.argmax(j < indices)].append(ref)
    return sq_splits


def extend_sqs(primary_contigs, mateinfo, orig_sq):
    """Adds to the list of contigs that were a reference for a read the
       contigs that just appeared in the mates.
    Args:
        primary_contigs(list): list of SQ entries in current split
        metainfo(dict): dictionary listing all contigs to which mates were
                        aligned to.
        orig_sq(list): list of all SQ entries from the SAM header
    Returns:
        list: list of SQ entries of SAM header
    """
    contigs_present = set([c['SN'] for c in primary_contigs])
    additional_contigs = []
    for p in primary_contigs:
        if p['SN'] in mateinfo:
            for c in mateinfo[p['SN']]:
                if c not in contigs_present:
                    if c not in additional_contigs:
                        additional_contigs.append(c)
    additional_contigs = set(additional_contigs)
    contigs = primary_contigs.copy()
    for c in orig_sq:
        if c['SN'] in additional_contigs:
            contigs.append(c)
    return contigs
