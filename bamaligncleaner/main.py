import pysam
import logging
from tqdm import tqdm
import sys
from pysam.utils import SamtoolsError


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


def filter_bam(bam, method, output):
    """Filter bam file to remove unaligned references

    Args:
        bam (str): Path to alignement file
        output (str): Path to output alignment file
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
    header = alignment.header.to_dict()
    header["SQ"] = [{"SN": r, "LN": rl} for (r, rl) in zip(refs, reflens)]
    header = pysam.AlignmentHeader.from_dict(header)

    logging.info("Step 4/4: Writing alignment file")
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
    logging.info(
        f"* {total_refs - len(present_refs)} references with unaligned reads were removed from index"
    )
    logging.info(f"* Output bam file written to: {output_message}")
