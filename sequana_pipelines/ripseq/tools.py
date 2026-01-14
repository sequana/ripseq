import pandas as pd
import pylab
from sequana import GFF3
from tqdm import tqdm
import numpy as np

from easydev import do_profile


# Get info to normalised data. correltion will not change but mmIP and mmIn will
def build_normalisation(multiqc_fastqc, outfile="normalisation.txt"):
    """
    Build a normalisation table from MultiQC FastQC output.

    The normalisation factor is computed as the number of reads per sample
    divided by the mean number of reads across all samples. This factor can
    later be used to normalise coverage or signal tracks.

    Parameters
    ----------
    multiqc_fastqc : str
        Path to a MultiQC FastQC summary file (tab-separated) containing
        at least the columns 'Sample' and 'Total Sequences'.
        Typically something like ``analysis/multiqc_fastqc.txt``.
    outfile : str, optional
        Output filename for the normalisation table (CSV format).

    Returns
    -------
    None
        The normalisation table is written to ``outfile``.
    """
    # input is e.g. analysis/multiqc_fastqc.txt
    data = pd.read_csv(multiqc_fastqc, sep="\t")[["Sample", "Total Sequences"]]
    data = data[["_R1" in x for x in data["Sample"]]]
    data["Sample"] = [
        x.replace("_R1", "").replace("_sample_", "_") for x in data["Sample"]
    ]
    data["Total Sequences"] = [round(x / 1000000) for x in data["Total Sequences"]]
    data.columns = ["sample", "reads"]
    data["norm"] = data["reads"] / data["reads"].mean()
    datanorm = data.set_index("sample")
    # divide all coverage by norm
    datanorm.to_csv(outfile, sep=",")


def get_attribute(chrom, position, step, gff):
    """
    Retrieve GFF features overlapping a genomic window.

    A feature is considered overlapping if it:
    - starts within the window
    - ends within the window
    - fully contains the window

    Parameters
    ----------
    chrom : str
        Chromosome or contig name.
    position : int
        Start position of the window.
    step : int
        Window size.
    gff : sequana.GFF3
        GFF3 object containing genome annotations.

    Returns
    -------
    pandas.DataFrame
        Subset of the GFF dataframe containing overlapping features.
    """
    # consider a position-position+step chunk to be in a gene if it overlaps
    # or is included inside
    return gff.df.query(
        "seqid == @chrom and ( \
                            (start >= @position and start <=(@position+@step)) or \
                            (stop >= @position and stop <=(@position+@step)) \
                            or (start<@position and (@position+@step)< stop))"
    )


def get_events(
    infile, chrom, step, gff=None, CIP=0.95, min_IP=10, ratio=1.1, mycase=None
):
    """
    Annotate correlation-based signal windows with genomic features.

    For a given chromosome and experimental condition, this function:
    - loads precomputed window-level statistics
    - filters windows based on correlation, IP signal, and IP/IN ratio
    - annotates each retained window using GFF features
    - returns a final annotated dataframe

    Parameters
    ----------
    chrom : str or int
        Chromosome identifier.
    step : int
        Window size used during signal aggregation.
    gff : sequana.GFF3, optional
        GFF3 annotation object used to retrieve gene attributes.
    CIP : float, optional
        Minimum IP correlation threshold.
    min_IP : float, optional
        Minimum mean IP signal.
    ratio : float, optional
        Minimum ratio between mean IP and mean IN signals.
    mycase : str
        Name of the experimental case (used to locate input files).

    gff = GFF3(...)
    gff.df['gene_name'] = gff.df['Name']

    Returns
    -------
    pandas.DataFrame
        Annotated dataframe containing filtered windows, signal statistics,
        genomic coordinates, and gene annotations.
    """
    df = pd.read_csv(infile)
    # "results/{mycase}/{chrom}.csv")

    df.index = [x * step for x in df.index]
    # this is to avoid division by zero
    # df.mmIP += 1
    df.mmIN += 1
    dd = df.query("CIP>=@CIP and mmIP>=@ratio*mmIN and mmIP>=@min_IP").copy()

    gene_names = []
    gene_ids = []
    starts = []
    stops = []
    for pos in dd.index:
        if gff:
            print(pos)
            attrs = get_attribute(str(chrom), pos, step, gff)[
                ["gene_name", "ID", "start", "stop"]
            ].fillna("")
            print(attrs)
            gene_names.append("-".join(attrs["gene_name"]))
            gene_ids.append("-".join(attrs["ID"]))
        else:
            pass
        starts.append(";".join([str(x) for x in attrs["start"]]))
        stops.append(";".join([str(x) for x in attrs["stop"]]))

    dd["ratio"] = (1 + dd["mmIP"]) / (1 + dd["mmIN"])
    dd["position"] = dd.index
    dd["start"] = starts
    dd["stop"] = stops
    dd["gene_name"] = gene_names
    dd["gene_id"] = gene_ids
    dd["chrom"] = chrom

    return dd


def plot_grouped(IPs, INs, chrom, start, stop, datanorm, logy=False, case=None):
    """
    Plot grouped IP and IN signals over a genomic region.

    Signals are optionally normalised using precomputed normalisation
    factors and can be displayed on a log scale.

    Parameters
    ----------
    IPs : list
        List of IP signal objects supporting a ``values(chrom, start, stop)``
        method.
    INs : list
        List of input (IN) signal objects.
    chrom : str or int
        Chromosome identifier.
    start : int
        Start genomic coordinate.
    stop : int
        Stop genomic coordinate.
    datanorm : pandas.DataFrame
        Normalisation table indexed by sample name with a ``norm`` column.
    logy : bool, optional
        If True, apply log transformation to the signal.
    case : str, optional
        Case name used to retrieve normalisation factors.

    Returns
    -------
    None
        The function produces a matplotlib plot.
    """
    from pylab import plot, log, title, legend, ylim

    patts = ["-", "--", "-.", ":"]
    for i, x in enumerate(INs):
        data = x.values(chrom, start, stop)
        if case:
            data /= datanorm.loc[f"{case}_INPUT_{i+1}", "norm"]
        if logy:
            data = [log(datum + 1) for datum in data]
        plot(data, "y" + patts[i], alpha=0.5, label=f"Input {i+1}")

    for i, x in enumerate(IPs):
        data = x.values(chrom, start, stop)
        if case:
            data /= datanorm.loc[f"{case}_IP_{i+1}", "norm"]
        if logy:
            data = [log(datum + 1) for datum in data]

        plot(data, "r" + patts[i], alpha=0.5, label=f"IP {i+1}")
    title(f"chr{chrom}:{start}-{stop}")
    legend()


def preload_bigwig_arrays(IPs, INs, chrom, chrom_len):
    """
    Load full-chromosome coverage arrays from BigWig files.
    """
    IP_arrays = [np.nan_to_num(bw.values(chrom, 0, chrom_len)) for bw in IPs]
    IN_arrays = [np.nan_to_num(bw.values(chrom, 0, chrom_len)) for bw in INs]
    return np.asarray(IP_arrays), np.asarray(IN_arrays)


def get_data(
    IPs,
    INs,
    chrom_name,
    case,
    start=0,
    stop=None,
    step=1000,
    norm_file="normalisation.txt",
    randomize=False,
):
    """
    Compute RIP-seq reproducibility and enrichment metrics along a chromosome
    using a sliding window approach.

    For each genomic window, this function extracts normalized coverage values
    from IP and INPUT BigWig files, computes replicate-to-replicate correlations,
    and summarizes signal intensity metrics. The resulting values are used to
    identify reproducibly enriched RIP-seq regions.

    Parameters
    ----------
    IPs : list of pyBigWig.BigWig
        List of open BigWig objects corresponding to IP replicates.
    INs : list of pyBigWig.BigWig
        List of open BigWig objects corresponding to INPUT replicates.
    chrom_name : str
        Chromosome or contig name for which the analysis is performed.
    case : str
        Experimental condition identifier (e.g. "24h", "48h"), used for
        normalization lookup.
    start : int, optional
        Genomic start coordinate (0-based) of the analysis region.
        Default is 0.
    stop : int or None, optional
        Genomic end coordinate of the analysis region. If None, the full
        chromosome length is used.
    step : int, optional
        Window size (in base pairs) used for signal extraction.
        Default is 1000 bp.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by genomic window start positions and containing
        the following columns:

        - CIP : float
            Mean pairwise Pearson correlation across IP replicates.
        - CIN : float
            Mean pairwise Pearson correlation across INPUT replicates.
        - mmIP : float
            Mean of the maximum normalized IP signal across replicates.
        - mmIN : float
            Mean of the maximum normalized INPUT signal across replicates.

        Each row corresponds to one genomic window.

    Notes
    -----
    - Coverage values are normalized using sequencing depth factors computed
      from FastQC statistics.
    - Windows are processed independently and do not overlap.
    - Correlation metrics quantify reproducibility across biological replicates
      and are central to the identification of enriched RIP-seq regions.
    """

    CIP, CIN, mmIP, mmIN = [], [], [], []

    print(IPs)
    print(IPs[0])
    chroms = IPs[0].chroms()
    chrom_len = chroms[chrom_name]

    if stop is None:
        stop = chroms[chrom_name]
    assert stop > start

    norm = pd.read_csv(norm_file, index_col="sample")

    # stop /= 1000
    # stop = int(stop)
    for x in tqdm(range(start, stop, step)):
        a, b, c, d = correlation_fast(
            IPs, INs, chrom_name, x, step, case, norm, chrom_len
        )
        CIP.append(a)
        CIN.append(b)
        mmIP.append(c)
        mmIN.append(d)

    df = pd.DataFrame({"CIP": CIP, "CIN": CIN, "mmIP": mmIP, "mmIN": mmIN})
    df.index = list(range(start, stop, step))
    return df


def correlation(IPs, INs, chrom_name, start, L, CASE, normalisation, randomize=False):
    """
    Compute replicate correlation and signal enrichment metrics for a genomic window.

    This function extracts coverage values from multiple IP and INPUT BigWig
    replicates over a specified genomic interval, applies sequencing-depth
    normalization, and computes summary statistics that quantify both
    reproducibility and enrichment of RIP-seq signal.

    Parameters
    ----------
    IPs : list of pyBigWig.BigWig
        List of open BigWig objects corresponding to IP replicates.
    INs : list of pyBigWig.BigWig
        List of open BigWig objects corresponding to INPUT replicates.
    chrom_name : str
        Chromosome or contig name.
    start : int
        Genomic start coordinate (0-based) of the window.
    L : int
        Window length in base pairs.
    CASE : str
        Experimental condition identifier (e.g. "24h", "48h"), used to retrieve
        normalization factors.

    Returns
    -------
    tuple of float
        A 4-element tuple containing:

        - corrIP : float
            Mean pairwise Pearson correlation across IP replicates.
        - corrIN : float
            Mean pairwise Pearson correlation across INPUT replicates.
        - mean_maxIP : float
            Mean of the maximum normalized IP signal across replicates.
        - mean_maxIN : float
            Mean of the maximum normalized INPUT signal across replicates.

    Notes
    -----
    - Coverage values are extracted directly from BigWig files using pyBigWig.
    - If the requested window extends beyond the chromosome boundary, the
      window size is truncated to the chromosome end.
    - Normalization factors are read from the file ``normalisation.txt`` and
      are derived from sequencing depth estimates.
    - Pearson correlation is computed across positions within the window,
      and correlations are averaged across all replicate pairs.
    - Missing values are replaced by zeros prior to correlation computation.
    """

    data = []
    colnames = []
    chroms = IPs[0].chroms()

    # for the randomization
    N1 = int(L / 4)
    N2 = int(2 * L / 4)
    N3 = int(3 * L / 4)

    def shift(datum, i):
        if i == 0:
            datum = datum[N1:] + datum[0:N1]
        elif i == 1:
            datum = datum[N2:] + datum[0:N2]
        elif i == 2:
            datum = datum[N3:] + datum[0:N3]
        return datum

    for i, x in enumerate(IPs):
        if start + L > chroms[chrom_name]:
            L_last_chunk = chroms[chrom_name] - start
            datum = x.values(chrom_name, start, start + L_last_chunk)
            if randomization:
                datum = shift(datum, i)
            data.append(datum)
        else:
            datum = x.values(chrom_name, start, start + L)
            if randomization:
                datum = shift(datum, i)
            data.append(datum)
        colnames.append(f"IP{i+1}")

    for i, x in enumerate(INs):
        if start + L > chroms[chrom_name]:
            L_last_chunk = chroms[chrom_name] - start
            datum = x.values(chrom_name, start, start + L_last_chunk)
            if randomization:
                datum = shift(datum, i)
            data.append(datum)
        else:
            datum = x.values(chrom_name, start, start + L)
            if randomization:
                datum = shift(datum, i)
            data.append(datum)
        colnames.append(f"IN{i+1}")

    df = pd.DataFrame(data, index=colnames).T.fillna(0)

    for i, x in enumerate(INs):
        df[f"IN{i}"] = df[f"IN{i}"] / normalisation.loc[f"{CASE}_INPUT_{i}", "norm"]
        df[f"IP{i}"] = df[f"IP{i}"] / normalisation.loc[f"{CASE}_IP_{i}", "norm"]

    corrIP = (
        df[[x for x in colnames if x.startswith("IP")]].corr().fillna(0).mean().mean()
    )
    corrIN = (
        df[[x for x in colnames if x.startswith("IN")]].corr().fillna(0).mean().mean()
    )

    mean_maxIP = df[[x for x in colnames if x.startswith("IP")]].max().mean()
    mean_maxIN = df[[x for x in colnames if x.startswith("IN")]].max().mean()

    return corrIP, corrIN, mean_maxIP, mean_maxIN


def correlation_fast(IPs, INs, chrom, start, L, CASE, norm, chrom_len):
    end = min(start + L, chrom_len)

    arrays = []
    norms = []

    for i, bw in enumerate(IPs):
        arrays.append(bw.values(chrom, start, end))
        norms.append(norm.loc[f"{CASE}_IP_{i+1}", "norm"])

    for i, bw in enumerate(INs):
        arrays.append(bw.values(chrom, start, end))
        norms.append(norm.loc[f"{CASE}_INPUT_{i+1}", "norm"])

    X = np.nan_to_num(np.array(arrays)) / np.array(norms)[:, None]

    IP = X[: len(IPs)]
    IN = X[len(IPs) :]

    corrIP = np.nan_to_num(np.corrcoef(IP), nan=0.0).mean()
    corrIN = np.nan_to_num(np.corrcoef(IN), nan=0.0).mean()

    mmIP = IP.max(axis=1).mean()
    mmIN = IN.max(axis=1).mean()

    return corrIP, corrIN, mmIP, mmIN
