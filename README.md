# 🧬🔍 INSPECT pipeline

Nextflow pipeline to process demultiplexed Illumina paired-end 
FASTQ files from INSPECT experiments.

- [Processing steps](#processing-steps)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Outputs](#outputs)

## Processing steps

Per set of reference files:

- Build reference FASTA for forward and reverse reads from CSV files.

Per sample:

1. Filter and trim reads to adapters using `cutadapt`. This ensures reads used downstream 
have the expected features and are trimmed so that the features are in predictable places.
2. Extract clone barcodes and UMIs using `umitools extract`.
3. Identify reference feature combinations using `cutadapt`.
4. Generate tables of counts per UMI and per clone.
5. Count UMIs per reference per clone using `umitools count_tab`.
6. Generate tables of UMI counts per reference.
7. Plot histograms of UMI, clone, and reference count distributions.

### Other steps

1. Get FASTQ quality metrics with `fastqc`.
2. Compile the logs of processing steps into an HTML report with `multiqc`.

## Requirements

### Software

You need to have Nextflow and either `conda` or `mamba` installed on your system. 
If possible, use `mamba` because it will be faster.

#### First time using Nextflow?

If you're at the Crick or your shared cluster has it already installed, try:

```bash
module load Nextflow
```

Otherwise, if it's your first time using Nextflow on your system, you can install it using `conda`:

```bash
conda install -c bioconda nextflow 
```

You may need to set the `NXF_HOME` environment variable. For example,

```bash
mkdir -p ~/.nextflow
export NXF_HOME=~/.nextflow
```

To make this a permanent change, you can do something like the following:

```bash
mkdir -p ~/.nextflow
echo "export NXF_HOME=~/.nextflow" >> ~/.bash_profile
source ~/.bash_profile
```

### References

You also need three files listing your reference sequences (A, promoter barcode, and B) 
in comma-separated (CSV) format, with the sequence name in the left column and the sequence 
in the right column:

| sequence_name | sequence | 
| ------------- | -------- | 
| G1            | ATCCGAGA | 
| G2            | GTCTTAGA | 

The sequence should only contain the unique sequences, i.e. barcodes or spacer sequence. 

## Quick start

Make a [sample sheet (see below)](#sample-sheet)  and, optionally, a [`nextflow.config` file](#inputs) 
in the directory where you want the pipeline to run. Then run Nextflow.

```bash 
nextflow run scbirlab/nf-inspect
```

Each time you run the pipeline after the first time, Nextflow will use a locally-cached version which 
will not be automatically updated. If you want to ensure that you're using the very latest version of 
the pipeline, use the `-latest` flag.

```bash 
nextflow run scbirlab/nf-inspect -latest
```

For help, use `nextflow run scbirlab/nf-inspect --help`.

The first time you run the pipeline on your system, the software dependencies in `environment.yml` 
will be installed. This can take around 10 minutes.

If your run is unexpectedly interrupted, you can restart from the last completed step using the `-resume` flag.

```bash 
nextflow run scbirlab/nf-inspect -resume
```

## Inputs

The following parameters are required:

- `sample_sheet`: path to a CSV containing sample IDs matched with FASTQ filenames, references, 
and adapter sequences
- `fastq`: path to directory containing the FASTQ files (optionally GZIPped)
- `reference`: path to directory containing [reference CSV files](#references)

The following parameters have default values can be overridden if necessary.

- `trim_qual = 10` : For `cutadapt`, the minimum Phred score for trimming 3' calls
- `min_length = 105` : For `cutadapt`, the minimum trimmed length of a read. Shorter reads will be discarded
- `umitools_error = 2`: For `umitools`, the number of errors allowed to correct cell barcodes

The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.

Here is an example of the `nextflow.config` file:

```nextflow
params {
   
    sample_sheet = "/path/to/sample-sheet.csv"
    fastq = "/path/to/fastqs"
    reference = "/path/to/reference"

    // Optional
    trim_qual = 15
    min_length = 90

}
```

Alternatively, you can provide these on the command line:

```bash
nextflow run scbirlab/nf-sbrnaseq -latest \
    --sample_sheet /path/to/sample_sheet.csv \
    --fastq /path/to/fastqs \
    --reference /path/to/reference \
    --trim_qual 15 --min_length 90
``` 

### Sample sheet

The sample sheet is a CSV file providing information about which demultiplexed FASTQ files belong 
to which sample, which references each sample should be mapped to, and the UMI and clone barcode 
scheme for each sample.

The file must have a header with the column names below, and one line per sample to be processed.

- `sample_id`: the unique name of the sample
- `fastq_pattern`: the search glob to find FASTQ files for each sample in `fastq` (see [config](#inputs)). 
The pipleine will look for files matching `<fastq>/*<fastq_pattern>*`, and should match only two files, 
corresponding to paired reads.
- `guideA`: filename of the [CSV containg the references](#references) for this feature. Note that these 
should be in the orientation of the amplicon, i.e. probably reverse complemented. Usually this is how they are 
ordered, so these sequences can be pulled from those files.
- `promoters`: filename of the [CSV containg the references](#references) for this feature
- `guideB`: filename of the [CSV containg the references](#references) for this feature
- `adapter3_read1`: the 3' adapter on the forward read to trim to in [`cutadapt` format](https://cutadapt.readthedocs.io/en/stable/guide.html#specifying-adapter-sequences). Sequence _matching the adapter and everything to the right_ will be removed.
- `adapter3_read2`:  the 3' adapter on the reverse read to trim to in [`cutadapt` format](https://cutadapt.readthedocs.io/en/stable/guide.html#specifying-adapter-sequences). Sequence _matching the adapter and everything to the right_ will be removed.
- `adapter5_read1`: the 5' adapter on the forward read to trim to in [`cutadapt` format](https://cutadapt.readthedocs.io/en/stable/guide.html#specifying-adapter-sequences). Sequence _to the left_ will be removed, but the adapters themselves will be retained.
- `adapter5_read2`:  the adapter on the reverse read to trim to in [`cutadapt` format](https://cutadapt.readthedocs.io/en/stable/guide.html#specifying-adapter-sequences). Sequence _to the left_ will be removed, but the adapters themselves will be retained.
- `umi_read1`: the cell barcode and UMI pattern in [`umitools` regex format](https://umi-tools.readthedocs.io/en/latest/regex.html#regex-regular-expression-mode) for the forward read
- `umi_read2`: the cell barcode and UMI pattern in [`umitools` regex format](https://umi-tools.readthedocs.io/en/latest/regex.html#regex-regular-expression-mode) for the reverse read

Here is an example of the sample sheet:

| sample_id | fastq_pattern | guideA     | promoters     | guideB     | adapter3_read1     | adapter3_read2                      | adapter5_read1                     | adapter5_read2        | umi_read1                                                                           | umi_read2                                                                                            | 
| --------- | ------------- | ---------- | ------------- | ---------- | ------------------ | ----------------------------------- | ---------------------------------- | ---------------------- | ---------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| lib001    | FAU6865A1     | guideA.csv | promoters.csv | guideB.csv | ACAGTN{10}GAGCTCAT | TCTGACCAGGGAAAATAGCCCTCTGACCTGGGGAT | NNNNNNNNTTGTAGCTTCTTTCGAGTACAAAAAC | CGACCGTCTGGAGTACAAAAAC | ^(?P<umi_1>.{8})(?P<discard_1>.{26}).*(?P<discard_2>.{71})(?P<promoterBC_1>.{10})$ | ^(?P<discard_3>.{22}).*(?P<discard_4>.{45})(?P<cell_1>.{10})(?P<discard_5>.{5})(?P<promoterBC_2>.{10})$ |
| lib002    | FAU6865A2     | guideA.csv | promoters.csv | guideB.csv | ACAGTN{10}GAGCTCAT | TCTGACCAGGGAAAATAGCCCTCTGACCTGGGGAT | NNNNNNNNTTGTAGCTTCTTTCGAGTACAAAAAC | CGACCGTCTGGAGTACAAAAAC | ^(?P<umi_1>.{8})(?P<discard_1>.{26}).*(?P<discard_2>.{71})(?P<promoterBC_1>.{10})$ | ^(?P<discard_3>.{22}).*(?P<discard_4>.{45})(?P<cell_1>.{10})(?P<discard_5>.{5})(?P<promoterBC_2>.{10})$ |

## Outputs

Outputs are saved in the same directory as `sample_sheet`. They are organised under 
three directories:

- `processed`: FASTQ files and logs resulting from trimming and UMI extraction
- `mapped`: FASTQ files and logs resulting mapping features
- `counts`: tables and plots corresponding to clone $\times$ reference counts
- `multiqc`: HTML report on processing steps

## Issues, problems, suggestions

Add to the [issue tracker](https://www.github.com/scbirlab/nf-inspect/issues).

## Further help

Here are the help pages of the software used by this pipeline.

- [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html)
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [multiqc](https://multiqc.info/)
- [nextflow](https://www.nextflow.io/docs/latest/index.html)
- [python](https://www.python.org/doc/)
- [matplotlib](https://matplotlib.org/stable/)
- [umitools](https://umi-tools.readthedocs.io/en/latest/index.html)