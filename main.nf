#!/usr/bin/env nextflow

/*
========================================================================================
   INSPECT Workflow
========================================================================================
   Github   : https://github.com/scbirlab/nf-inspect
   Contact  : Eachan Johnson <eachan.johnson@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2
pipeline_title = """\
                 S C B I R   I N S P E C T   P I P E L I N E
                 ===========================================
                 """
                 .stripIndent()

/*
========================================================================================
   Help text
========================================================================================
*/
if ( params.help ) {
   println pipeline_title + """\

         S C B I R   I N S P E C T   P I P E L I N E
         ===========================================

         Nextflow pipeline to process demultiplexed Illumina paired-end 
         FASTQ files from INSPECT experiments.

         Usage:
            nextflow run sbcirlab/nf-inspect --help
            nextflow run sbcirlab/nf-inspect --sample_sheet <csv> --fastq <dir> --reference <dir>
            nextflow run sbcirlab/nf-inspect -c <config-file>

         Required parameters:
            sample_sheet         Path to a CSV containing sample IDs matched with FASTQ filenames, reference filenames, adapter sequences.
            fastq                Path to directory containing the FASTQ file.
            reference            Path to directory containing references to map to.

         Optional parameters (with defaults):
            trim_qual = ${params.trim_qual}             For `cutadapt`, the minimum Phred score for trimming 3' calls
            min_length = ${params.min_length}           For `cutadapt`, the minimum trimmed length of a read. Shorter reads will be discarded
            umitools_error = ${params.umitools_error}         For `umitools`, the number of errors allowed to correct cell barcodes

         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """.stripIndent()
   System.exit(0)
}

/*
========================================================================================
   Check parameters
========================================================================================
*/
if ( !params.sample_sheet ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to sample_sheet")
}
if ( !params.fastq ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to directory containing FASTQ files")
}
if ( !params.reference ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to directiory containing reference files")
}

wd = file(params.sample_sheet)
working_dir = wd.getParent()

processed_o = "${working_dir}/outputs/processed"
mapped_o = "${working_dir}/outputs/mapped"
counts_o = "${working_dir}/outputs/counts"
multiqc_o = "${working_dir}/outputs/multi_qc"

dirs_to_make = [processed_o, mapped_o, counts_o, multiqc_o]

log.info pipeline_title + """\
         inputs
            sample sheet   : ${params.sample_sheet}
            fastq directory: ${params.fastq}
            reference dir. : ${params.reference}
         trimming 
            quality        : ${params.trim_qual}
            minimum length : ${params.min_length}
         UMI detection
            Error number   : ${params.umitools_error}
         output            
            Processed      : ${processed_o}
            Mapped         : ${mapped_o}
            Counts         : ${counts_o}
            MultiQC        : ${multiqc_o}
         """
         .stripIndent()

log.info  """
          Making directories: 
          """.stripIndent()

dirs_to_make.each { 
   log.info "${it}: " 
   log.info file(it).mkdirs() ? "OK" : "Cannot create directory: ${it}"
}

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

   Channel.fromPath( params.sample_sheet, 
                     checkIfExists: true )
      .splitCsv( header: true )
      .set { csv_ch }

   csv_ch
      .map { tuple( it.sample_id, 
                    "${it.guideA}-${it.promoters}-${it.guideB}",  // Reference ID
                    tuple( it.adapter3_read1, 
                           it.adapter3_read2,
                           it.adapter5_read1, 
                           it.adapter5_read2 ),
                    tuple( it.umi_read1, 
                           it.umi_read2 ),
                    file( "${params.fastq}/*${it.fastq_pattern}*",
                           checkIfExists: true ).sort() ) }
      .set { sample_ch }

   csv_ch
      .map { tuple( "${it.guideA}-${it.promoters}-${it.guideB}",  // Reference ID
             file( "${params.reference}/${it.guideA}",
                   checkIfExists: true),
             file( "${params.reference}/${it.promoters}",
                   checkIfExists: true),
             file( "${params.reference}/${it.guideB}",
                   checkIfExists: true) ) }
      .unique()
      .set { reference_ch }

   sample_ch | FASTQC 

   sample_ch | TRIM_CUTADAPT 
   TRIM_CUTADAPT.out.main | UMITOOLS_WHITELIST
   // TRIM_CUTADAPT.out.main.join( UMITOOLS_WHITELIST.out.main, 
   //                              by: 1 ) | UMITOOLS_EXTRACT

   TRIM_CUTADAPT.out.main | UMITOOLS_EXTRACT

   reference_ch | BUILD_CUTADAPT_REF 
   UMITOOLS_EXTRACT.out.main
      .combine( BUILD_CUTADAPT_REF.out.fastas, 
                by: 0 ) 
   | CUTADAPT_DEMUX

   CUTADAPT_DEMUX.out.main 
   | FASTQ2TAB 
   | READS_PER_UMI_AND_PER_CLONE_AND_PER_GUIDE
   
   FASTQ2TAB.out | UMITOOLS_COUNT_TAB
   UMITOOLS_COUNT_TAB.out.main
      .combine( BUILD_CUTADAPT_REF.out.combos, 
                by: 0 ) 
   | COUNTS_PER_GUIDE

   COUNTS_PER_GUIDE.out.main
      .combine( READS_PER_UMI_AND_PER_CLONE_AND_PER_GUIDE.out, 
                by: 0 )
      .combine( UMITOOLS_COUNT_TAB.out.main_to_plot,
                by: 0 )
      .set { plots_input }
   plots_input | PLOTS
   plots_input
      .map { it[0..1] } 
   | PLOT_READS_VS_UMIS

   TRIM_CUTADAPT.out.multiqc_logs
      .concat(
         CUTADAPT_DEMUX.out.multiqc_logs,
         FASTQC.out.multiqc_logs 
      )
      .flatten()
      .unique()
      .collect()
   | MULTIQC

}


// Do quality control checks
process FASTQC {

   label 'med_mem'

   tag "${sample_id}" 

   input:
   tuple val( sample_id ), val( reference_id ), val( adapter ), val( umis ), path( reads )

   output:
   tuple val( sample_id ), path ( "*.zip" ), emit: logs
   path "*.zip", emit: multiqc_logs

   script:
   """
   zcat ${reads} > ${sample_id}.fastq
   fastqc --noextract --memory 10000 --threads ${task.cpus} ${sample_id}.fastq
   rm ${sample_id}.fastq
   """
   stub:
   """
   zcat ${reads} | head -n1000 > ${sample_id}.fastq
   fastqc --noextract --memory 10000 --threads ${task.cpus} ${sample_id}.fastq
   rm ${sample_id}.fastq
   """

}


// Trim adapters from reads
process TRIM_CUTADAPT {

   tag "${sample_id}"

   label 'big_cpu'

   publishDir( processed_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), val( reference_id ), val( adapters ), val( umis ), path( reads, stageAs: "???/*" )
   
   output:
   tuple val( sample_id ), val( reference_id ), val( adapters ), val( umis ), path( "*.trimmed.fastq.gz" ), emit: main
   tuple val( sample_id ), path( "*.log" ), emit: logs
   path "*.log", emit: multiqc_logs

   script:
   """
   for i in \$(seq 1 2)
   do
      cat */*_R"\$i"*.fastq.gz > ${sample_id}_3prime_R"\$i".fastq.gz
   done

   cutadapt \
		-a '${adapters[0]}' \
      -A '${adapters[1]}' \
		-q ${params.trim_qual} \
      --no-indels \
		--nextseq-trim ${params.trim_qual} \
		--minimum-length ${params.min_length} \
		--report full \
      --action trim \
		--discard-untrimmed \
      -j ${task.cpus} \
		-o ${sample_id}_5prime_R1.fastq.gz \
      -p ${sample_id}_5prime_R2.fastq.gz \
		${sample_id}_3prime_R?.fastq.gz \
      > ${sample_id}.3prime.cutadapt.log

   cutadapt \
		-g '${adapters[2]}' \
      -G '${adapters[3]}' \
      --no-indels \
		--report full \
      --action retain \
		--discard-untrimmed \
      -j ${task.cpus} \
		-o ${sample_id}_R1.trimmed.fastq.gz \
      -p ${sample_id}_R2.trimmed.fastq.gz \
		${sample_id}_5prime_R?.fastq.gz \
      > ${sample_id}.5prime.cutadapt.log

   rm *_?prime_R?.fastq.gz

   """
}


// Identify cell barcodes
process UMITOOLS_WHITELIST {

   label 'big_mem'
   time '24h'

   tag "${sample_id}"

   publishDir( processed_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), val( reference_id ), val( adapter ), val( umis ), path( reads )

   output:
   tuple path( "*.txt" ), val( sample_id ), emit: main
   tuple path( "*.png" ), val( "*.tsv" ), emit: plots
   path "*.log" , emit: logs

   script:
   """
   umi_tools whitelist \
      --knee-method distance \
      --method umis \
		--bc-pattern '${umis[0]}' \
		--bc-pattern2 '${umis[1]}' \
      --extract-method regex \
      --error-correct-threshold ${params.umitools_error} \
      --ed-above-threshold correct \
      --plot-prefix ${sample_id}.whitelist \
		--log ${sample_id}.whitelist.log \
		--stdin ${reads[0]} \
		--read2-in ${reads[1]} \
      --stdout ${sample_id}.whitelist.txt
   """
}


// Extract cell barcodes and UMIs
process UMITOOLS_EXTRACT {

   time '6h'

   tag "${sample_id}"

   publishDir( processed_o, 
               mode: 'copy', 
               pattern: "*.extract{.log,ed.fastq.gz}" )

   input:
   tuple val( sample_id ), val( reference_id ), val( adapter ), val( umis ), path( reads )
   // tuple val( sample_id ), val( reference_id ), val( adapter ), val( umis ), path( reads ), path ( whitelist )
   output:
   tuple val( reference_id ), val( sample_id ), path( "*.extracted.fastq.gz" ), emit: main
   tuple val( sample_id ), path( "*.log" ), emit: logs

   script:
   """
   umi_tools extract \
		--bc-pattern "${umis[0]}" \
		--bc-pattern2 "${umis[1]}" \
      --extract-method regex \
      --quality-filter-mask ${params.trim_qual} \
      --quality-encoding phred33 \
      --log ${sample_id}.extract.log \
		--stdin ${reads[0]} \
		--read2-in ${reads[1]} \
		--stdout ${sample_id}_R1.extracted.fastq.gz \
		--read2-out ${sample_id}_R2.extracted.fastq.gz

   """
}

// Build a reference genome for use by Cutadapt.
process BUILD_CUTADAPT_REF {

   tag "${reference_id}"
   
   input:
   tuple val( reference_id ), path( guideA ), path( promoters ), path( guideB )

   output:
   tuple val( reference_id ), path( "*.fasta" ), emit: fastas
   tuple val( reference_id ), path( "*.txt" ), emit: combos

   script:
   """
   for f in ${promoters} ${guideB}
   do
      fclean=\$f.clean
      cat \$f | tr \$'\r' \$'\n' | grep -v '^\$' > \$fclean
      cat \
         <(head -n1 \$fclean) \
         <(paste -d, \
            <(cut -d, -f1 <(tail -n+2 \$fclean)) \
            <(cut -d, -f2 <(tail -n+2 \$fclean) \
               | tr ATCGatcg TAGCtagc | rev)) \
         > \$f.rc.csv
   done   

   for f in ${guideA} ${promoters} ${promoters}.rc.csv ${guideB}.rc.csv
   do
      cat \$f \
      | tr \$'\r' \$'\n' \
      | grep -v '^\$' \
      | sed 's/^/lib_1,/' \
      > \$f.x.tsv
   done

   join -t, --header ${guideA}.x.tsv ${promoters}.x.tsv \
   | cut -d, -f 2- \
   > referenceA.csv

   awk 'BEGIN{ FS=","; OFS=ARGV[1]; ARGV[1]="" }; NR>1 { print ">"\$1"|"\$3,\$2\$4 }' \$'\\n' \
      referenceA.csv \
   > referenceA.fasta

   join -t, --header ${promoters}.rc.csv.x.tsv ${guideB}.rc.csv.x.tsv \
      | cut -d, -f 2- \
   > referenceB.csv
   
   awk 'BEGIN{ FS=","; OFS=ARGV[1]; ARGV[1]="" }; NR>1 { print ">"\$1"|"\$3,\$4\$2 }' \$'\\n' \
      referenceB.csv \
   > referenceB.fasta

   # generate all unique combos
   join -t, --header ${guideA}.x.tsv ${promoters}.x.tsv \
   | join -t, --header - ${guideB}.rc.csv.x.tsv \
   | cut -d, -f 2- \
   > reference-combo.csv

   cat \
      <(printf 'guide_name\\n') \
      <(awk 'BEGIN{ FS="," }; NR>1 { print \$1"|"\$3"@"\$3"|"\$5 }' reference-combo.csv) \
   > reference-combo.txt

   """
}

// Trim adapters from reads
process CUTADAPT_DEMUX {

   tag "${sample_id}-${reference_id}" 

   label 'big_cpu'

   publishDir( path: mapped_o, 
               mode: 'copy', 
               pattern: "*.{fastq.gz,log}" )

   input:
   tuple val( reference_id ), val( sample_id ), path( reads ), path( fastas )

   output:
   tuple val( reference_id ), val( sample_id ), path( "*.matched.fastq.gz" ), emit: main
   tuple path( "*.log" ),  path( "*.json" ), emit: logs
   path "*.log", emit: multiqc_logs

   script:
   """
   cutadapt \
		-g '^file:${fastas[0]}' \
      -G '^file:${fastas[1]}' \
      -e 1 \
      -j ${task.cpus} \
      --no-indels \
		--report full \
      --action lowercase \
      --rename '{id} {r1.adapter_name}@{r2.adapter_name} {comment}' \
		--discard-untrimmed \
		-o ${sample_id}_R1.matched.fastq.gz \
      -p ${sample_id}_R2.matched.fastq.gz \
      --json ${sample_id}.matched.cutadapt.json \
		${reads} \
      > ${sample_id}.matched.cutadapt.log

   """
}


process FASTQ2TAB {

   tag "${sample_id}"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( reference_id ), val( sample_id ), path( fastqs )

   output:
   tuple val( reference_id ), val( sample_id ), path( "*.tab.tsv" )

   script:
   """
   zcat ${fastqs[0]} \
      | awk '(NR + 3) % 4 == 0' \
      | tr ' ' \$'\t' \
      | cut -f1-2 \
      | sort -k2 \
      > ${sample_id}.tab.tsv
   
   ## Hack because of bug in `umitools count_tab`. It expects read_id_UMI_CB 
   ## when `umitools extract` makes read_id_CB_UMI (!!!)
   #f=${sample_id}.tab0.tsv
   #paste <(paste -d_ <(cut -d_ -f1 \$f) <(cut -f1 \$f | cut -d_ -f3) <(cut -d_ -f2 \$f)) \
   #   <(cut -f2 \$f) \
   #   > ${sample_id}.tab.tsv

   #rm \$f
   """
}


// Count unique UMIs per cell per gene
process UMITOOLS_COUNT_TAB {

   tag "${sample_id}"

   label 'big_mem'
   time '48h'

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( reference_id ), val( sample_id ), path( tabfile )

   output:
   tuple val( reference_id ), val( sample_id ), path( "*.umitools_count.tsv" ), emit: main
   tuple val( sample_id ), path( "*.umitools_count.tsv" ), emit: main_to_plot
   tuple val( reference_id ), val( sample_id ), path( "*.umitools_count.log" ), emit: logs

   script:
   """
   umi_tools count_tab \
		--per-cell \
		--stdin ${tabfile} \
      --stdout ${sample_id}.umitools_count0.tsv \
      --log ${sample_id}.umitools_count.log

   ## Hack for older versions of UMI-tools
   #tail -n+2 ${sample_id}.umitools_count0.tsv \
   #   | sed 's/^b'\\''//;s/'\\''\\t/\\t/' \
   #   > ${sample_id}.umitools_count0.tsv.tail
   NLINES=\$(cat ${sample_id}.umitools_count0.tsv | wc -l)

   printf 'clone_bc.guide_name\\tsample_id\\tclone_bc\\tguide_name\\tumi_count\\n' \
   > ${sample_id}.umitools_count-a.tsv
   paste \
      <(cut -f1-2 ${sample_id}.umitools_count0.tsv | tr \$'\\t' .) \
      <(yes ${sample_id} | head -n \$NLINES) \
      ${sample_id}.umitools_count0.tsv \
   | sort -k1 \
   >> ${sample_id}.umitools_count-a.tsv

   paste -d. \
      <(cut -f1 ${tabfile} | cut -d_ -f2) \
      <(cut -f2 ${tabfile}) \
   > ${sample_id}.read_count0.tsv

   printf 'clone_bc.guide_name\\tread_count\\n' \
   > ${sample_id}.read_count.tsv
   sort ${sample_id}.read_count0.tsv \
   | uniq -c \
   | awk -F' ' -v OFS=\$'\\t' '{ print \$2,\$1 }' \
   | sort -k1 \
   >> ${sample_id}.read_count.tsv

   join --header ${sample_id}.umitools_count-a.tsv ${sample_id}.read_count.tsv \
      | tr ' ' \$'\\t' \
      | cut -f2- \
   > ${sample_id}.umitools_count.tsv

   rm ${sample_id}.read_count0.tsv
   """
}


process READS_PER_UMI_AND_PER_CLONE_AND_PER_GUIDE {

   tag "${sample_id} - ${reference_id}"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( reference_id ), val( sample_id ), path( tabfile )

   output:
   tuple val( sample_id ), path( "*.umi.tsv" ), path( "*.clone_bc.tsv" ), path( "*.guide_name.tsv" )

   script:
   """
   cut -f1 ${tabfile} | cut -d _ -f2 > umi.tsv
   cut -f1 ${tabfile} | cut -d _ -f3 > clone_bc.tsv
   cut -f2 ${tabfile} > guide_name.tsv

   for f in umi.tsv clone_bc.tsv guide_name.tsv
   do
      BASENAME=\$(basename \$f .tsv)
      NLINES=\$(cat \$f | wc -l)

      printf 'sample_id\\t'\$BASENAME'\\t'\$BASENAME'_read_count\\n' \
         > ${sample_id}.\$BASENAME.tsv

      sort \$f | uniq -c \
         | awk -F' ' '{ print "${sample_id}\\t"\$2"\\t"\$1 }' \
         | sort -k3 -n \
         >> ${sample_id}.\$BASENAME.tsv
   done

   """
}


process COUNTS_PER_GUIDE {

   tag "${sample_id}-${reference_id}"

   label 'big_mem'

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( reference_id ), val( sample_id ), path( countfile ), path( combos )

   output:
   tuple val( sample_id ), path( "*.per-guide.tsv" ), emit: main
   tuple val( sample_id ), path( "*.unexpected-guides.tsv" ), emit: unexpected

   script:
   """
   #!/usr/bin/env python

   import sys

   import pandas as pd
   
   df = pd.read_csv("${countfile}", sep='\\t')
   combos = pd.read_csv("${combos}", sep='\\t')

   grouped =  df.groupby('guide_name')
   results = []

   for guide_name, data in grouped:
      count_df = pd.DataFrame(dict(guide_name=[guide_name],
                                   clone_count=[data.shape[0]],
                                   umi_count=[data['umi_count'].sum()],
                                   read_count=[data['read_count'].sum()]))
      results.append(count_df)

   results = pd.concat(results, axis=0)
   print(results.head(), file=sys.stderr)
   unmerged_guides = results['guide_name'].copy().unique()

   results = (results.merge(combos, on='guide_name', how='right')
                .fillna(0.)
                .assign(mean_umis_per_clone=lambda x: x['umi_count'] / x['clone_count'],
                        mean_reads_per_umi=lambda x: x['read_count'] / x['umi_count'],
                        mean_reads_per_clone=lambda x: x['read_count'] / x['clone_count'],
                        sample_id="${sample_id}")
                .sort_values('mean_umis_per_clone', ascending=False))
   results.to_csv('${sample_id}.per-guide.tsv', sep='\\t', index=False)
   merged_guides = results['guide_name'].copy().unique()

   guide_diff = sorted(set(unmerged_guides) - set(merged_guides))

   print(f"There were {len(guide_diff)} unexpected guides!", 
         file=sys.stderr)
   (pd.DataFrame(dict(unexpected_guides=guide_diff))
      .to_csv('${sample_id}.unexpected-guides.tsv', sep='\\t', index=False))

   """
}


process PLOTS {

   tag "${sample_id}"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( guide_umi_counts ), path( umi_reads ), path( clone_reads ), path( guide_reads ), path( clone_guide_counts )

   output:
   tuple val( sample_id ), path( "*.png" )

   script:
   """
   #!/usr/bin/env python

   import sys

   import matplotlib.pyplot as plt
   import numpy as np
   import pandas as pd

   def loghist(df, col, bins=40, *args, **kwargs):

      # df = df[df[col].astype(str).str.isnumeric()].copy()
      # df[col] = df[col].astype(float)
      _min = df.query(f"{col} > 0.")[col].min()
      _max = max(_min + 1, df[col].max())

      bins = np.geomspace(_min, _max, num=bins)
      counts, bins, patches = ax.hist(df[col], 
                                      bins=bins, 
                                      histtype='stepfilled', 
                                      color='#33BBEE')

      ax.set(xlabel=col, 
             ylabel="Frequency",
             xscale='log', *args, **kwargs)

      max_count = max(counts)
      zero_count = (df[col] == 0.).sum()
      if zero_count > 0 and (zero_count / max_count) < 5.:
         ax.axhline(zero_count, 
                    color='lightgrey',
                    zorder=-5)
      
      ax.set_title(ax.get_title() + f'\\nNo. zeroes: {zero_count}')

      return ax

   
   files = dict(per_guides="${guide_umi_counts}", 
                per_umis="${umi_reads}", 
                per_clone_bc="${clone_reads}",
                per_guides_per_clone="${clone_guide_counts}")
   dfs = {name: pd.read_csv(f, sep='\\t') 
          for name, f in files.items()}

   panel_size = 3.5
   n_bins = 40
   n_rows, n_cols = 1, len(dfs)
   fig, axes = plt.subplots(n_rows, n_cols, 
                            figsize=(n_cols * panel_size, n_rows * panel_size),
                            layout='constrained')

   for ax, (name, df) in zip(axes, dfs.items()):

      print(f"Working on {name}...", file=sys.stderr)

      count_col = [c for c in df if 'count' in c][0]
      print(f"Plotting {name} : {count_col}...", file=sys.stderr)

      ax = loghist(df, count_col, bins=n_bins, 
                   title=f"${sample_id}\\n{name}")
      
   fig.savefig('${sample_id}.histograms.png', dpi=600, bbox_inches='tight')

   name = 'per_guides'
   df = dfs[name]
   cols_to_plot = ['clone_count', 'umi_count', 'mean_umis_per_clone'] #df.columns[-3:].tolist()
   print(df[~df['mean_umis_per_clone'].astype(str).str.isnumeric()], file=sys.stderr)
   n_rows, n_cols = 1, len(cols_to_plot)
   fig, axes = plt.subplots(n_rows, n_cols, 
                            figsize=(n_cols * panel_size, n_rows * panel_size),
                            layout='constrained')

   for ax, count_col in zip(axes, cols_to_plot):

      print(f"Plotting {name} : {count_col}...", file=sys.stderr)
      ax = loghist(df, count_col, bins=n_bins, 
                   title=f"${sample_id}\\n{name}")

   fig.savefig('${sample_id}.histograms-per-guide-per-clone.png', dpi=600, bbox_inches='tight')
   
   """
}


process PLOT_READS_VS_UMIS {

   tag "${sample_id}"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( guide_umi_counts )

   output:
   tuple val( sample_id ), path( "*.png" )

   script:
   """
   #!/usr/bin/env python

   from carabiner.mpl import figsaver, scattergrid
   import pandas as pd
   
   df = (
      pd.read_csv(
         "${guide_umi_counts}", 
         sep='\\t',
      )
      .query("read_count > 0 and umi_count > 0")
   )

   fig, axes = scattergrid(
      df,
      grid_columns=["read_count", "umi_count"],
      log=["read_count", "umi_count"]
   )
   figsaver()(
      fig=fig,
      name='${sample_id}.umi-vs-reads',
   )
   
   """
}


// Make log report
process MULTIQC {

   publishDir( multiqc_o, 
               mode: 'copy' )

   input:
   path '*'

   output:
   tuple path( "*.html" ), path( "multiqc_data" )

   script:
   """
   multiqc .
   """
}

/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/