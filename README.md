# Denovo-transcriptome

De novo pipeline for the identification of novel L1 containing transcripts in **genome guided** mode. To run this pipeline the user can create a conda environment from this [environment.yml](https://github.com/miladansari/Transcriptomic/blob/master/environment.yaml).
Concerning gffread is downloadable from [here](https://github.com/gpertea/gffread).

The pipeline includes: reads quality control, preprocessing (trimming and removing rRNA) , two-pass alignment that improves novel splice junction quantification, genome-guided transcriptome assembly approach, and finally StringTie takes as input the list of GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcripts.

## Nextflow Installation
---
Nextflow can be installed in any Conda Environment with `conda install -c bioconda nextflow=22.08.`.
Alternatively,`curl https://get.nextflow.io | bash`. In this case, run Nextflow as `./nextflow run`.

## Getting Started:
---
To run Nextflow, samples and their location must be specified in a csv file. 

The **input.csv** must be generated, comma separated, in the following way, having as header:

| name | fastq1 | fastq2 | strandedness|
| ---- | ------ | ------ | ------------|
| sample01 | /path/to/S01_R1.fastq |/path/to/S01_R2.fastq|forward / reverse / unstranded|
| sample02 | /path/to/S02_R1.fastq |/path/to/S02_R2.fastq|forward / reverse / unstranded|
| sample03 | /path/to/S03_R1.fastq |/path/to/S03_R2.fastq|forward / reverse / unstranded|
| sample ... | ... |(*if single-end, leave this field blank*)|forward / reverse / unstranded|



## Which and Where:

| Directory    | file  | function  |
| ------------- |-------------| -----|
| ./assets      | [human_ribosomal.fa](https://github.com/miladansari/Transcriptomic/blob/master/assets/human_ribosomal.fa)| Common ribosomal 31-mers for BBDuk |
| ./assets      | [multiqc_config.yaml](https://github.com/miladansari/Transcriptomic/blob/master/assets/multiqc_config.yaml)    |MultiQC configuration file|
| ./conf | [base.config](https://github.com/miladansari/Transcriptomic/blob/master/conf/base.config)     |    Configuration file to execute the pipeline with HPC Slurm:Workload Manager  |
| ./conf | [fastq_screen.conf](https://github.com/miladansari/Transcriptomic/blob/master/conf/fastq_screen.config)| Fastqscreen configuration files and references |
| ./conf | [genom.config](https://github.com/miladansari/Transcriptomic/blob/master/conf/genomes.config) | Configuration file containing the path to annotations and fasta files
| ./modules | all the modules.nf | Modules to run the pipeline, the workflow is described below |
|./ | [main.nf](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/main.nf) | Main nextflow script
| ./ | [nextflow.config](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/nextflow.config) | Nextflow configuration file containing all the parameters and profile needed to run the pipeline



## Workflow:
1(1). **parse items emitted by a channel** ([SplitCsv](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/modules/splitCsv.nf)), Default = Run

1(2). **Reads concatenation** ([CONCAT](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/modules/concat_reads.nf)) Default = Run, skipped automatically if not needed


2. **FASTQ data pre-processing** ([FASTP](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/modules/fastp.nf)) quality control, trimming of adapters, filtering by quality, and read pruning. Default = Run


3. **Decontamination Using Kmers** ([BBDuk](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/modules/bbduk.nf)) Search for contaminant sequences, part of the BBTools package. Default = Run


4(1). **Detection of library contamination** ([FASTQSCREEN](https://github.com/miladansari/Denovo-transcriptreconstruction/blob/main/modules/fastqscreen.nf)). Default = Run


4(2). **MultiQC** ([MULTIQC](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/modules/multiqc.nf)) Default = Run 


5. **STAR indexing** (STAR_INDEX) Default = Skip 
STAR Index already provided. To create a custom STAR Index: 
provide genome (fasta) and annotation (gtf) in ./main.nf, and --star_idxist false.


6. **STAR alignment 1st_pass** ([STAR_1PASS](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/modules/star_1pass.nf)) Default = Run.

Edit star_1pass.nf to define custom parameters.


7. **Splice junctions filtering** ([SJ_1FILTER](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/modules/sj_1filter.nf) Default = Run.

Edit sj_1filter.nf to define custom filtering parameters. 


8. **STAR alignment 2st_pass** ([STAR_2PASS](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/modules/star_2pass.nf)) Default = Run.

Edit star_2pass.nf to define custom parameters.


9. **SAMTOOLS filterin BAM** ([SAMTOOLS](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/modules/samtools.nf)) Default = Run 


10. **STRINGTIE Reference Guided Transcript Assembly** ([STRINGTIE](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/modules/stringTie.nf)) Default = Run 

11. **STRINGTIE MERGE** ([MERGE_STRINGTIE](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/modules/merge_stringtie.nf)) Default = Run

12. **MultiQC** ([MULTIQC](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/modules/multiqc.nf)) Default = Run 


## How to Run the Pipeline:

Activate your nextflow conda environment with `conda activate nameCondaEnv`, then:

```
nextflow -C nextflow.config run main.nf \
         --input /path/to/input.csv \
         -w /path/to/work_dir \
         --outdir /path/to/your/results \
         -profile singularity \
         -resume -bg 
```

**Note:**
The Pipeline assumes by default a **Paired-End** library. To work with Single-End files: `--single_end true`.

## List of Params:

Parameter | Default Value | Alternative Value | function
--- | --- | --- | --- 
-C | when defined, any other config file will be overwritten | can be a file in the root directory e.g. [nextflow.config](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/nextflow.config) | defines the path of the main nextflow.config file
--input | must be defined  | can be defined inside [main.nf](https://github.com/miladansari/Denovo-transcript-reconstruction/blob/main/main.nf) (line 51 of the script)| defines the path of the input.csv samplesheet
-w | must be defined | nextflow defines it in root directory of the pipeline | defines the path of the Nextflow work directory
--outdir | ./results | - | defines the path where results will be saved separately from work directory
-bg | optional, but recommended | - | parameter to run Nextflow in background, prevents a broken pipeline in case of disconnection
-resume | - | - | allows for the continuation of a workflow execution

**Important:**  disabling certain processes could cause fatal errors (downstream steps may require previously generated files). 
For instance, STAR is required to perform most of the subsequent processes.
