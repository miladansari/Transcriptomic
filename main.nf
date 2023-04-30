#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl = 2
/*
  ----------------------------------------------------------------------------
                       Main Nextflow script:
*/

log.info """\
    TRANSCRIPTOME RECONSTRUCTION
    Novel short reads assembly improving transcripts reconstruction.
    ================================================================
    Genome assembly: ${params.genome}
    reads          : ${params.input}
    outdir         : ${params.outdir}
    """
    .stripIndent()
 
params.fasta          = params.genomes[ params.genome ]?.fasta
params.gtf            = params.genomes[ params.genome ]?.gtf
params.star_index     = params.genomes[ params.genome ]?.star
params.transcriptome  = params.genomes[ params.genome ]?.transcriptome
params.fastqS_conf= "${baseDir}/conf/fastq_screen.conf"


log.info """\
   
    TRANSCRIPTOME RECONSTRUCTION
    Novel short reads assembly improving transcripts reconstruction.
    ================================================================
    Genome assembly: ${params.genome}
    reads          : ${params.input}
    outdir         : ${params.outdir}
 
    """
    .stripIndent()

include { CONCAT_READS as CONCAT } from './modules/concat_reads'
include { FASTP } from './modules/fastp'
include {STAR_1PASS} from './modules/star_1pass'
include {SJ_1FILTER} from './modules/sj_1filter'
include {STAR_2PASS} from './modules/star_2pass'
include {STRINGTIE}  from './modules/stringTie'
include {SAMTOOLS}   from './modules/samtools'
include {MERGE_STRINGTIE}  from './modules/merge_stringtie' 

workflow{
    //ch_samplesheet = file(params.input, checkIfExists: true)
    reads = Channel
        .fromPath(params.input)  //CSV header: name,fastq[,fastq2],strandness
        .splitCsv(header: true)
        .tap { ch_samplesheet } //
        .map { sample ->
             [ sample.name, file(sample.fastq), file(sample.fastq2) ]
         }
        .groupTuple()


    //store samplesheet metadata in a hash table
    samplesheet = [:]
    ch_samplesheet.subscribe { sample ->
        samplesheet.put(
            sample.name,
            [
                'name': sample.name,
                'strandness': sample.strandness
            ]
        )
    }
star = file(params.star_index)
gtf = file( params.gtf)

CONCAT(reads)
FASTP(CONCAT.out.concated)
STAR_1PASS(star,gtf,FASTP.out.trimmed)
SJ_1FILTER(STAR_1PASS.out.sj.collect())
STAR_2PASS(star,gtf,SJ_1FILTER.out.filtered,FASTP.out.trimmed)
SAMTOOLS(STAR_2PASS.out.pass2)
STRINGTIE (SAMTOOLS.out.paisortprop,gtf)
MERGE_STRINGTIE ( STRINGTIE.out.guided_,gtf)
}
