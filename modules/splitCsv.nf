nextflow.enable.dsl=2

workflow{
    //ch_samplesheet = file(params.input, checkIfExists: true)
    ch_reads = Channel
        .fromPath(params.input)  //CSV header: name,fastq[,fastq2],strandness
        .splitCsv(header: true)
        .tap { ch_samplesheet } //ch_samplesheet mi dà errore
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
 ch_reads.view()
 ch_samplesheet.view()

}