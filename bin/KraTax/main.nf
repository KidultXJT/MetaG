#! /usr/bin/env nextflow

/* ------------ Input parameters ------------ */
params.query      = false
params.format     = ""
params.ident      = ""
params.smp        = ""
params.inPath     = ""
params.binOutPath = "04.BinGenome"
params.taxOutPath = "04.KraTax"

/* ------------ MaxBin ------------ */
params.bin     = false
params.Soft    = "maxbin"
params.minLen  = 1000

/* ------------ Kraken ------------ */
params.Krakendb = "/Bio/Database/Kraken_db/forPipe/kraken_1803_mini_10b"

/* ------------ Slurm params ------------ */
params.help  = false
params.xN    = false
params.fN    = false
params.local = false

if (params.help) {
    log.info "                                                                                                                                                "
    log.info "====================================================================                                                                            "
    log.info "          Sequence Taxonomy Annot by Kraken  (Version 0.1)                                                                                      "
    log.info "====================================================================                                                                            "
    log.info "                                                                                                                                                "
    log.info "USAGE:                                                                                                                                          "
    log.info "                                                                                                                                                "
    log.info "nextflow run main.nf                                                                                                                            "
    log.info " --query      input (contigs)fasta file [None]                                                                                                  "
    log.info " --smp        sample list. Example A,B,C [None]                                                                                                 "
    log.info " --ident      indentify  [None]                                                                                                                 "
    log.info " --format     input File Format for Quant [Default : bam | sam | fq ]                                                                           "
    log.info " --inPath     input FASTQ file / Sam File [None]                                                                                                "
    log.info " --binOutPath Bin Output Path    [04.BinGenome]                                                                                                 "
    log.info " --taxOutPath Kraken Output Path [04.KraTax]                                                                                                    "
    log.info "                                                                                                                                                "
    log.info " ============================================                                                                                                   "
    log.info "                 Slurm Args                                                                                                                     "
    log.info " ============================================                                                                                                   "
    log.info " --fN         compute_R830 | compute_R430 [compute_R430]                                                                                        "
    log.info " --xN         Except Slurm-R430-Node[1,2,3,4,5,6,7,8] [Slurm-R430-Node[?]]                                                                      "
    log.info "                                                                                                                                                " 
    log.info " ============================================                                                                                                   "
    log.info "               NextFlow Args                                                                                                                    "
    log.info " ============================================                                                                                                   "
    log.info " -resume      |NextFlow Args| :: Execute the script using the cached results, useful to continue executions that was stopped by an error [false]"
    log.info " -qs          |NextFlow Args| :: Max number of processes that can be executed in parallel by each executor [false]                              "
    log.info " -with-report |NextFlow Args| :: Create processes execution html report [Add HTML output FileName]                                              "
    log.info " -with-trace  |NextFlow Args| :: Create processes execution tracing file [Add TEXT output FileName]                                             "
    log.info "                                                                                                                                                " 
    exit 1
}


/* -------------- Necessary Parameter ------------- */
/*
if (!params.smp) {
    log.info""
    log.info"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    log.info"x Sample ID for 'smp'                        x"
    log.info"x                                            x"
    log.info"x inPath/{smp}_R{1,2}.fq.gz.${params.ident}  x"
    log.info"x inPath/{smp}.sam.gz                        x"
    log.info"x                                            x"
    log.info"x Example:                                   x"
    log.info"x ~/A_R1.fq.gz                               x"
    log.info"x smp = 'A' HERE !!!                         x"
    log.info"x ident = ''                                 x"
    log.info"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    exit 1
}

if (!params.inPath) {
    log.info""
    log.info"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    log.info"x Your Input File Path for 'inPath'     x"
    log.info"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    exit 1
}
*/
/* ---------------------------------- Input ------------------------------------- */
/* variables */
smplst = params.smp.tokenize(',')
Channel.fromPath("$params.query").into {binfa;krafa}
/*
if (params.format=="fq"){
    Channel
        .fromFilePairs("${params.inPath}/*_{1,2}.fq.gz${params.ident}")
        .filter {it[0] in smplst}.into {quantBin;quantTax}
}


if (params.format=="bam"){
    Channel
        .fromPath("${params.inPath}/*.bam")
        .filter {it[0] in smplst}.into {quantBin;quantTax}
}

if (params.format=="sam"){
    Channel
        .fromPath("${params.inPath}/*.sam.gz")
        .filter {it[0] in smplst}.into {quantBin;quantTax}
}
*/
/*quantTax.println()*/
/*
process bin {
    tag {params.pject}
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J Kt_BIN"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J Kt_BIN"
    }
    cpus 30
    memory 30
    errorStrategy 'retry'
    maxErrors 1
    
    input:
        file quant from quantBin.map {it[1]}.collect()
        file fa  from binfa
    output:
        file("work.sh") into binwork
    
    when:
        params.bin
 
        script:
        if (params.format=="bam") {
        """
        # Make Abundance Information
        echo "# Format is Bam. Script :" >> work.sh
        cat .command.sh >> work.sh
        # Run MaxBin
        """
        }
}
binwork.println()
*/
process  Classify{
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J BG_CLS"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J BG_CLS"
    }
    cpus 10
    memory 10
    errorStrategy 'retry'
    maxErrors 1
    
    publishDir = "${params.taxOutPath}/raw"

    input:
        file fa from krafa
    output:
        file "Taxonomy" into tax
        file("Assigned.fasta")
        file("UnAssigned.fasta")
        file("assign.log") 

    script:
        """
        kraken --threads 5 --fasta-input ${fa} --output output.kraken -db $params.Krakendb > assign.log 2>&1
        kraken-translate --db $params.Krakendb --mpa-format output.kraken >> Taxonomy 
        fish_in_winter.pl          -bf table -ff fasta --bcolumn 1 --fcolumn 1 -gene Taxonomy ${fa} > Assigned.fasta
        fish_in_winter.pl --except -bf table -ff fasta --bcolumn 1 --fcolumn 1 -gene Taxonomy ${fa} > UnAssigned.fasta
        # --except :: get things not in the bait file
        """
}
