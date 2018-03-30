#! /usr/bin/env nextflow

source = "$workflow.projectDir"
dir = "$workflow.launchDir"

/* ------------ Input parameters ------------ */
params.inPath  = "01.CleanData"
params.outPath = "00.MetaGPrp"
params.smp     = false

/* ------------ Assembly params ------------ */
params.trim   = true
params.trimR  = 8
params.trimL  = 12
params.ext    = false
params.extcnt = 10000000000
params.insert = 130

/* ------------ Slurm params ------------ */
params.help   = false
params.local  = false
params.fN     = false
params.xN     = false

/* ------------- HELP Infos ------------ */
if (params.help) {
    log.info "                                                                                                                                                "
    log.info "================================================                                                                                                "
    log.info "  MetaGenome Denovo Assembly Prp (Version 0.1)                                                                                                  "
    log.info "================================================                                                                                                "
    log.info "                                                                                                                                                "
    log.info "USAGE:                                                                                                                                          "
    log.info "                                                                                                                                                "
    log.info "nextflow run main.nf                                                                                                                            "
    log.info " --inPath     input Direction   [01.CleanData]                                                                                                  "
    log.info " --outPath    output Direction  [00.MetaGPrp]                                                                                                   "
    log.info " --smp        input Sample Name (Prefix) [None]                                                                                                 "
    log.info " --trimL      trim Sequence LEFT [12]                                                                                                           "
    log.info " --trimR      trim Sequence RIGHT [8]                                                                                                           "
    log.info " --ext        extract Sequence to Assembly [false]                                                                                              "
    log.info " --extcnt     extract Nums of Sequence to Assembly [10000000000]                                                                                "
    log.info " --insert     [130=150-12-8]                                                                                              "
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
    exit 1
}

/* -------------- Necessary Parameter ------------- */
if (!params.smp) {
    log.info""
    log.info"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    log.info"x Sample ID for 'smp'                   x"
    log.info"x                                       x"
    log.info"x inPath/{smp}_R{1,2}.fq.gz             x"
    log.info"x                                       x"
    log.info"x Example:                              x"
    log.info"x ~/A_R1.fq.gz                          x"
    log.info"x smp = 'A' HERE !!!                    x"
    log.info"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    exit 1
}

if (!params.inPath) {
    log.info""
    log.info"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    log.info"x Your Input File Path for 'inPath'     x"
    log.info"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    exit 1
}

/* ------------- Input ------------ */
/* variables */
smplst = params.smp.tokenize(',')

Channel
    .fromFilePairs("${params.inPath}/*_{1,2}.fq.gz")
    .filter {it[0] in smplst}.into {reads;Reads}

process trim {
    tag {SMP}
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J Trim"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J Trim"
    }
    cpus 21
    memory 40G
    errorStrategy 'retry'
    maxErrors 1
    
    publishDir params.outDir, mode: 'link', overwrite: true
    
    input:
        set val(SMP),file(reads) from reads 
    output:
        set val(SMP),file("*_1.fq.gz.trim") into trimReads1
        set val(SMP),file("*_2.fq.gz.trim") into trimReads2
    
    script:
    """
    seqtk trimfq -b $params.trimL -e $params.trimR ${reads[0]} | pigz -c > ${reads[0]}.trim
    seqtk trimfq -b $params.trimL -e $params.trimR ${reads[1]} | pigz -c > ${reads[1]}.trim
    """
}
trimReads1.into {extReads1;renameReads1}
trimReads2.into {extReads2;renameReads2}

process extract {
    tag {SMP}
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J Extract"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J Extract"
    }
    cpus 21
    memory 40G
    errorStrategy 'retry'
    maxErrors 1
    
    publishDir params.outDir , mode: 'link', overwrite: true
    
    input:
        set val(SMP),file(reads1) from extReads1
        set val(SMP),file(reads2) from extReads2
    output:
        set val(SMP),file("*.ext")
    
    when:
        params.ext

    script:
    Cont    = params.extcnt/2/params.insert
    extCont = (int)Cont
    """ 
    seqtk sample -s100 ${reads1} $extCont | pigz -c > ${reads1}.ext
    seqtk sample -s100 ${reads2} $extCont | pigz -c > ${reads2}.ext
    """
}
