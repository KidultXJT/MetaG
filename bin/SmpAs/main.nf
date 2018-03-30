#! /usr/bin/env nextflow

source = "$workflow.projectDir"
dir = "$workflow.launchDir"

/* ------------ Input parameters ------------ */
params.inPath  = "00.MetaGPrp/"
params.outDir = "02.SmpAssembly"
params.smp     = false
params.ident   = ""

/* ------------ Assembly ------------ */
params.asSoft = "megahit"
params.minlen = 1000
params.mincnt = 2
params.kmer   = "21,29,39,59,79,99,119,141"

/* ------------ Cluster ------------ */
params.clstSoft = "cd-hit"
params.clst     = true
params.Memory   = 100000
params.c        = 0.9
params.LenWord  = 5
params.aL       = 0
params.aS       = 0.9
params.Global   = 0

/* ------------ Slurm params ------------ */
params.help   = false
params.local  = false
params.fN     = false
params.xN     = false

/* ------------- HELP Infos ------------ */
if (params.help) {
    log.info "                                                                                                                                                "
    log.info "=============================================                                                                                                   "
    log.info "  MetaGenome Denovo Assembly (Version 0.1)                                                                                                      "
    log.info "=============================================                                                                                                   "
    log.info "                                                                                                                                                "
    log.info "USAGE:                                                                                                                                          "
    log.info "                                                                                                                                                "
    log.info "nextflow run main.nf                                                                                                                            "
    log.info " --inPath     input Direction  [00.MetaGPrp]                                                                                                    "
    log.info " --outPath    output Direction [02.SmpAssembly]                                                                                                 "
    log.info " --smp        input Sample Name (Prefix) [None]                                                                                                 "
    log.info " --ident      indentify  [None]                                                                                                                 "
    log.info " --minlen     minimum length of contigs to output [1000]                                                                                        "
    log.info " --mincnt     minimum multiplicity for filtering (k_min+1)-mers [2]                                                                             "
    log.info " --kmer       comma-separated list of kmer size all must be odd, in the range 15-255, increment <= 28) [21,29,39,59,79,99,119,141]              "
    log.info "                                                                                                                                                "
    log.info " --clst       cluster or not ? [true]                                                                                                           "
    log.info " --Memory     memory limit (in MB) [100000 == 100G]                                                                                             "
    log.info " --c          similarity [0.9]                                                                                                                  "
    log.info " --LenWord    Word Length [5]                                                                                                                   "
    log.info " --Global     Global(0) or Local(1)  [0]                                                                                                        "
    log.info " --aL         Align Longer Sequence (for Global --G)  [0.8]                                                                                     "
    log.info " --aS         Align Shorter Sequence (for Global --G) [0.9]                                                                                     "
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
    log.info"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    log.info"x Sample ID for 'smp'                        x"
    log.info"x                                            x"
    log.info"x inPath/{smp}_R{1,2}.fq.gz.${params.ident}  x"
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


/* ------------- Input ------------ */
/* variables */
smplst = params.smp.tokenize(',')

Channel
    .fromFilePairs("${params.inPath}/*_{1,2}.fq.gz${params.ident}")
    .filter {it[0] in smplst}.into {reads;unmapreads}

process asm {
    tag {SMP}
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J DAS_AS"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J DAS_AS"
    }
    cpus 21
    memory 50G
    errorStrategy 'retry'
    maxErrors 1

    input:
        set val(SMP),file(reads) from reads
    output:
        set val(SMP),file("*.fasta") into smpcontigs
        set val(SMP),file("*_Length*")
        set val(SMP),file("*_GC*")
    
    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} --min-count $params.mincnt --k-list $params.kmer --min-contig-len $params.minlen --num-cpu-threads 20 --memory 0.8
    mv ./megahit_out/*.fa ./megahit_out/${SMP}.megahit.contigs.fa
    python ${source}/bin/FADealEr.py ./megahit_out/${SMP}.megahit.contigs.fa ${SMP}.fasta ${SMP}
    # infa outfa seqID(_num)
    """
}

/* ------------- cluster or Not ------------ */
smpcontigs.into {mergecontigs;mergeclustcontigs}

process MergeClst {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J DAS_MC"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J DAS_MC"
    }
    cpus 20
    memory 100G
    errorStrategy 'retry'
    maxErrors 1

    publishDir params.outDir, mode: 'link', overwrite: true
    
    input:
        file SamplesContigs from mergecontigs.map {it[1]}.collect()
    output:
        file("Samples.fasta")
    
    script:
    if (params.clst) { 
    """
    cat *.fasta > Samples.fa
    cd-hit-est -i Samples.fa -c $params.c -n $params.LenWord -G $params.Global -aS $params.aS -g 1 -d 0 -o Samples.fasta -T 20 -M $params.Memory 
    # -c similarity; -G by the length of the alignment; -aS Coverage; -g accurate; -T 20 threads 20 CPUs; -M unlimitted memory
    """
    }else{
    """
    cat *.fasta > Samples.fasta
    """
    }
}
