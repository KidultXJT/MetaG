#! /usr/bin/env nextflow

source = "$workflow.projectDir"
dir = "$workflow.launchDir"

/* ------------ Input parameters ------------ */
params.inPath  = "00.MetaGPrp/"
params.outDir  = "02.MixAssembly"
params.Ref     = false
params.smp     = false
params.ident   = ""

/* ------------ Mapping ------------ */
params.mapSoft = "bowtie2"

/* Index */
params.buildbt2T = 20

/* Unmap */
params.umbt2I  = 0 
params.umbt2X  = 800
params.umbt2L  = 31 
params.umbt2p  = 16
params.umbt2N  = 1 
params.umbt2k  = 1 

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
    log.info " --inPath     input Direction  [None]                                                                                                           "
    log.info " --outPath    output Direction [02.MixAssembly]                                                                                                 "
    log.info " --smp        input Sample Name (Prefix) [None]                                                                                                 "
    log.info " --Ref        INDEX Ref [None]                                                                                                                  "
    log.info " --ident      indentify  [None]                                                                                                                 "
    log.info "                                                                                                                                                "
    log.info " --mapSoft    Find Unmap Reads Software [bowtie2]                                                                                               "
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

if (!params.Ref) {
    log.info""
    log.info"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    log.info"x Find Unmap Ref Sequence(fasta format)      x"
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
    .fromPath(params.Ref)
    .into {RefIndex;RMergeClst}

Channel
    .fromFilePairs("${params.inPath}/*_{1,2}.fq.gz${params.ident}")
    .filter {it[0] in smplst}.into {unmapreads;quantreads}

process INDEX {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J DAS_UnmapINDEX"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J DAS_UnmapINDEX"
    }
    cpus 20
    memory 100G
    errorStrategy 'retry'
    maxErrors 1

    input:
        file Ref from RefIndex
    output:
        file("Ref*") into INDEX
    
    script:
    """
    bowtie2-build --threads $params.buildbt2T --bmaxdivn 80 $Ref Ref
    """
}
unmapreads.combine(INDEX).set {un} 
/*un.println()*/

process unmap {
    tag {SMP}
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J DAS_Unmap"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J DAS_Unmap"
    }
    cpus 20
    memory 40G
    errorStrategy 'retry'
    maxErrors 1
     
    publishDir "UnmapData", mode: 'link', overwrite: true
    
    input:
        set val(SMP),file(reads),file(Ref1),file(Ref2),file(Ref3),file(Ref4),file(Ref5),file(Ref6),file(Ref7) from un
    output:
        set val(SMP),file("*unmapped_*fq.gz") into unmapReads
        set val(SMP),file("*_align.log") into unmaplog

    script:
    """
    mv ${reads[0]} reads1.fq.gz
    mv ${reads[1]} reads2.fq.gz
    bowtie2 --mm --very-fast -q -I 0 -X 800 -L 31 -p 12 -x Ref -N 1 -k 1 -1 reads1.fq.gz -2 reads2.fq.gz -S /dev/null --un-conc-gz ${SMP}_unmapped_%.fq.gz > ${SMP}_align.log 2>&1
    """
}

/* ------ compute_R830 !!!! Only !!!! ------ */
process mixasm {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "-p compute_R830 -x $params.xN -J DAS_MAS"
    }
    else {
        executor "slurm"
        clusterOptions "-p compute_R830 -J DAS_MAS"
    }
    cpus 45
    memory 200G
    errorStrategy 'retry'
    maxErrors 1
    
    input:
        file reads from unmapReads.map {it[1]}.collect()
    output:
        file("Mix.fa") into mixcontigs
   
    script: 
    """
    zcat *unmapped_1.fq.gz | pigz -c > R1.unmap.fq.gz
    zcat *unmapped_2.fq.gz | pigz -c > R2.unmap.fq.gz
    megahit -1 R1.unmap.fq.gz -2 R2.unmap.fq.gz --min-count $params.mincnt --k-list $params.kmer --min-contig-len $params.minlen --num-cpu-threads 45 -m 0.8
    mv ./megahit_out/*.fa ./megahit_out/mix.megahit.contigs.fa
    python ${source}/bin/FADealEr.py ./megahit_out/mix.megahit.contigs.fa Mix.fa Mix
    # infa outfa seqID(_num)
    """
}

/* ------------- cluster or Not ------------ */
mixcontigs.into {mergecontigs;mergeclustcontigs}

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
        file MixContigs from mergecontigs
    output:
        file("Mix.fa")
        file("Mix.fasta")
        file("Mix.fasta.clstr")
        file("*_Length*")
        file("*_GC*")

    script:
    if (params.clst) { 
    """
    cd-hit-est -i Mix.fa -c $params.c -n $params.LenWord -G $params.Global -aS $params.aS -aL $params.aL -g 1 -d 0 -o Mix.fasta -T 20 -M $params.Memory 
    # -c similarity; -G by the length of the alignment; -aS Coverage; -g accurate; -T 20 threads 20 CPUs; -M unlimitted memory
    python ${source}/bin/FADealEr.py Mix.fa Mix.fasta Mix
    """
    }else{
    """
    python ${source}/bin/FADealEr.py Mix.fa Mix.fasta Mix
    """
    }
}
