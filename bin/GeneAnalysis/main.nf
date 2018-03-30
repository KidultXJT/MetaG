#! /usr/bin/env nextflow

source = "$workflow.projectDir"
dir = "$workflow.launchDir"

/* ------------ Input parameters ------------ */
params.tax        = "Taxonomy"
params.smp        = false
params.grp        = false
params.inPath     = false
params.OutPath    = "03.GeneAnalysis"

params.identcov   = ".cov.txt"
params.prefixcov  = "Abundance"
params.valuescov  = 1
params.keycov     = 0
params.headercov  = "None"

params.identpkm   = ".pkm.txt"
params.prefixpkm  = "Reads"
params.valuespkm  = 1
params.keypkm     = 0
params.headerpkm  = "None"

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
    log.info " --smp        Samples list. Example A1,A2,A3,B1,B2,B3,C1,C2,C3 [None]                                                                           "
    log.info " --grp        Groups  list. Example A,A,A,B,B,B,C,C,C          [None]                                                                           "
    log.info " --inPath     Quant files :: A1.cov.txt A2.cov.txt...          [None]                                                                           "
    log.info " --OutPath    Quant Output Path                                [03.GeneAnalysis]                                                                "
    log.info "                                                                                                                                                "
    log.info " --identcov   indentify(by Default Can Use :: '.cov.txt')  [.cov.txt]                                                                           "
    log.info " --valuescov  values Column [1 :: AvgDepth]                                                                                                     "
    log.info " --keycov     Key Column    [0]                                                                                                                 "
    log.info " --headercov  Header or Not [None]                                                                                                              "
    log.info "                                                                                                                                                "
    log.info " --identpkm   indentify(by Default Can Use :: '.pkm.txt')  [.pkm.txt]                                                                           "
    log.info " --valuespkm  values Column [4 :: Reads]                                                                                                        "
    log.info " --keypkm     Key Column    [0]                                                                                                                 "
    log.info " --headerpkm  Header or Not [None]                                                                                                              "
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
if (!params.smp) {
    log.info""
    log.info"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    log.info"x Sample ID for 'smp'                        x"
    log.info"x                                            x"
    log.info"x inPath/{smp}${params.identcov}             x"
    log.info"x inPath/{smp}${params.identpkm}             x"
    log.info"x                                            x"
    log.info"x Example:                                   x"
    log.info"x ~/A.txt                                    x"
    log.info"x smp = 'A' HERE !!!                         x"
    log.info"x ident = '.txt'                             x"
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

/* ---------------------------------- Input ------------------------------------- */
/* variables */
smplst = params.smp.tokenize(',')
grplst = params.grp.tokenize(',')

Channel
    .fromPath("${params.inPath}/*${params.identcov}")
    .filter {it.baseName.tokenize('.')[0] in smplst}.into {covExp;covLEfSe}

Channel
    .fromPath("${params.inPath}/*${params.identpkm}")
    .filter {it.baseName.tokenize('.')[0] in smplst}.into {pkmExp;pkmLEfSe}

lefsetaxcovfiles=smplst.join("${params.identcov},")+"${params.identcov}"
process  lefseCov{
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J LEfSe_COV"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J LEfSe_COV"
    }
    cpus 1
    memory 5
    errorStrategy 'retry'
    maxErrors 1
    
    publishDir "${params.OutPath}", mode: 'link', overwrite: true

    input:
        file lefsecovs from covLEfSe.collect()
        file lefsetax  from lefseTaxcov 
    output:
        file("LEfSe_${params.prefixcov}")

    script:
    """
    mkdir LEfSe_${params.prefixcov}
    # ------------------- LEfSe Prp >>>>>>>>>>>>
    python ${source}/bin/lefsefrom.py $lefsetax $lefsetaxcovfiles "${params.delimited}" $params.smp $params.grp $params.valuescov $params.keycov 0 LEfSe_${params.prefixcov} ${params.prefixcov}
    # Taxonomy From Kraken
    # quantfiles
    # "|"
    # smp
    # grp
    # cov Value Column
    # key Value Column
    # Header or Not
    # OutDir 
    # Out Prefix
    # -------------------- LEfSe Analysis >>>>>>>>>>>>>
    python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/format_input.py   LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.txt LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.in -c 1 -u 2 # no normalization
    python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/run_lefse.py      LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.in  LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.res
    python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_cladogram.py LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.res LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.cladogram.png --format png --dpi 100
    python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_cladogram.py LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.res LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.cladogram.pdf --format pdf
    python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_res.py       LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.res LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.png --format png --dpi 100
    python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_res.py       LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.res LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.pdf --format pdf
    # -------------------- Make Filter Result >>>>>>>>>>>>>>>
    echo "#cat LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.res | sed -e '/-/d' | sed '/Other/d' | sed '/Ambiguous/d' | sed '/unclassified/d' | sed -e '/\\.\t/d' | sort -k3 -n > ${params.prefixcov}.LEfSe.FLT.res" >> log
    echo "#python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_cladogram.py LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.FLT.res LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.FLT.cladogram.png --format png --dpi 100" >> log
    echo "#python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_cladogram.py LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.FLT.res LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.FLT.cladogram.png --format pdf" >> log
    echo "#python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_res.py       LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.FLT.res LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.FLT.png --format png --dpi 100" >> log
    echo "#python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_res.py       LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.FLT.res LEfSe_${params.prefixcov}/${params.prefixcov}.LEfSe.FLT.png --format pdf" >> log
    """
}

lefsetaxpkmfiles=smplst.join("${params.identpkm},")+"${params.identpkm}"
process  lefsePKM{
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J LEfSe_PKM"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J LEfSe_PKM"
    }
    cpus 1
    memory 5
    errorStrategy 'retry'
    maxErrors 1
    
    publishDir "${params.OutPath}", mode: 'link', overwrite: true

    input:
        file lefsepkms from pkmLEfSe.collect()
        file lefsetax  from lefseTaxpkm    
    output:
        file("LEfSe_${params.prefixcov}")

    script:
    """
    mkdir LEfSe_${params.prefixpkm}
    # ------------------- LEfSe Prp >>>>>>>>>>>>
    python ${source}/bin/lefsefrom.py $lefsetax $lefsetaxpkmfiles "${params.delimited}" $params.smp $params.grp $params.valuespkm $params.keypkm 0 LEfSe_${params.prefixpkm} ${params.prefixpkm}
    # Taxonomy From Kraken
    # quantfiles
    # "|"
    # smp
    # grp
    # cov Value Column
    # key Value Column
    # Header or Not
    # OutDir 
    # Out Prefix
    # -------------------- LEfSe Analysis >>>>>>>>>>>>>
    python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/format_input.py   LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.txt LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.in -c 1 -u 2 # no normalization
    python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/run_lefse.py      LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.in  LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.res
    python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_cladogram.py LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.res LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.cladogram.png --format png --dpi 100
    python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_cladogram.py LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.res LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.cladogram.png --format pdf
    python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_res.py       LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.res LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.png --format png --dpi 100
    python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_res.py       LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.res LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.png --format pdf
    # -------------------- Make Filter Result >>>>>>>>>>>>>>>
    echo "#cat LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.res | sed -e '/-/d' | sed '/Other/d' | sed '/Ambiguous/d' | sed '/unclassified/d' | sed -e '/\\.\t/d' | sort -k3 -n > ${params.prefixpkm}.LEfSe.FLT.res" >> log
    echo "#python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_cladogram.py LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.FLT.res LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.FLT.cladogram.png --format png --dpi 100" >> log
    echo "#python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_cladogram.py LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.FLT.res LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.FLT.cladogram.png --format pdf" >> log
    echo "#python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_res.py       LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.FLT.res LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.FLT.png --format png --dpi 100" >> log
    echo "#python ${source}/bin/LEfSe-NextFlow/nsegata-lefse-82605a2ae7b7/plot_res.py       LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.FLT.res LEfSe_${params.prefixpkm}/${params.prefixpkm}.LEfSe.FLT.png --format pdf" >> log
    """
}

krataxcovfiles=smplst.join("${params.identcov},")+"${params.identcov}"
process  krataxCOV{
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J Kratax_COV"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J Kratax_COV"
    }
    cpus 10
    memory 10
    errorStrategy 'retry'
    maxErrors 1
    
    publishDir "${params.OutPath}", mode: 'link', overwrite: true

    input:
        file covs   from covExp.collect()
        file covtax from covTax
    
    output:
        file("${params.prefixcov}")
        file("Differential_${params.prefixcov}")
    
    script:
    """
    python ${source}/bin/KraTaxQuant.py $covtax "${params.delimited}" $krataxcovfiles $params.smp $params.valuescov $params.keycov 0 ./ ${params.prefixcov} 
    # Rscript
    # Heatmap StackPlot BoxPlot(All Samples)
    Rscript ${source}/bin/phmstackbox.r ${params.prefixcov}/RelativeAbundance/Kingdom.RelativeAbundance.xls ${params.prefixcov}/RelativeAbundance/ Kingdom
    Rscript ${source}/bin/phmstackbox.r ${params.prefixcov}/RelativeAbundance/Phylum.RelativeAbundance.xls  ${params.prefixcov}/RelativeAbundance/ Phylum
    Rscript ${source}/bin/phmstackbox.r ${params.prefixcov}/RelativeAbundance/Class.RelativeAbundance.xls   ${params.prefixcov}/RelativeAbundance/ Class
    Rscript ${source}/bin/phmstackbox.r ${params.prefixcov}/RelativeAbundance/Order.RelativeAbundance.xls   ${params.prefixcov}/RelativeAbundance/ Order
    Rscript ${source}/bin/phmstackbox.r ${params.prefixcov}/RelativeAbundance/Family.RelativeAbundance.xls  ${params.prefixcov}/RelativeAbundance/ Family
    Rscript ${source}/bin/phmstackbox.r ${params.prefixcov}/RelativeAbundance/Genus.RelativeAbundance.xls   ${params.prefixcov}/RelativeAbundance/ Genus
    Rscript ${source}/bin/phmstackbox.r ${params.prefixcov}/RelativeAbundance/Species.RelativeAbundance.xls ${params.prefixcov}/RelativeAbundance/ Species
    # BetaDiv
    # Heatmap PCoA/PCA
    Rscript ${source}/bin/betadiv.r ${params.prefixcov}/RelativeAbundance/Kingdom_Dominant.xls ${params.prefixcov}/RelativeAbundance $params.smp $params.grp ${params.prefixcov}_Kingdom_Dominant
    Rscript ${source}/bin/betadiv.r ${params.prefixcov}/RelativeAbundance/Phylum_Dominant.xls  ${params.prefixcov}/RelativeAbundance $params.smp $params.grp ${params.prefixcov}_Phylum_Dominant
    Rscript ${source}/bin/betadiv.r ${params.prefixcov}/RelativeAbundance/Class_Dominant.xls   ${params.prefixcov}/RelativeAbundance $params.smp $params.grp ${params.prefixcov}_Class_Dominant
    Rscript ${source}/bin/betadiv.r ${params.prefixcov}/RelativeAbundance/Order_Dominant.xls   ${params.prefixcov}/RelativeAbundance $params.smp $params.grp ${params.prefixcov}_Order_Dominant
    Rscript ${source}/bin/betadiv.r ${params.prefixcov}/RelativeAbundance/Family_Dominant.xls  ${params.prefixcov}/RelativeAbundance $params.smp $params.grp ${params.prefixcov}_Family_Dominant
    Rscript ${source}/bin/betadiv.r ${params.prefixcov}/RelativeAbundance/Genus_Dominant.xls   ${params.prefixcov}/RelativeAbundance $params.smp $params.grp ${params.prefixcov}_Genus_Dominant
    Rscript ${source}/bin/betadiv.r ${params.prefixcov}/RelativeAbundance/Species_Dominant.xls ${params.prefixcov}/RelativeAbundance $params.smp $params.grp ${params.prefixcov}_Species_Dominant
    # PCA Euclidean Div
    Rscript ${source}/bin/taxpca.r ${params.prefixcov}/RelativeAbundance/Kingdom_Dominant.xls ${params.prefixcov}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Kingdom_Dominant
    Rscript ${source}/bin/taxpca.r ${params.prefixcov}/RelativeAbundance/Phylum_Dominant.xls  ${params.prefixcov}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Phylum_Dominant
    Rscript ${source}/bin/taxpca.r ${params.prefixcov}/RelativeAbundance/Class_Dominant.xls   ${params.prefixcov}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Class_Dominant
    Rscript ${source}/bin/taxpca.r ${params.prefixcov}/RelativeAbundance/Order_Dominant.xls   ${params.prefixcov}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Order_Dominant
    Rscript ${source}/bin/taxpca.r ${params.prefixcov}/RelativeAbundance/Family_Dominant.xls  ${params.prefixcov}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Family_Dominant
    Rscript ${source}/bin/taxpca.r ${params.prefixcov}/RelativeAbundance/Genus_Dominant.xls   ${params.prefixcov}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Genus_Dominant
    Rscript ${source}/bin/taxpca.r ${params.prefixcov}/RelativeAbundance/Species_Dominant.xls ${params.prefixcov}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Species_Dominant
    # Differential Analysis
    mkdir Differential_${params.prefixcov}
    mkdir Differential_${params.prefixcov}/LC_Report
    python ${source}/bin/LC-Report-master/LCGDF/ezDiff.py -d ${params.prefixcov}/RelativeAbundance/Kingdom_Dominant.xls -t 1 -g $params.grp --rn 1 --gdir Differential_${params.prefixcov}/ --pdir Differential_${params.prefixcov} --ptest -p ${params.prefixcov}_Kingdom_Dominant --keep_rmd --html 
    python ${source}/bin/LC-Report-master/LCGDF/ezDiff.py -d ${params.prefixcov}/RelativeAbundance/Phylum_Dominant.xls  -t 1 -g $params.grp --rn 1 --gdir Differential_${params.prefixcov}/ --pdir Differential_${params.prefixcov} --ptest -p ${params.prefixcov}_Phylum_Dominant  --keep_rmd --html 
    python ${source}/bin/LC-Report-master/LCGDF/ezDiff.py -d ${params.prefixcov}/RelativeAbundance/Class_Dominant.xls   -t 1 -g $params.grp --rn 1 --gdir Differential_${params.prefixcov}/ --pdir Differential_${params.prefixcov} --ptest -p ${params.prefixcov}_Class_Dominant   --keep_rmd --html 
    python ${source}/bin/LC-Report-master/LCGDF/ezDiff.py -d ${params.prefixcov}/RelativeAbundance/Order_Dominant.xls   -t 1 -g $params.grp --rn 1 --gdir Differential_${params.prefixcov}/ --pdir Differential_${params.prefixcov} --ptest -p ${params.prefixcov}_Order_Dominant   --keep_rmd --html 
    python ${source}/bin/LC-Report-master/LCGDF/ezDiff.py -d ${params.prefixcov}/RelativeAbundance/Family_Dominant.xls  -t 1 -g $params.grp --rn 1 --gdir Differential_${params.prefixcov}/ --pdir Differential_${params.prefixcov} --ptest -p ${params.prefixcov}_Family_Dominant  --keep_rmd --html 
    python ${source}/bin/LC-Report-master/LCGDF/ezDiff.py -d ${params.prefixcov}/RelativeAbundance/Genus_Dominant.xls   -t 1 -g $params.grp --rn 1 --gdir Differential_${params.prefixcov}/ --pdir Differential_${params.prefixcov} --ptest -p ${params.prefixcov}_Genus_Dominant   --keep_rmd --html 
    python ${source}/bin/LC-Report-master/LCGDF/ezDiff.py -d ${params.prefixcov}/RelativeAbundance/Species_Dominant.xls -t 1 -g $params.grp --rn 1 --gdir Differential_${params.prefixcov}/ --pdir Differential_${params.prefixcov} --ptest -p ${params.prefixcov}_Species_Dominant --keep_rmd --html 
    # Print Venn 
    """
}

krataxpkmfiles=smplst.join("${params.identpkm},")+"${params.identpkm}"
process  krataxPKM{
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J Kratax_PKM"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J Kratax_PKM"
    }
    cpus 10
    memory 10
    errorStrategy 'retry'
    maxErrors 1
    
    publishDir "${params.OutPath}", mode: 'link', overwrite: true

    input:
        file pkms   from pkmExp.collect()
        file pkmtax from pkmTax
    
    output:
        file("${params.prefixpkm}")
    
    script:
    """
    python ${source}/bin/KraTaxQuant.py $pkmtax "${params.delimited}" $krataxpkmfiles $params.smp $params.valuespkm $params.keypkm 0 ./ ${params.prefixpkm}
    # Rscript
    # Heatmap StackPlot BoxPlot(All Samples)
    Rscript ${source}/bin/phmstackbox.r ${params.prefixpkm}/RelativeAbundance/Kingdom.RelativeAbundance.xls ${params.prefixpkm}/RelativeAbundance/ Kingdom
    Rscript ${source}/bin/phmstackbox.r ${params.prefixpkm}/RelativeAbundance/Phylum.RelativeAbundance.xls  ${params.prefixpkm}/RelativeAbundance/ Phylum
    Rscript ${source}/bin/phmstackbox.r ${params.prefixpkm}/RelativeAbundance/Class.RelativeAbundance.xls   ${params.prefixpkm}/RelativeAbundance/ Class
    Rscript ${source}/bin/phmstackbox.r ${params.prefixpkm}/RelativeAbundance/Order.RelativeAbundance.xls   ${params.prefixpkm}/RelativeAbundance/ Order
    Rscript ${source}/bin/phmstackbox.r ${params.prefixpkm}/RelativeAbundance/Family.RelativeAbundance.xls  ${params.prefixpkm}/RelativeAbundance/ Family
    Rscript ${source}/bin/phmstackbox.r ${params.prefixpkm}/RelativeAbundance/Genus.RelativeAbundance.xls   ${params.prefixpkm}/RelativeAbundance/ Genus
    Rscript ${source}/bin/phmstackbox.r ${params.prefixpkm}/RelativeAbundance/Species.RelativeAbundance.xls ${params.prefixpkm}/RelativeAbundance/ Species
    # BetaDiv
    # Heatmap PCoA/PCA
    # PCoA
    Rscript ${source}/bin/betadiv.r ${params.prefixpkm}/RelativeAbundance/Kingdom_Dominant.xls ${params.prefixpkm}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Kingdom_Dominant
    Rscript ${source}/bin/betadiv.r ${params.prefixpkm}/RelativeAbundance/Phylum_Dominant.xls  ${params.prefixpkm}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Phylum_Dominant
    Rscript ${source}/bin/betadiv.r ${params.prefixpkm}/RelativeAbundance/Class_Dominant.xls   ${params.prefixpkm}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Class_Dominant
    Rscript ${source}/bin/betadiv.r ${params.prefixpkm}/RelativeAbundance/Order_Dominant.xls   ${params.prefixpkm}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Order_Dominant
    Rscript ${source}/bin/betadiv.r ${params.prefixpkm}/RelativeAbundance/Family_Dominant.xls  ${params.prefixpkm}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Family_Dominant
    Rscript ${source}/bin/betadiv.r ${params.prefixpkm}/RelativeAbundance/Genus_Dominant.xls   ${params.prefixpkm}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Genus_Dominant
    Rscript ${source}/bin/betadiv.r ${params.prefixpkm}/RelativeAbundance/Species_Dominant.xls ${params.prefixpkm}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Species_Dominant
    # PCA Euclidean Div
    Rscript ${source}/bin/taxpca.r ${params.prefixpkm}/RelativeAbundance/Kingdom_Dominant.xls ${params.prefixpkm}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Kingdom_Dominant
    Rscript ${source}/bin/taxpca.r ${params.prefixpkm}/RelativeAbundance/Phylum_Dominant.xls  ${params.prefixpkm}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Phylum_Dominant
    Rscript ${source}/bin/taxpca.r ${params.prefixpkm}/RelativeAbundance/Class_Dominant.xls   ${params.prefixpkm}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Class_Dominant
    Rscript ${source}/bin/taxpca.r ${params.prefixpkm}/RelativeAbundance/Order_Dominant.xls   ${params.prefixpkm}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Order_Dominant
    Rscript ${source}/bin/taxpca.r ${params.prefixpkm}/RelativeAbundance/Family_Dominant.xls  ${params.prefixpkm}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Family_Dominant
    Rscript ${source}/bin/taxpca.r ${params.prefixpkm}/RelativeAbundance/Genus_Dominant.xls   ${params.prefixpkm}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Genus_Dominant
    Rscript ${source}/bin/taxpca.r ${params.prefixpkm}/RelativeAbundance/Species_Dominant.xls ${params.prefixpkm}/RelativeAbundance $params.smp $params.grp ${params.prefixpkm}_Species_Dominant
    # Differential Analysis
    mkdir Differential_${params.prefixpkm}
    mkdir Differential_${params.prefixpkm}/LC_Report
    python ${source}/bin/LC-Report-master/LCGDF/ezDiff.py -d ${params.prefixpkm}/RelativeAbundance/Kingdom_Dominant.xls -t 1 -g $params.grp --rn 1 --gdir Differential_${params.prefixpkm}/ --pdir Differential_${params.prefixpkm} --ptest -p ${params.prefixpkm}_Kingdom_Dominant --keep_rmd --html 
    python ${source}/bin/LC-Report-master/LCGDF/ezDiff.py -d ${params.prefixpkm}/RelativeAbundance/Phylum_Dominant.xls  -t 1 -g $params.grp --rn 1 --gdir Differential_${params.prefixpkm}/ --pdir Differential_${params.prefixpkm} --ptest -p ${params.prefixpkm}_Phylum_Dominant  --keep_rmd --html 
    python ${source}/bin/LC-Report-master/LCGDF/ezDiff.py -d ${params.prefixpkm}/RelativeAbundance/Class_Dominant.xls   -t 1 -g $params.grp --rn 1 --gdir Differential_${params.prefixpkm}/ --pdir Differential_${params.prefixpkm} --ptest -p ${params.prefixpkm}_Class_Dominant   --keep_rmd --html 
    python ${source}/bin/LC-Report-master/LCGDF/ezDiff.py -d ${params.prefixpkm}/RelativeAbundance/Order_Dominant.xls   -t 1 -g $params.grp --rn 1 --gdir Differential_${params.prefixpkm}/ --pdir Differential_${params.prefixpkm} --ptest -p ${params.prefixpkm}_Order_Dominant   --keep_rmd --html 
    python ${source}/bin/LC-Report-master/LCGDF/ezDiff.py -d ${params.prefixpkm}/RelativeAbundance/Family_Dominant.xls  -t 1 -g $params.grp --rn 1 --gdir Differential_${params.prefixpkm}/ --pdir Differential_${params.prefixpkm} --ptest -p ${params.prefixpkm}_Family_Dominant  --keep_rmd --html 
    python ${source}/bin/LC-Report-master/LCGDF/ezDiff.py -d ${params.prefixpkm}/RelativeAbundance/Genus_Dominant.xls   -t 1 -g $params.grp --rn 1 --gdir Differential_${params.prefixpkm}/ --pdir Differential_${params.prefixpkm} --ptest -p ${params.prefixpkm}_Genus_Dominant   --keep_rmd --html 
    python ${source}/bin/LC-Report-master/LCGDF/ezDiff.py -d ${params.prefixpkm}/RelativeAbundance/Species_Dominant.xls -t 1 -g $params.grp --rn 1 --gdir Differential_${params.prefixpkm}/ --pdir Differential_${params.prefixpkm} --ptest -p ${params.prefixpkm}_Species_Dominant --keep_rmd --html 
    """
}


