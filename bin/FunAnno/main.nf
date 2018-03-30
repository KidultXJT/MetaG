#! /usr/bin/env nextflow

/* define constant */
source = "$workflow.projectDir"
dir = "$workflow.launchDir"

/* define database path */
nrDBPath     = "/Bio/User/yinsijun/database/anno_database/NR/current/"
keggDBPath   = "/Bio/User/yinsijun/database/anno_database/KEGG/current/"
goDBPath     = "/Bio/User/yinsijun/database/anno_database/GO/current/"
cogDBPath    = "/Bio/Database/eggNOG-mapper/"
spDBPath     = "/Bio/User/xiejunting/DataBase/SwissProt/"
cazyDBPath   = "/Bio/User/yinsijun/database/FuncAnno/dbCAN/current/"
ardbDBPath   = "/Bio/User/xiejunting/DataBase/ARDB/"
vfdbDBPath   = "/Bio/User/xiejunting/DataBase/VFDB/"
phiDBPath    = "/Bio/User/xiejunting/DataBase/PHI/"
pfamDBPath   = "/Bio/User/xiejunting/DataBase/PFam/"
phastDBPath  = "/Bio/User/xiejunting/DataBase/PHAST/"
figfamDBPath = "/Bio/User/xiejunting/InterPro/figFam/"

/*spDBPath   = "/Bio/User/yinsijun/database/anno_database/SwissProt/current"*/
/*ardbDBPath = "/Bio/User/yinsijun/database/FuncAnno/ARDB/for_pipe/"*/
/*vfdbDBPath = "/Bio/User/yinsijun/database/FuncAnno/VFDB/for_pipe/"*/

/* input parameters */
params.query          = false
params.pjectID        = "SG"
params.seqtype        = "aa"
params.BasicOutPath   = "06.BasicFun"
params.EnhanceOutPath = "07.EnhanceFun"
params.chunkSize      = 20000

/* BLAST Paramsment */
params.evalue         = "1e-6"
params.identity       = 50
params.maxTarget      = 1

/* DataBase */
params.KEGG    = false
params.KGtax   = "Meta"
params.GO      = false
params.GOtax   = "Meta"
params.SP      = false
params.NR      = false
params.NRtax   = "Meta"
params.COG     = false
params.COGtax  = "bact"
params.COGmode = "diamond"
params.CAZy    = false
params.ARDB    = false
params.VFDB    = false
params.VFDBset = "SetA"
params.PHI     = false
params.PFam    = false
params.PFamset = "SetA"
params.PHAST   = false
params.FIGFam  = false

/* SLURM */
params.xN       = false
params.fN       = false
params.local    = false
params.help     = false
params.version  = false

/* help info */
if (params.help) {
    log.info " "
    log.info "===================================================================="
    log.info " Ez Annotation pipeline (Version 0.5)                               "
    log.info "===================================================================="
    log.info "                                                                    "
    log.info " USAGE:                                                             "
    log.info "                                                                    "
    log.info "nextflow run main.nf \\                                             "
    log.info " --pjectID   Project ID, Example: SGM0000. [SG] \\"
    log.info " --query     Input fasta file. [None] \\                            "
    log.info " --seqtype   Sequence type: aa/nt. [aa] \\                          "
    log.info " --evalue    Evalue threshold. [1e-6] \\                            "
    log.info " --identity  BLAST identity threshold. [50] \\                      "
    log.info " --chunkSize Seq number per chunk. [10000] \\                       "
    log.info "                                                                    "
    log.info " --BasicoutPath     BasicFun Output path. [06.BasicFun] \\          "
    log.info " --EnhanceoutPath   EnhanceFun Output path. [07.EnhanceFun] \\      "
    log.info "                                                                    "
    log.info " ## SLURM ##                                                        "
    log.info " --fN        Force task running on SLURM fat node.overrides --xN [false] \\"
    log.info " --xN        Forbid task running on these SLURM nodes. [false] \\"
    log.info " --local     Force task running on local computer.overrides --xN/--fN [false] \\"
    log.info " --version   Show all databases version. [false] \\                 "
    log.info "                                                                    "
    log.info " ## DataBase ##                                                     "
    log.info " --KEGG      Annotate query with KEGG database. [false] \\          "
    log.info " --KGtax     Specify taxonomy type for KEGG annotation:  \\         "
    log.info "             All/Animal/Archaea/Bacteria/Fungi/Meta/Micro/Plants/Protists.                     [Meta] \\"
    log.info " --GO        Annote query with GO database. [false] \\              "
    log.info " --GOtax     Specify taxonomy type for GO annotation: \\            "
    log.info "             Eukaryota/Fungi/Bacteria/Archaea/Metazoa/Viridiplantae/Viruses/Meta/All.          [Meta] \\"
    log.info " --SP        Annotate query with SwissProt database. [false] \\     "
    log.info " --NR        Annotate query with NR database. [false] \\            "
    log.info " --NRtax     Specify taxonomy type for NR annotation: \\            "
    log.info "             Animals/Archaea/Bacteria/Fungi/Invertebrate/Meta/Microbe/Plants/Vertebrate/Virus. [Meta] \\"
    log.info " --COG       Annotate query with EggNOG database. [false] \\        "
    log.info " --COGtax    Specify taxonomy type for COG annotation: \\           "
    log.info "             euk/bact/arch or path to local hmmpressed database.                               [bact] \\"
    log.info " --COGmode   Method used for search database: hmmer/diamond. [diamond] \\"
    log.info " --CAZy      Annotate query with [CAZy]   database. [false] \\          "
    log.info " --ARDB      Annotate query with [ARDB]   database. [false] \\          "
    log.info " --VFDB      Annotate query with [VFDB]   database. [false] \\          "
    log.info " --PHI       Annotate query with [PHI]    database. [false] \\          "
    log.info " --PFam      Annotate query with [PFam]   database. [false] \\          "
    log.info " --PHAST     Annotate query with [PHAST]  database. [false] \\          "
    log.info " --FIGFam    Annotate query with [FIGFam] database. [false] \\          "
    exit 1
}


if (!params.query) {
    log.info " "
    log.info "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    log.info "A query file/path is necessary! check --help for more information."
    log.info "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    exit 1
}

/* ------------------ Input ----------------- */
Channel
.fromPath(params.query)
.splitFasta(by: params.chunkSize)
.into {Qcog;Qnr}

Channel
.fromPath(params.query)
.into {Qgo;Qsp;Qkegg;Qcazy;Qardb;Qvfdb;Qphi;Qpfam;Qphast;Qfigfam}

process KEGG {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J KEGG"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J KEGG"
    }
    cpus 10
    memory 10G
    errorStrategy 'retry'
    maxErrors 2
    publishDir "$params.BasicOutPath"

    input:
	file(query) from Qkegg
    output:
	file "KEGG/"

    when:
    	params.KEGG

    script:
        /* ---------------- Select Mode ---------------- */
    	mode = (params.seqtype == 'aa')?"blastp":"blastx"
    	db = "${keggDBPath}/fasta/${params.KGtax}.dmnd"
    	"""
    	diamond $mode -d $db --evalue $params.evalue -p 8 --max-target-seqs $params.maxTarget -- --id $params.identity -q $query -o $params.pjectID \
    	--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle
    	python ${source}/bin/keggAnnot.py $params.pjectID $keggDBPath
    	mkdir -p KEGG
        # Rscript
        cat ${params.pjectID}.KEGG.Path.xls | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk -F '\/' '{print $1}' > ${params.pjectID}.KEGG.Path.txt
        Rscript ${source}/bin/keggfreq.r ${params.pjectID}.KEGG.Path.txt ./ ${params.pjectID}.KEGG 
        # avg[1]=intable avg[2]=outDir avg[3]=prefix
        # Krona
        cat ${params.pjectID}.KEGG.Path.xls |awk -F "\t" '{print $1"\t"$2"\t"$3}' |sed '/KEGG_A_class\tKEGG_B_class\tPathway/d' > ${params.pjectID}.KEGG.Path.Text
        perl /Bio/User/xiejunting/Software/Krona/KronaTools/scripts/ImportText.pl ${params.pjectID}.KEGG.Path.Text -o ${params.pjectID}.KEGG.Path.krona.HTML -q
    	mv ${params.pjectID}.KEGG* KEGG
    	mv ${params.pjectID}_map/ KEGG
    	"""
}


process GO {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J GO"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J GO"
    }
    cpus 10
    memory 5G
    errorStrategy 'retry'
    maxErrors 2
    publishDir "$params.BasicOutPath"

    input:
    	file(query) from Qgo

    output:
    	file "GO/"

    when:
    	params.GO

    script:
 	/* --------------- Select Mode ----------------- */ 
    	mode = (params.seqtype == 'aa')?"blastp":"blastx"
    	db = "${goDBPath}/fasta/${params.GOtax}.dmnd"
    	obo = "${goDBPath}/go-basic.obo"
    	"""
    	diamond $mode -d $db --evalue $params.evalue -p 10 --max-target-seqs $params.maxTarget --id $params.identity -q $query -o ${params.pjectID}.out \
    	--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle
    	python ${source}/bin/goAnnot.py ${params.pjectID}.out $obo
    	python ${source}/bin/goLst.py $query
        cat *.glst > ${params.pjectID}.glst
    	/usr/bin/perl ${source}/bin/EnrichGO.pl -fg ${params.pjectID}.glst -bg ${params.pjectID}.glst -a ${params.pjectID} -op ${params.pjectID}.GO -ud nodiff
    	mkdir -p GO
    	mv ${params.pjectID}.GO* GO/
    	"""
}


process NRblast {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''}"
    }
    cpus 10
    memory 10G
    maxForks 5
    errorStrategy 'retry'
    maxErrors 2

    input:
    	file(query) from Qnr

    output:
    	file "NR.tbl" into NRblast

    when:
    	params.NR

    script:
        /* --------------- Select Mode --------------- */
    	mode = (params.seqtype == 'aa')?"blastp":"blastx"
    	db = "${nrDBPath}/${params.NRtax}.fa.dmnd"
    	"""
        diamond $mode -d $db -p 8 --evalue $params.evalue --max-target-seqs $params.maxTarget --id $params.identity -q $query -o NR.tbl \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle
    	echo "diamond $mode -d $db -q $query -o NR.tbl " > log
    	"""
}
NRmerge = NRblast.collectFile()

process NRanno {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J NR_fm"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J NR_fm"
    }
    cpus 1
    memory 1G
    errorStrategy 'retry'
    maxErrors 2
    publishDir "$params.BasicOutPath"

    input:
        file "nr" from NRmerge
    output:
        file "NR/"
    script:
    	"""
    	echo 'query	subject	identity(%)	align_len	mismatch	gap	start	end	subStart	subEnd	Evalue	score	description' | cat - nr > ${params.pjectID}.NR.Anno.xls
    	mkdir NR
    	cp ${params.pjectID}* NR
    	"""
}


process EggNOG {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''}"
    }
    cpus 10
    memory 15G

    input:
    	file(query) from Qcog

    output:
    	file "rmhead" into COGhmm

    when:
    	params.COG

    script:
    	"""
    	emapper.py -i $query -o out --cpu 10 -d $params.COGtax --seed_ortholog_evalue $params.evalue -m $params.COGmode
    	grep ^[^#] out.emapper.annotations > rmhead
    	"""
}
COGmerge = COGhmm.collectFile()

process COGanno {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J COG"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J COG"
    }
    cpus 1
    memory 1G
    errorStrategy 'retry'
    maxErrors 2
    publishDir "$params.BasicOutPath"

    input:
        file "${params.pjectID}" from COGmerge
    output:
        file "COG/"
    script:
    	if (params.COGmode=='hmmer') {
    		"""
    		python ${source}/bin/EggNOG_formatter.py ${params.pjectID}
              	mkdir COG
    		cp ${params.pjectID}.COG.pdf ${params.pjectID}.COG.png COG/
                cp ${params.pjectID}.COG.Anno.xls COG/${params.pjectID}.COG.Anno.xls
                cp ${params.pjectID} COG/${params.pjectID}.EggNOG.Anno.xls
    		"""
    	}else {
    		"""
   		python ${source}/bin/COGdiamond.py ${params.pjectID} ${cogDBPath}/eggnog.db
    		mkdir COG
    		cp ${params.pjectID}.* COG/
    		"""
        }
}


process SwissProt {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J SwProt"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J SwProt"
    }
    cpus 8
    memory 10G
    errorStrategy 'retry'
    maxErrors 2
    publishDir "$params.EnhanceOutPath"

    input:
    	file(query) from Qsp

    output:
    	file "SwissProt/"

    when:
    	params.SP

    script:
        /* ------------- Select Mode ------------- */
    	mode = (params.seqtype == 'aa')?"blastp":"blastx"
    	db = "${spDBPath}/uniprot_sprot.dmnd"
       	"""
    	diamond $mode -d $db --evalue $params.evalue -p 8 --max-target-seqs $params.maxTarget --id $params.identity -q $query -o ${params.pjectID}.out \
    	--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle
    	python ${source}/bin/spFormatter.py ${params.pjectID}.out
    	mkdir SwissProt/
    	mv ${params.pjectID}.* SwissProt/
        mv SwissProt/${params.pjectID}.out ./
    	"""
}

process CAZy {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J CAZy"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J CAZy"
    }
    cpus 10
    memory 10G
    errorStrategy 'retry'
    maxErrors 2
    publishDir "$params.EnhanceOutPath"

    input:
        file fa from Qcazy     
    output:
	file "CAZy/"

    when:
    	params.CAZy

    script:
        /* ---------------- Select Mode ---------------- */
    	db = "${cazyDBPath}/"
    	"""
        hmmsearch --cpu 8 -o out --tblout tsv ${db}/dbCAN-fam-HMMs.txt $fa
        grep ^[^#] tsv > ${params.pjectID}
        python ${source}/bin/CAZy_formatter.py ${params.pjectID} ${db}/FamInfo.txt $params.evalue
        mkdir CAZy
        cp ${params.pjectID}.* CAZy
    	"""
}

process ARDB {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J ARDB"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J ARDB"
    }
    cpus 10
    memory 10G
    errorStrategy 'retry'
    maxErrors 2
    publishDir "$params.EnhanceOutPath"

    input:
        file query from Qardb
    output:
	file "ARDB/"

    when:
    	params.ARDB

    script:
    	mode = (params.seqtype == 'aa')?"blastp":"blastx"
    	db = "${ardbDBPath}/resisGenes.pfasta.dmnd"
       	"""
    	diamond $mode -d $db --evalue $params.evalue -p 8 --max-target-seqs $params.maxTarget --id $params.identity -q $query -o ${params.pjectID} \
    	--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore #stitle
        python ${source}/bin/ARDB_formatter.py ${params.pjectID} ${source}/bin/ARDB.mapping $params.evalue
        mkdir ARDB
        mv ${params.pjectID}.* ARDB
    	"""
}

process VFDB {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J VFDB"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J VFDB"
    }
    cpus 10
    memory 5G
    errorStrategy 'retry'
    maxErrors 2
    publishDir "$params.EnhanceOutPath"

    input:
        file query from Qvfdb
    output:
	file "VFDB/"

    when:
    	params.VFDB

    script:
    	mode = (params.seqtype == 'aa')?"blastp":"blastx"
    	db = (params.VFDBset == "SetA")?"${vfdbDBPath}vfdb.A.dmnd":"${vfdbDBPath}vfdb.B.dmnd"
        map = "${vfdbDBPath}/VFGs.map"
       	"""
    	diamond $mode -d $db --evalue $params.evalue -p 8 --max-target-seqs $params.maxTarget --id $params.identity -q $query -o ${params.pjectID} \
    	--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore #stitle
        python ${source}/bin/VFDB_formatter.py ${params.pjectID} $map $params.evalue
        mkdir VFDB
        mv ${params.pjectID}.* VFDB
    	"""
}


process PHI {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J PHI"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J PHI"
    }
    cpus 10
    memory 5G
    errorStrategy 'retry'
    maxErrors 2
    publishDir "$params.EnhanceOutPath"

    input:
        file query from Qphi
    output:
	file "PHI/"

    when:
    	params.PHI

    script:
    	mode = (params.seqtype == 'aa')?"blastp":"blastx"
    	db = "${phiDBPath}/phi_accessions.dmnd"
       	"""
    	diamond $mode -d $db --evalue $params.evalue -p 8 --max-target-seqs $params.maxTarget --id $params.identity -q $query -o ${params.pjectID} \
    	--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        python ${source}/bin/PHI_formatter.py ${params.pjectID} $params.evalue $params.identity
        mkdir PHI
        mv ${params.pjectID}.* PHI
    	"""
}


process PFam {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J PFam"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J PFam"
    }
    cpus 10
    memory 5G
    errorStrategy 'retry'
    maxErrors 2
    publishDir "$params.EnhanceOutPath"

    input:
        file query from Qpfam
    output:
	file "PFam/"

    when:
    	params.PFam

    script:
    	mode = (params.seqtype == 'aa')?"blastp":"blastx"
        map  = "${pfamDBPath}/pdb_pfam_mapping.txt"
    	db = (params.PFamset == "SetA")?"${pfamDBPath}/Pfam-A.hmm":"${pfamDBPath}/Pfam-B.hmm"
       	"""
        hmmsearch --cpu 8 -o ${params.pjectID}.out --tblout ${params.pjectID}.tsv $db $query
        python ${source}/bin/PFAM_formatter.py ${params.pjectID}.tsv $map $params.evalue 
        mkdir PFam
        mv ${params.pjectID}.* PFam
    	"""
}


process PHAST {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J PHAST"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J PHAST"
    }
    cpus 10
    memory 10G
    errorStrategy 'retry'
    maxErrors 2
    publishDir "$params.EnhanceOutPath"

    input:
        file query from Qphast
    output:
	file "PHAST/"

    when:
    	params.PHAST

    script:
    	mode = (params.seqtype == 'aa')?"blastp":"blastx"
    	db = "${phastDBPath}/prophage_virus.dmnd"
       	"""
    	diamond $mode -d $db --evalue $params.evalue -p 8 --max-target-seqs $params.maxTarget --id $params.identity -q $query -o ${params.pjectID}.blastout \
    	--outfmt 5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle
        python ${source}/bin/PHAST_formatter.py ${params.pjectID}.blastout $params.evalue
        awk \'FNR==1 { header = \$0; print } \$0 != header\' ${params.pjectID}.tbl > ${params.pjectID}.PHAST.xls
        mkdir PHAST
        mv ${params.pjectID}.PHAST.* PHAST
    	"""
}


process FIGfam {
    if (params.local) {
        executor "local"
    }
    else if (params.xN) {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -x $params.xN -J FIGfam"
    }
    else {
        executor "slurm"
        clusterOptions "${params.fN?'-p compute_R830':''} -J FIGfam"
    }
    cpus 10
    memory 10G
    errorStrategy 'retry'
    maxErrors 2
    publishDir "$params.EnhanceOutPath"

    input:
        file query from Qfigfam
    output:
	file "SEED/"

    when:
    	params.FIGFam

    script:
    	mode = (params.seqtype == 'aa')?"blastp":"blastx"
    	db = "${figfamDBPath}/patric.fams/families.dmnd"
       	"""
    	diamond $mode -d $db --evalue $params.evalue -p 8 --max-target-seqs $params.maxTarget --id $params.identity -q $query -o ${params.pjectID} \
    	--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle
        mkdir SEED
        mv ${params.pjectID} SEED
    	"""
}


