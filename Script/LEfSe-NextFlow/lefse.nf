#!/usr/bin/env nextflow


/******* Parameters *******/

/***** Input and Output *****/
params.pfx = "0.001"
params.inputfile = "/Bio/User/yangliyan/test/MetaSSU/03.Analysis/01.OTUTable/0.001/0.001_otu_table.biom"
params.outdir = "03.Analysis/04.Lefse/"
params.mapfile = "03.Analysis/mapping_file.txt"

/***** LEfSe *****/
params.level = 7
params.cl = "BarcodeSequence"
params.split = "LinkerPrimerSequence"
params.p = 0.05
params.score = 2
params.strict = 0

/***** Plot *****/
params.dpi = 100

params.silent = false

log.info "========= [MetaSSU] LefSe Analysis ========="
log.info ""
log.info "The parameters are:"
log.info ""
log.info "==> Analysis Prefix:\t\t$params.pfx"
log.info "==> Input File:\t\t\t$params.inputFile"
log.info "==> Output Direcotry:\t\t$params.outDir"
log.info "==> Mapping File:\t\t$params.mapFile"
log.info "==> Level:\t\t\t$params.level"
log.info "==> p-value threashold:\t\t$params.p"
log.info "==> LDA score threashold:\t$params.score"
log.info "==> More strict:\t\t$params.strict"
log.info ""


inputOTU = Channel.fromPath(params.inputFile)
inputMAP = Channel.fromPath(params.mapFile)
Channel.from(params.pfx).into{lefseOut;plotOut}
level = Channel.from(params.level)
cl = Channel.from(params.cl)
split = Channel.from(params.split)
p = Channel.from(params.p)
score = Channel.from(params.score)
strict = Channel.from(params.strict)

process lefse{

	input:
	file inputOTU
	file inputMAP
	val level
	val cl
	val split
	val p
	val score
	val strict
	val pfx from lefseOut

	publishDir "$params.outDir"

	output:
	file "${pfx}_lefse_output.txt.formatted.forplot" into plotInput
	file "${pfx}_lefse_output.txt"

	"""
	koeken.py -i $inputOTU -o "./" -m $inputMAP -l $level -cl $cl --split $split -p $p -e $score -str $strict
	mv lefse_output/run_lefse/*.txt ${pfx}_lefse_output.txt
	lefse_format.py ${pfx}_lefse_output.txt
	rm -rf ${pfx}_lefse_output.txt
	mv ${pfx}_lefse_output.txt.formatted ${pfx}_lefse_output.txt
	"""

}

process plot{

	input:
	file plotInput
	val pfx from plotOut

	output:
	file "${pfx}_LEfSe_barplot.pdf"
	file "${pfx}_LEfSe_barplot.png"
	file "${pfx}_LEfSe_cladogram.pdf"
	file "${pfx}_LEfSe_cladogram.png"

	"""
	plot_res.py $plotInput ${pfx}_LEfSe_barplot.pdf --dpi 100 --format pdf
	plot_res.py $plotInput ${pfx}_LEfSe_barplot.png --dpi 100 --format png
	plot_cladogram.py $plotInput ${pfx}_LEfSe_cladogram.pdf --dpi 100 --format pdf
	plot_cladogram.py $plotInput ${pfx}_LEfSe_cladogram.png --dpi 100 --format png
	"""

}

workflow.onComplete{
	log.info ""
	log.info "All done."
	log.info ""
	log.info "========= [MetaSSU] LefSe Analysis ========="
}

workflow.onError{
	log.info "Workflow execution stopped with the following message:"
	log.info " " + workflow.errorMessage
}
