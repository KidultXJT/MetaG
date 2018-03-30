#!perl
use warnings;
use strict;
use File::Basename qw(basename);
use Getopt::Long;
use FileHandle;

###########################################################
use FindBin qw($Bin);
my $drawclass = "$Bin/Draw_GoClass.pl";

###########################################################
my %opts;
GetOptions(\%opts,"fg=s","bg=s","a=s","op=s", "p=f", "q=f", "ud=s");

my $usage = <<"USAGE";
        Usage:  perl $0 -fg <gene.list> -bg <backgroundgene.list> -a <annotprefix> -op <test> [options]
        Notes:  Need R
        Options:
                -fg     files     gene list
                -bg     files     background gene list
                -a      files     annotation prefix (*).[P|F|C]
                -op     prefix    output prefix
                -ud     strings   updown table(2 col) or only gene(1 col), default 'nodiff', or 'diff'
                -p      float     pvalue, default 0.05
                -q      float     qvalue, default 1
        eg:
                perl $0 -fg gene.list -bg backgroundgene.list -a annotprefix -op test
USAGE

die $usage if(!$opts{fg} or !$opts{bg} or !$opts{a} or !$opts{op});
$opts{p}=$opts{p}?$opts{p}:0.05;
$opts{q}=$opts{q}?$opts{q}:1;
$opts{ud}=$opts{ud}?$opts{ud}:"nodiff";

my @afs;
push @afs, "$opts{a}.P" if(-s "$opts{a}.P");
push @afs, "$opts{a}.F" if(-s "$opts{a}.F");
push @afs, "$opts{a}.C" if(-s "$opts{a}.C");

my %namespace = (
	"P" => "Biological Process",
	"F" => "Molecular Function",
	"C" => "Cellular Component"
);

my %gene_ud;
if($opts{ud} eq "diff"){
	open GENE_LIST, "$opts{fg}" || die $!;
	while(<GENE_LIST>){
		chomp;
		my @tmp = split /\t/;
		$gene_ud{$tmp[0]} = $tmp[1];
	}
	close GENE_LIST;
}

my %slgo;
my %level_file;
my $title = basename($opts{op});
for my $i (2..4){
	my $sl = "/Bio/User/yinsijun/database/anno_database/GO/current/GO_level$i.xls";

	open SL, $sl or die $!;
	while(<SL>)
	{
		chomp;
		my @tmp = split /\t/;
		$slgo{$i}{$tmp[1]} = "$tmp[0]\t$tmp[2]";
	}

	open $level_file{$i}, "> $opts{op}.Level$i.xls" or die $!;
	if($opts{ud} eq "nodiff"){
		$level_file{$i}->print("Ontology\tClass\tnumber_of_$title\tgenes_of_$title\n");
	}
	elsif($opts{ud} eq "diff"){
		$level_file{$i}->print("Ontology\tClass\tnumber_of_${title}_up\tnumber_of_${title}_down\tgenes_of_${title}_up\tgenes_of_${title}_down\n");
	}

}

my @xls;
foreach(@afs)
{
	my ($t) = $_ =~ /\.(\S)$/;
	my $type;
	if($t eq "C")
	{
		$type = "CC";
	}elsif($t eq "P"){
		$type = "BP";
	}elsif($t eq "F"){
		$type = "MF";
	}else{
		die $usage;
	}

	open CMD, "> $opts{op}.plot.r" or die $!;
	print CMD "
	library(clusterProfiler)
	annot <- read.table(\"$_\", header=F, sep=\"\\t\", quote=\"\")
	allg <- read.table(\"$opts{bg}\", header=F)
	deg <- read.table(\"$opts{fg}\", header=F)
	g2g <- annot[, c(\"V1\", \"V2\")]
	g2t <- annot[, c(\"V1\", \"V3\")]
	ego <- enricher(deg\$V1, universe=allg\$V1, pAdjustMethod = \"none\", minGSSize = 1, maxGSSize = Inf, TERM2GENE= g2g, TERM2NAME = g2t, pvalueCutoff =1, qvalueCutoff =1)
	write.table(summary(ego)[,-6][,-8], \"$opts{op}.$t.xls\", sep=\"\\t\", row.names=F, quote=F)
	ego <- enricher(deg\$V1, universe=allg\$V1, pAdjustMethod = \"none\", minGSSize = 1, maxGSSize = Inf, TERM2GENE= g2g, TERM2NAME = g2t, pvalueCutoff =0.05, qvalueCutoff =1)
	ego\@ontology = \"$type\"
	if(length(rownames(summary(ego))) != 0)
	{
	pdf(\"$opts{op}.$t.pdf\")
	if(length(rownames(summary(ego))) > 9)
	{
		plotGOgraph(ego)
	}else{
		plotGOgraph(ego, firstSigNodes = length(rownames(summary(ego))))
	}
	dev.off()
	}
	";
	close CMD;
	
	`Rscript $opts{op}.plot.r`;
	`convert -density 300 $opts{op}.$t.pdf $opts{op}.$t.png` if(-s "$opts{op}.$t.pdf");

	open XLS,      "< $opts{op}.$t.xls" || die $!;
	open HEAD_XLS, "> $opts{op}.$t.xls.head" || die $!;
	open OUT_XLS,  "> $opts{op}.$t.xls.tmp" || die $!;
	my ($fg, $bg) = (0, 0);
	my $head = <XLS>; chomp $head;
	my @head = split /\t/, $head;
	while(<XLS>){
		chomp;
		if($_ =~ s#(\d+)/(\d+)\t(\d+)/(\d+)#$1\t$3#){
			$fg = $2;
			$bg = $4;
		}
		my @tmp = split /\t/;
		my @modify = split /\//, $tmp[-1];
		$tmp[-1] = join ";", @modify;
		$tmp[4] = 0 if($tmp[4] < 1e-308);
		$tmp[5] = 0 if($tmp[5] < 1e-308);
		my $txt = join "\t", @tmp;
		print OUT_XLS "$txt\n";
	}
	$head[0] = "GO ID";
	$head[2] .= " ($fg)";
	$head[3] .= " ($bg)";
	$head = join "\t", @head;
	print HEAD_XLS "$head\n";
	close XLS;
	close HEAD_XLS;
	close OUT_XLS;
	`cat $opts{op}.$t.xls.head $opts{op}.$t.xls.tmp > $opts{op}.$t.xls`;
	`rm $opts{op}.$t.xls.head $opts{op}.$t.xls.tmp -rf`;

#my %ud;
	open PFC, "$opts{op}.$t.xls" or die $!;
	<PFC>;
	if($opts{ud} eq "nodiff")
	{
		while(<PFC>)
		{
			chomp;
			my @tmp = split /\t/;
			for my $i (2..4){
				if(exists $slgo{$i}{$tmp[0]})
				{
#				my $n = (split /\//, $tmp[2])[0];
					my ($n) = $tmp[2] =~ /^(\d+)\D*/;
					$level_file{$i}->print("$slgo{$i}{$tmp[0]}\t$n\t$tmp[6]\n");
				}
			}
		}
	}elsif($opts{ud} eq "diff"){
		while(<PFC>)
		{
			chomp;
			my @tmp = split /\t/;
			for my $i (2..4){
				if(exists $slgo{$i}{$tmp[0]})
				{
					my @genes = split /;/, $tmp[6];
					my (@up_gene, @down_gene);
					my ($up, $down) = (0, 0);
					foreach my $i(@genes)
					{
						if($gene_ud{$i} > 0)
						{
							$up ++;
							push(@up_gene, $i);
						}else{
							$down ++;
							push(@down_gene, $i);
						}
					}
					my $up_gene = join ";", @up_gene;
					my $down_gene = join ";", @down_gene;
					$level_file{$i}->print("$slgo{$i}{$tmp[0]}\t$up\t$down\t$up_gene\t$down_gene\n");
				}
			}
		}
	}else{
		die $usage;
	}
	close PFC;
	
	push @xls, "$opts{op}.$t.xls";
}

for my $i (2..4){
	close $level_file{$i};
	`perl $drawclass -i $opts{op}.Level$i.xls -t $opts{ud} -p $opts{op}.Level$i`;
	`rm $opts{op}.plot.r -rf`;
}

foreach(@xls)
{
	my $i = $_;
	$i =~ s/\.xls$//;
	my ($type) = $i =~ /\.(\w+?)$/;
	my $pi = basename($i);
	open HTML, "> $i.html" or die $!;
	print HTML <<CT;
	<html>
	<head>
	<meta http-equiv="content-type" content="text/html; charset=utf-8">
	<title>GO Enrichment Analysis</title>
	<style type="text/css">
		a {
			text-decoration: none;
			outline: none;
		}
		a:hover {
			color:#FF0000;
			text-decoration:none;
			outline:none;
		}
		body {
			font-family: "Microsoft YaHei","微软雅黑","雅黑宋体","新宋体","宋体","Microsoft JhengHei","华文细黑",STHeiti,MingLiu;
			background-color: #FFFFFF;
			padding-left: 8%;
			padding-right: 8%;
		}
		table {
			font-size: 14px;
			width: 100%;
			max-width: 100%;
			background-color: transparent;
			border-collapse: collapse;
			border-spacing: 0;
			display: table;
			border: 1px solid #dddddd;
		}
		th {
			background-color: rgba(2,79,101,1);
			color: #fff;
		}
		td, th {
			text-align: center;
			padding: 5px;
			border: 1px solid #dddddd;
		}
		tr:hover{
			background-color: #f5f5f5;
		}
		table caption{
			font-weight: bold;
			color: rgba(2,79,101,1);
			font-size: 1.5em;
			padding-bottom: 12px;
		}
		.go_graph {
			text-align: center;
			font-weight: bold;
			font-size: 1.5em;
			color: rgba(2,79,101,1);
			padding-bottom: 12px;
		}
		#bt {
			font-size: 14px;
			position: fixed;
			right: 2%;
			bottom: 5%;
		}
	</style>
	<script type="text/javascript">
	<!--
	function reSize2() {
		try {
			parent.document.getElementsByTagName("iframe")[0].style.height = document.body.scrollHeight + 10;
			parent.parent.document.getElementsByTagName("iframe")[0].style.height = parent.document.body.scrollHeight;
		} catch(e) {}
	}

	preRow = null;
	preColor = null;
	function colorRow(trObj) {
		if (preRow != null) {
			preRow.style.backgroundColor = preColor;
		}
		preRow = trObj;
		preColor = trObj.style.backgroundColor;
		trObj.style.backgroundColor = "f2dede";
	}

	function diffColor(tables) {
		color = ["#FFFFFF", "#CCFF99"];
		for (i = 0; i < tables.length; i++) {
			trObj = tables[i].getElementsByTagName("tr");
			for (j = 1; j < trObj.length; j++) {
				trObj[j].style.backgroundColor = color[j % color.length];
			}
		}
	}

	function markColor(table) {
		trs = table.getElementsByTagName("tr");
			for (i = 1; i < trs.length; i++) {
				if(table.rows[i].cells[6].innerHTML < 0.05){
					//trs[i].style.fontWeight = "500";
					table.rows[i].cells[6].style.color = "#FF0000";
					table.rows[i].cells[6].style.fontWeight = "900";
				}
				if(table.rows[i].cells[5].innerHTML < 0.05){
					table.rows[i].cells[5].style.color = "#FF0000";
				}
			}
	}

	function showPer(tableObj) {
		trObj = tableObj.getElementsByTagName("tr");
		if (trObj.length < 2) {
			return;
		}
		sum1 = trObj[0].cells[3].innerHTML.replace(/^.*\\(([\\d]+)\\).*\$/, "\$1");
		sum2 = trObj[0].cells[4].innerHTML.replace(/^.*\\(([\\d]+)\\).*\$/, "\$1");
		/*
		if (trObj[0].cells.length > 4) {
		}
		if (trObj[0].cells.length > 4) {
			trObj[0].cells[2].innerHTML = "DEGs genes with pathway annotation (" + sum1 + ")";
			trObj[0].cells[3].innerHTML = "All genes with pathway annotation (" + sum2 + ")";
		}else{
			trObj[0].cells[2].innerHTML = "All genes with pathway annotation (" + sum1 + ")";
		}
		*/
		for (i = 1; i < trObj.length; i++) {
			trObj[i].cells[3].innerHTML += " (" + (Math.round(trObj[i].cells[3].innerHTML * 10000/ sum1) / 100) + "%)";
			trObj[i].cells[4].innerHTML += " (" + (Math.round(trObj[i].cells[4].innerHTML * 10000/ sum2) / 100) + "%)";
		}
	}


	window.onload = function() {
		setTimeout("reSize2()", 1);
	}
	//-->
	</script>
	
	</head>
		<body>
CT
	
	open FA, $_ or die $!;
	my $head = <FA>; chomp $head;
	my @head = split /\t/, $head;
	my $t1 = "<table><caption>$title GO Enrichment ($namespace{$type})</caption><tr><th>#</th><th>$head[0]</th><th>$head[1]</th><th>$head[2]</th></tr>";  #<th>$head[3]</th><th>$head[4]</th><th>$head[5]</th></tr>";
	my $t2 = "<table><caption>$title GO Enrichment ($namespace{$type}) Gene Details</caption><tr><th>#</th><th>$head[0]</th><th>$head[6]</th></tr>";
	my $index = 0;
	while(<FA>)
	{
		chomp;
		my @tmp = split /\t/;
		$index ++;
		my $genes = pop @tmp;
		$genes =~ s/\;/ /g;
		$tmp[-1] = sprintf("%.6f", $tmp[-1]);
		$tmp[-2] = sprintf("%.6f", $tmp[-2]);
		$t1 .= "<tr><td>$index</td><td><a href='#gene$index' title='click to view genes' onclick='javascript: colorRow(document.getElementsByTagName(\"table\")[1].rows[$index]);'>$tmp[0]</a></td><td style=\"text-align: left;\">$tmp[1]</td><td>$tmp[2]</td></tr>";# . (join "</td><td>", @tmp[2..$#tmp]) . "</td></tr>";
		$t2 .= "<tr><td>$index</td><td><a href='http://amigo.geneontology.org/amigo/medial_search?q=$tmp[0]' title='click to view GO' target='_blank'>$tmp[0]</a></td><td style=\"text-align: left;\"><a name='gene$index'></a>$genes</td></tr>";
	}
	close FA;
	$t1 .= "</table>";
	$t2 .= "</table><script type='text/javascript'>showPer(document.getElementsByTagName('table')[0]);markColor(document.getElementsByTagName('table')[0]);</script>";
	print HTML $t1."\n<br><hr>\n".$t2."\n";
	print HTML <<CT;
		<a id="bt" href="#">Back Top</a>
		</body>
	</html>
CT
}
