#!/usr/bin/env perl

###########################################################
## Author: junting xie
## Email : junting.xie@sagene.com.cn
## Date  : 2018-02-11
###########################################################

use utf8;
use strict;
use warnings;

use lib "/Bio/User/xiejunting/MetaG/ReportFactory3";
use FindBin qw/$Bin $Script/;
use Cwd qw/abs_path/;
use ReportFactory;
use ReportFactory::Basic qw/file2array dir2array/;

my $output = abs_path($ARGV[0]);
my $report = ReportFactory->new();
my $body   = $report->init(
    out_path => $output,
    out_name => "Report",
    type     => "宏基因组报告",
);

$report->append_presets("$Bin/MetaG_Report", "MetaG");

my $presets    = $report->{presets}->{MetaG}[0][0];
my $tree_seed  = $report->{presets}->{MetaG}[0][1];

###########################################################

## Set Values
my @sample_list = ("hsm2-1","msm2-1","ksm2-1","ssm2-1","hsm2-2","ksm2-2","msm2-2","ssm2-2","hsm2-3","ksm2-3","msm2-3","ssm2-3","hsm3-1","ksm3-1","msm3-1","ssm3-1","hsm3-2","ksm3-2","msm3-2","ssm3-2","hsm3-3","ksm3-3","msm3-3","ssm3-3");
my @tax_list = ("Kingdom","Phylum","Class","Order","Family","Genus","Species");
my @beta_list = ("BrayCurtis","Euclidean","Mahalanobis","Jaccard","CorSpearman");
my @group_list  = ("hsm","msm","ksm","ssm");
my @sdiff_list  = ();
my @gdiff_list  = ();
my %groups;
@{$groups{hsm}} = ("hsm2-1", "hsm2-2", "hsm2-3", "hsm3-1", "hsm3-2", "hsm3-3");
@{$groups{msm}} = ("msm2-1", "msm2-2", "msm2-3", "msm3-1", "msm3-2", "msm3-3");
@{$groups{ksm}} = ("ksm2-1", "ksm2-2", "ksm2-3", "ksm3-1", "ksm3-2", "ksm3-3");
@{$groups{ssm}} = ("ssm2-1", "ssm2-2", "ssm2-3", "ssm3-1", "ssm3-2", "ssm3-3");
my $tmp;
my @tmp;
my %tmp;
my @tmp_image;
my @tmp_text;
my $tmp_images;
my $tmp_collapse;


###########################################################

$body->title(text => "项目概述", level => 1);

my @basic;
@{$basic[0]} = ("项目编号",  "SGM1052");
@{$basic[1]} = ("项目名称",  "宏基因组");
@{$basic[2]} = ("客户名称",  "杨小军");
@{$basic[3]} = ("公司/机构", "西北农林科技大学动物科技学院");
@{$basic[4]} = ("样品数目",  24);
$body->table(array => \@basic, caption => "基础信息", table_header => "none");

my @advan;
@{$advan[0]} = ("测序平台",   "Illumina");
@{$advan[1]} = ("测序仪型号", "NovaSeq");
$body->table(array => \@advan, caption => "高级信息", table_header => "none");

my @sample_info;
$tmp = join " <span class=\"project_info_separator\">&brvbar;</span> ", @sample_list;
@{$sample_info[0]} = ("样品名称", $tmp);
@tmp = @sdiff_list;
$tmp = join " <span class=\"project_info_separator\">&brvbar;</span> ", map { s#:# <span class="project_info_diff">&harr;</span> #g; $_; } @tmp;
@{$sample_info[1]} = ("样品间差异方案(对照/处理)", $tmp);
$tmp = join " <span class=\"project_info_separator\">&brvbar;</span> ", @group_list;
foreach my $group (@group_list){
my $g = join " <span class=\"project_info_and\">&amp;</span> ", @{$groups{$group}};
$tmp =~ s#$group# <span class="project_info_group_lab">$group : </span>$g#;
}
@{$sample_info[2]} = ("分组方案", $tmp);
@tmp = @gdiff_list;
$tmp = join " <span class=\"project_info_separator\">&brvbar;</span> ", map { s#:# <span class="project_info_diff">&harr;</span> #g; $_; } @tmp;
@{$sample_info[3]} = ("分组间差异方案(对照/处理)", $tmp);
$body->table(array => \@sample_info, caption => "样品信息", table_header => "none");

###########################################################
@tmp=();
$body->title(text => "宏基因组概念", level => 1);
$body->text(text => $presets->{metagenome}->{intro});
$body->text(text => $presets->{metagenome}->{diff});

$body->title(text => "建库测序流程", level => 1);
$body->text(text => $presets->{experiment_workflow}->{intro});
$body->image(image => $presets->{experiment_workflow}->{image}, text => $presets->{experiment_workflow}->{image_text});

$body->title(text => "分析流程", level => 1);
$body->text(text => $presets->{analysis_workflow}->{intro});
$body->image(image => $presets->{analysis_workflow}->{image}, text => $presets->{analysis_workflow}->{image_text});

$body->simList(text => $presets->{analysis_workflow}->{text}, list => $presets->{analysis_workflow}->{list});

############################################################
$body->title(text => "测序数据统计", level => 1);
$body->simList(text => $presets->{clean_data}->{intro}, list => $presets->{clean_data}->{list});

$body->title(text => "原始数据FastQC", level => 2);
$body->text(text => $presets->{fastqc}->{result});
%{$tmp[0]} = (title => "fastQC", desc => "结果文件夹", href => $presets->{fastqc}->{list});
$body->advList(array => \@tmp);

$body->title(text => "数据过滤(统计)", level => 2);
my $table = file2array(file => $presets->{Summary}->{dat_table});
my %hf = $body->table(array => $table, caption => "过滤前后Reads数据统计", file_link => $presets->{Summary}->{dat_table}, footer => 1);
$hf{footer}->simList(list => $presets->{Summary}->{footer});

$table = file2array(file => $presets->{Filter}->{flt_table});
%hf = $body->table(array => $table, caption => "过滤前后数据统计", file_link => $presets->{Filter}->{flt_table}, footer => 1);
$hf{footer}->simList(list => $presets->{Filter}->{footer});

$table = file2array(file => $presets->{FltCompose}->{cp_table});
%hf = $body->table(array => $table, caption => "Reads(过滤)组成统计", file_link => $presets->{FltCompose}->{cp_table}, footer => 1);
$hf{footer}->simList(list => $presets->{FltCompose}->{footer});
$body->image(image => $presets->{FltCompose}->{image});

$body->title(text => "样品数据过滤(统计)", level => 3);
$body->text(text => $presets->{smpFlt}->{text});
%{$tmp[0]} = (title => "样品过滤统计", desc => "结果文件夹", href => $presets->{smpFlt}->{href});
$body->advList(array => \@tmp);
@tmp = ("%NAME%/%NAME%.old.png", "%NAME%/%NAME%.new.png");
@tmp_text = ("%NAME% 过滤前", "%NAME% 过滤后");
$tmp_images = dir2array(samples => \@sample_list, path => $presets->{smpFlt}->{href}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, header => 1, footer => 1);
$hf{header}->text(text => $presets->{smpFlt}->{header});
$hf{footer}->simList(list => $presets->{smpFlt}->{footer});

##########################################################################################
$body->title(text => "宏基因组组装", level => 1);
$body->text(text => $presets->{metaG_asm}->{intro});
$tmp_collapse = $body->collapse(title => "Assembly原理", god => $report);
%tmp = $tmp_collapse->panel();
$tmp{body}->simList(text => $presets->{metaG_asm}->{asm_list_text},list => $presets->{metaG_asm}->{asm_list}, order =>0);
$tmp{body}->image(image => $presets->{metaG_asm}->{asm_image});
$tmp{body}->text(text => $presets->{metaG_asm}->{asm_text});

$body->title(text => "单样品组装", level => 2);
$body->simList(text => $presets->{smp_asm}->{intro}, list => $presets->{smp_asm}->{list});
@tmp = ();
%{$tmp[0]} = (title => "单样品组装", desc => "结果文件夹", href => $presets->{smp_asm}->{result}->{href});
$body->advList(array => \@tmp);
$table = file2array(file => $presets->{smp_asm}->{result}->{eval_table}, cols => "1-10");
%hf = $body->table(array => $table, caption => "组装结果统计", file_link => $presets->{smp_asm}->{result}->{eval_table}, footer => 1);
$hf{footer}->simList(list => $presets->{smp_asm}->{result}->{list});
@tmp = ();
@tmp = ("%NAME%_Length.PNG");
@tmp_text = ("%NAME% 长度分布图");
$tmp_images = dir2array(samples => \@sample_list, path => $presets->{smp_asm}->{result}->{href}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, header => 1, footer => 1);
$hf{header}->text(text => $presets->{smp_asm}->{result}->{len_header});
$hf{footer}->simList(list => $presets->{smp_asm}->{result}->{len_footer});
@tmp = ();
@tmp = ("%NAME%_GC.PNG");
@tmp_text = ("%NAME% GC Ratio分布图");
$tmp_images = dir2array(samples => \@sample_list, path => $presets->{smp_asm}->{result}->{href}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, header => 1, footer => 1);
$hf{header}->text(text => $presets->{smp_asm}->{result}->{gc_header});
$hf{footer}->simList(list => $presets->{smp_asm}->{result}->{gc_footer});

$body->title(text => "混合组装", level => 2);
$body->simList(text => $presets->{mix_asm}->{intro}, list => $presets->{mix_asm}->{list});
@tmp = ();
%{$tmp[0]} = (title => "混合组装", desc => "结果文件夹", href => $presets->{mix_asm}->{result}->{href});
$body->advList(array => \@tmp);
$table = file2array(file => $presets->{mix_asm}->{result}->{eval_table});
%hf = $body->table(array => $table, caption => "组装结果统计", file_link => $presets->{mix_asm}->{result}->{eval_table}, footer => 1);
$hf{footer}->simList(list => $presets->{mix_asm}->{result}->{list});
@tmp = ();
@tmp = ("%NAME%_Length.PNG");
@tmp_text = ("%NAME% 长度分布图");
my @class = ("All","Samples","Mix");
$tmp_images = dir2array(samples => \@class, path => $presets->{mix_asm}->{result}->{href}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, header => 1, footer => 1);
$hf{header}->text(text => $presets->{mix_asm}->{result}->{len_header});
$hf{footer}->simList(list => $presets->{mix_asm}->{result}->{len_footer});
@tmp = ();
@tmp = ("%NAME%_GC.PNG");
@tmp_text = ("%NAME% GC Ratio分布图");
$tmp_images = dir2array(samples => \@class, path => $presets->{mix_asm}->{result}->{href}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, header => 1, footer => 1);
$hf{header}->text(text => $presets->{mix_asm}->{result}->{gc_header});
$hf{footer}->simList(list => $presets->{mix_asm}->{result}->{gc_footer});

#$body->title(text => "组装比对评价", level => 2);
########################################################################################
$body->title(text => "物种分析", level => 1);
$body->text(text => $presets->{classify}->{intro});

$body->title(text => "物种注释", level => 2);
$body->text(text => $presets->{assign}->{intro});
@tmp = ();
%{$tmp[0]} = (title => "物种注释", desc => "结果文件夹", href => $presets->{assign}->{result}->{href});
$body->advList(array => \@tmp);
$table = file2array(file => $presets->{assign}->{result}->{assign_table});
%hf = $body->table(array => $table, caption => "注释结果表格", file_link => $presets->{assign}->{result}->{assign_table}, footer => 1, table_header => "none" );
$hf{footer}->simList(list => $presets->{assign}->{result}->{assign_footer});

$body->text(text => $presets->{krona}->{result}->{krona_text});
@tmp = ();
%{$tmp[0]} = (title => "Krona", desc => "HTML", href => $presets->{krona}->{result}->{krona_html});
#%{$tmp[1]} = (title => "Krona", desc => "PDF", href => $presets->{krona}->{result}->{krona_pdf});
$body->advList(array => \@tmp);
$body->image(image => $presets->{krona}->{result}->{krona_image});
$body->simList(list => $presets->{krona}->{result}->{krona_footer_list});

#################################################################################
$body->title(text => "物种定量", level => 2);
$body->text(text => $presets->{quantify}->{text});
@tmp = ();
%{$tmp[0]} = (title => "物种定量", desc => "结果文件夹", href => $presets->{quantify}->{result}->{href});
$body->advList(array => \@tmp);
#$table = file2array(file => $presets->{smp_asm}->{result}->{eval_table}, cols => "1-10");
#%hf = $body->table(array => $table, caption => "组装结果统计", file_link => $presets->{smp_asm}->{result}->{eval_table}, footer => 1);
#$hf{footer}->simList(list => $presets->{smp_asm}->{result}->{list});

#$body->title(text => "丰度柱状图", level => 2);
#$body->text(text => $presets->{topbar}->{text});
#@tmp = ();
#@tmp = ("Phylum.%NAME%_Top10Bar.PNG");
#@tmp_text = ("%NAME% Top10 物种丰度柱状图");
#$tmp_images = dir2array(samples => \@sample_list, path => $presets->{topbar}->{result}->{emphref}, pattern => \@tmp, text => \@tmp_text);
#%hf = $body->images(array => $tmp_images, god => $report, header => 1, footer => 1);
#$hf{header}->text(text => $presets->{topbar}->{result}->{len_header});
#$hf{footer}->simList(list => $presets->{topbar}->{result}->{len_footer});

################################################
$body->title(text => "丰度韦恩图(Venn)", level => 2);
$body->text(text => $presets->{venn}->{intro});
@tmp = ();
@tmp = ("%NAME%_ksm-V.S-hsm-V.S-msm-V.S-ssm_Venn.png");
@tmp_text = ("%NAME% 物种丰度韦恩图");
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{venn}->{result}->{emphref}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, header => 1, footer => 1);
$hf{header}->text(text => $presets->{venn}->{result}->{len_header});
$hf{footer}->simList(list => $presets->{venn}->{result}->{len_footer});

################################################
$body->title(text => "物种相对丰度", level => 2);
$body->text(text => $presets->{rltabd}->{text});
@tmp = ();
%{$tmp[0]} = (title => "物种定量", desc => "结果文件夹", href => $presets->{rltabd}->{result}->{href});
$body->advList(array => \@tmp);

$body->title(text => "相对丰度热图", level => 2);
$body->text(text => $presets->{heatmap}->{text});
@tmp = ();
@tmp = ("%NAME%_Dominant_Heatmap.png",  "%NAME%_Dominant_cHeatmap.png", "%NAME%_Dominant_rHeatmap.png");
@tmp_text = ("%NAME% 物种相对丰度热图", "%NAME% 物种相对丰度热图(Column)", "%NAME% 物种相对丰度热图(Row)");
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{heatmap}->{result}->{emphref}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, header => 1, footer => 1);
$hf{header}->text(text => $presets->{heatmap}->{result}->{len_header});
$hf{footer}->simList(list => $presets->{heatmap}->{result}->{len_footer});
$body->text(text => $presets->{heatmap}->{result}->{explain});

$body->title(text => "相对丰度堆叠图", level => 2);
$body->text(text => $presets->{stack}->{text});
@tmp = ();
@tmp = ("%NAME%_Stack.png");
@tmp_text = ("%NAME% 物种相对丰度堆叠图");
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{stack}->{result}->{emphref}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, header => 1, footer => 1);
$hf{header}->text(text => $presets->{stack}->{result}->{len_header});
$hf{footer}->simList(list => $presets->{stack}->{result}->{len_footer});
$body->text(text => $presets->{stack}->{result}->{explain});

$body->title(text => "相对丰度箱形图", level => 2);
$body->text(text => $presets->{box}->{text});
@tmp = ();
@tmp = ("%NAME%_Box.png");
@tmp_text = ("%NAME% 物种相对丰度箱形图");
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{box}->{result}->{emphref}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, header => 1, footer => 1);
$hf{header}->text(text => $presets->{box}->{result}->{len_header});
$hf{footer}->simList(list => $presets->{box}->{result}->{len_footer});
$body->text(text => $presets->{box}->{result}->{explain});

################################################
#$body->title(text => "Alpha 多样性", level => 1);
#$body->text(text => $presets->{alpha}->{intro});

################################################
$body->title(text => "(物种)样品间差异", level => 1);
$body->text(text => $presets->{beta}->{intro});
$body->text(text => $presets->{beta}->{text});
@tmp = ();
%{$tmp[0]} = (title => "物种差异分析", desc => "结果文件夹", href => $presets->{beta}->{result}->{href});
$body->advList(array => \@tmp);
$table = ();
$table = file2array(file => $presets->{beta}->{result}->{div_table}, cols => "1-10");
$body->table(array => $table, caption => "(物种)样品差异距离矩阵(Bray Curtis)", file_link => $presets->{beta}->{result}->{div_table});

my %pillow = $body->pillow(tabs => \@beta_list, god => $report);
# BrayCurtis
@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant_BrayCurtis_div.png");
@tmp_text = ();
@tmp_text = ("%NAME% (物种)样品间差异距离矩阵热图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{beta}->{result}->{exphref}, pattern => \@tmp, text => \@tmp_text);
$pillow{"BrayCurtis"}->images(array => $tmp_images, god => $report);
@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant_BrayCurtis_pcoa.png");
@tmp_text = ();
@tmp_text = ("%NAME% (物种)样品间PCoA图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{beta}->{result}->{exphref}, pattern => \@tmp, text => \@tmp_text);
$pillow{"BrayCurtis"}->images(array => $tmp_images, god => $report);
@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant_BrayCurtis_ClusterDendrogram.png");
@tmp_text = ();
@tmp_text = ("%NAME% (物种)样品间聚类树图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{beta}->{result}->{exphref}, pattern => \@tmp, text => \@tmp_text);
$pillow{"BrayCurtis"}->images(array => $tmp_images, god => $report);

# Euclidean
@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant_Euclidean_div.png");
@tmp_text = ();
@tmp_text = ("%NAME% (物种)样品间差异距离矩阵热图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{beta}->{result}->{exphref}, pattern => \@tmp, text => \@tmp_text);
$pillow{"Euclidean"}->images(array => $tmp_images, god => $report);
@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant_Euclidean_pca.png");
@tmp_text = ();
@tmp_text = ("%NAME% (物种)样品间PCoA图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{beta}->{result}->{exphref}, pattern => \@tmp, text => \@tmp_text);
$pillow{"Euclidean"}->images(array => $tmp_images, god => $report);
@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant_Euclidean_ClusterDendrogram.png");
@tmp_text = ();
@tmp_text = ("%NAME% (物种)样品间聚类树图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{beta}->{result}->{exphref}, pattern => \@tmp, text => \@tmp_text);
$pillow{"Euclidean"}->images(array => $tmp_images, god => $report);

# Mahalanobis
@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant_Mahalanobis_div.png");
@tmp_text = ();
@tmp_text = ("%NAME% (物种)样品间差异距离矩阵热图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{beta}->{result}->{exphref}, pattern => \@tmp, text => \@tmp_text);
$pillow{"Mahalanobis"}->images(array => $tmp_images, god => $report);
@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant_Mahalanobis_pcoa.png");
@tmp_text = ();
@tmp_text = ("%NAME% (物种)样品间PCoA图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{beta}->{result}->{exphref}, pattern => \@tmp, text => \@tmp_text);
$pillow{"Mahalanobis"}->images(array => $tmp_images, god => $report);
@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant_Mahalanobis_ClusterDendrogram.png");
@tmp_text = ();
@tmp_text = ("%NAME% (物种)样品间聚类树图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{beta}->{result}->{exphref}, pattern => \@tmp, text => \@tmp_text);
$pillow{"Mahalanobis"}->images(array => $tmp_images, god => $report);

# Jaccard
@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant_Jaccard_div.png");
@tmp_text = ();
@tmp_text = ("%NAME% (物种)样品间差异距离矩阵热图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{beta}->{result}->{exphref}, pattern => \@tmp, text => \@tmp_text);
$pillow{"Jaccard"}->images(array => $tmp_images, god => $report);
@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant_Jaccard_pcoa.png");
@tmp_text = ();
@tmp_text = ("%NAME% (物种)样品间差异距离矩阵热图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{beta}->{result}->{exphref}, pattern => \@tmp, text => \@tmp_text);
$pillow{"Jaccard"}->images(array => $tmp_images, god => $report);
@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant_Jaccard_ClusterDendrogram.png");
@tmp_text = ();
@tmp_text = ("%NAME% (物种)样品间聚类树图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{beta}->{result}->{exphref}, pattern => \@tmp, text => \@tmp_text);
$pillow{"Jaccard"}->images(array => $tmp_images, god => $report);

# CorSpearman
@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant_CorSpearman_div.png");
@tmp_text = ();
@tmp_text = ("%NAME% (物种)样品间差异距离矩阵热图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{beta}->{result}->{exphref}, pattern => \@tmp, text => \@tmp_text);
$pillow{"CorSpearman"}->images(array => $tmp_images, god => $report);
################################################
$body->title(text => "物种相关性分析", level => 1);
$body->text(text => $presets->{cor}->{intro});
$body->text(text => $presets->{cor}->{text});
@tmp = ();
%{$tmp[0]} = (title => "物种相关性分析", desc => "结果文件夹", href => $presets->{cor}->{result}->{href});
$body->advList(array => \@tmp);
$table = ();
$table = file2array(file => $presets->{cor}->{result}->{emp_table}, cols => "1-8");
$body->table(array => $table, caption => "物种相关矩阵(Spearman)", file_link => $presets->{cor}->{result}->{emp_table});

@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant_CorSpearman.png");
@tmp_text = ();
@tmp_text = ("%NAME% 物种相关矩阵热图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{cor}->{result}->{emphref}, pattern => \@tmp, text => \@tmp_text);
$body->images(array => $tmp_images, god => $report);

###########################################################
$body->title(text => "物种差异分析", level => 1);
$body->simList(text => $presets->{taxdiff}->{intro}, list => $presets->{taxdiff}->{list});
@tmp = ();
%{$tmp[0]} = (title => "物种差异分析", desc => "结果文件夹", href => $presets->{taxdiff}->{result}->{href});
$body->advList(array => \@tmp);

############################################
$body->title(text => "多分组差异检验(ANOVA)", level => 2);
$body->text(text => $presets->{lcgd}->{intro});
$body->text(text => $presets->{lcgd}->{result}->{anova_text});
$table = ();
$table = file2array(file => $presets->{lcgd}->{result}->{anova_table}, cols => "1-3,10-12");
%hf = $body->table(array => $table, caption => "ANOVA统计表格", file_link => $presets->{lcgd}->{result}->{anova_table}, footer => 1);
$hf{footer}->simList(list => $presets->{lcgd}->{result}->{anova_footer});
@tmp = ();
@tmp = ("AvgDepth_%NAME%_Dominant.ANOVA.p.filtered.0.05.png");
@tmp_text = ();
@tmp_text = ("%NAME% ANOVA多分组差异比较箱形图");
$tmp_images = ();
$tmp_images = dir2array(samples => \@tax_list, path => $presets->{lcgd}->{result}->{exphref}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, header => 1, footer => 1);
$hf{header}->text(text => $presets->{lcgd}->{result}->{anova_header});
$hf{footer}->simList(list => $presets->{lcgd}->{result}->{anova_explain});

############################################
$body->title(text => "两两差异检验(T-Test)", level => 2);
$body->text(text => $presets->{lcgd}->{intro});
$body->text(text => $presets->{lcgd}->{result}->{ttest_text});
$table = ();
$table = file2array(file => $presets->{lcgd}->{result}->{anova_table}, cols => "1-3,10-12");
%hf = $body->table(array => $table, caption => "T-Test统计表格", file_link => $presets->{lcgd}->{result}->{ttest_table}, footer => 1);
$hf{footer}->simList(list => $presets->{lcgd}->{result}->{ttest_footer});
$body->image(image => $presets->{lcgd}->{result}->{ttest_expplot}, text => $presets->{lcgd}->{result}->{ttest_header});

############################################
$body->title(text => "LEfSe 分析与可视化", level => 2);
$body->text(text => $presets->{lefse}->{intro});
$body->text(text => $presets->{lefse}->{text});

@tmp = ();
%{$tmp[0]} = (title => "LEfSe分析", desc => "结果文件夹", href => $presets->{lefse}->{result}->{href});
$body->advList(array => \@tmp);
$table = ();
$table = file2array(file => $presets->{lefse}->{result}->{res_table});
%hf = $body->table(array => $table, caption => "预测结果统计", file_link => $presets->{lefse}->{result}->{href}, footer => 1);
$hf{footer}->simList(list => $presets->{lefse}->{result}->{res_table_footer});

$body->image(image => $presets->{lefse}->{result}->{res_image}, text => $presets->{lefse}->{result}->{res_image_text});
$body->simList(list => $presets->{lefse}->{result}->{res_image_footer});

############################################
#$body->title(text => "Metastats 分析与可视化", level => 2);
#$body->text(text => $presets->{lcgd}->{intro});
#$body->text(text => $presets->{lcgd}->{result}->{ttest_text});
###########################################################
$body->title(text => "基因预测与分析", level => 1);
$body->title(text => "基因序列统计", level => 2);
$body->simList(text => $presets->{genepred}->{intro}, list => $presets->{genepred}->{list});
@tmp = ();
%{$tmp[0]} = (title => "基因预测", desc => "结果文件夹", href => $presets->{genepred}->{result}->{href});
$body->advList(array => \@tmp);
$table = ();
$table = file2array(file => $presets->{genepred}->{eval_table});
%hf = $body->table(array => $table, caption => "预测结果统计", file_link => $presets->{genepred}->{eval_table}, footer => 1);
$hf{footer}->simList(list => $presets->{genepred}->{eval_list});
@tmp = ();
@tmp = ("%NAME%FNA_Length.PNG","%NAME%FAA_Length.PNG");
@tmp_text = ("%NAME% 核酸长度分布图","%NAME% 氨基酸长度分布");
my @gclass = ("Raw","FLT");
$tmp_images = dir2array(samples => \@gclass, path => $presets->{genepred}->{result}->{href}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, header => 1, footer => 1);
$hf{header}->text(text => $presets->{genepred}->{len_header});
$hf{footer}->simList(list => $presets->{genepred}->{len_footer});

$body->image(image => $presets->{genepred}->{type_image}, text => $presets->{genepred}->{type_text});

#############################################################
$body->title(text => "基因丰度以及差异分析", level => 1);
$body->text(text => $presets->{geneab}->{intro});
@tmp = ();
%{$tmp[0]} = (title => "基因丰度", desc => "结果文件夹", href => $presets->{geneab}->{result}->{href});
$body->advList(array => \@tmp);
$table = ();
$table = file2array(file => $presets->{geneab}->{result}->{table}, cols => "1-10");
$body->table(array => $table, caption => "基因丰度结果表格", file_link => $presets->{geneab}->{result}->{table});
$body->image(image => $presets->{geneab}->{result}->{venn_image}, text => $presets->{geneab}->{result}->{venn_text});

##############################################################
$body->title(text => "基因丰度热图", level => 2);
$body->text(text => $presets->{geneheatmap}->{intro});
@tmp = ();
%{$tmp[0]} = (title => "基因丰度热图", desc => "结果文件夹", href => $presets->{geneheatmap}->{result}->{href});
$body->advList(array => \@tmp);
@tmp = ();
@tmp = ("%NAME%_Dominant_Heatmap.png",  "%NAME%_Dominant_cHeatmap.png", "%NAME%_Dominant_rHeatmap.png");
@tmp_text = ("%NAME% 基因相对丰度热图", "%NAME% 基因相对丰度热图(Column)", "%NAME% 基因相对丰度热图(Row)");
my @geneclass = ("Gene");
$tmp_images = dir2array(samples => \@geneclass, path => $presets->{geneheatmap}->{result}->{href}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, footer => 1);
$hf{footer}->simList(list => $presets->{geneheatmap}->{result}->{footer});
#$body->text(text => $presets->{heatmap}->{result}->{explain});

##########################################################
$body->title(text => "(基因)样品间差异分析", level => 2);
$body->text(text => $presets->{genediv}->{intro});
@tmp = ();
%{$tmp[0]} = (title => "样品间差异", desc => "结果文件夹", href => $presets->{genediv}->{result}->{href});
$body->advList(array => \@tmp);
$table = ();
$table = file2array(file => $presets->{genediv}->{result}->{table}, cols=>"1-10");
$body->table(array => $table, caption => "样品间距离矩阵", file_link => $presets->{genediv}->{result}->{table});
@tmp = ();
@tmp = ("%NAME%_Dominant_Euclidean_ClusterDendrogram.png");
@tmp_text = ("%NAME% 聚类树图");
$tmp_images = dir2array(samples => \@geneclass, path => $presets->{genediv}->{result}->{href}, pattern => \@tmp, text => \@tmp_text);
$body->images(array => $tmp_images, god => $report);

##########################################################
$body->title(text => "(基因)样品间相关性分析", level => 2);
$body->text(text => $presets->{genecor}->{intro});
@tmp = ();
%{$tmp[0]} = (title => "样品间相关性", desc => "结果文件夹", href => $presets->{genecor}->{result}->{href});
$body->advList(array => \@tmp);
$table = ();
$table = file2array(file => $presets->{genecor}->{result}->{div_table}, cols=>"1-10");
$body->table(array => $table, caption => "样品间距离矩阵", file_link => $presets->{genecor}->{result}->{div_table});
@tmp = ();
@tmp = ("%NAME%_Dominant_CorSpearman_div.png","%NAME%_Dominant_CorSpearman_div.Half.png");
@tmp_text = ("%NAME% 矩阵热图","%NAME% 三角热图");
$tmp_images = dir2array(samples => \@geneclass, path => $presets->{genecor}->{result}->{href}, pattern => \@tmp, text => \@tmp_text);
$body->images(array => $tmp_images, god => $report);

##########################################################
$body->title(text => "基因间相关性分析", level => 2);
#$body->text(text => $presets->{genecor}->{intro});
@tmp = ();
%{$tmp[0]} = (title => "基因间相关性", desc => "结果文件夹", href => $presets->{genecor}->{result}->{href});
$body->advList(array => \@tmp);
$table = ();
$table = file2array(file => $presets->{genecor}->{result}->{gene_table}, cols=>"1-10");
$body->table(array => $table, caption => "样品间距离矩阵", file_link => $presets->{genecor}->{result}->{gene_table});
@tmp = ();
@tmp = ("%NAME%_Dominant_CorSpearman.png");
@tmp_text = ("%NAME% 矩阵热图");
$tmp_images = dir2array(samples => \@geneclass, path => $presets->{genecor}->{result}->{href}, pattern => \@tmp, text => \@tmp_text);
$body->images(array => $tmp_images, god => $report);

########################################################################
$body->title(text => "(预测)基因功能注释", level => 1);
$body->simList(text => $presets->{func}->{intro}, list => $presets->{func}->{list});
$body->text(text => $presets->{func}->{text});
@tmp = ();
%{$tmp[0]} = (title => "功能注释", desc => "结果文件夹", href => $presets->{func}->{result}->{href});
$body->advList(array => \@tmp);

$body->title(text => "功能注释基础统计", level => 2);
$body->text(text => $presets->{funstats}->{text});
@tmp = ();
%{$tmp[0]} = (title => "功能注释", desc => "结果文件夹", href => $presets->{funstats}->{result}->{href});
$body->advList(array => \@tmp);
$table = ();
$table = file2array(file => $presets->{funstats}->{result}->{table});
%hf = $body->table(array => $table, caption => "注释统计结果表格", file_link => $presets->{funstats}->{result}->{table}, footer => 1);
$hf{footer}->simList(list => $presets->{funstats}->{result}->{footer});
$body->image(image => $presets->{funstats}->{result}->{venn}, text => $presets->{funstats}->{result}->{header});

###########################################################
$body->title(text => "KEGG", level => 2);
$body->text(text => $presets->{kegg}->{intro});
$table = ();
$table = file2array(file => $presets->{kegg}->{infos});
$body->table(array => $table, caption => "");
$body->text(text => $presets->{kegg}->{text});
@tmp = ();
%{$tmp[0]} = (title => "KEGG注释", desc => "结果文件夹", href => $presets->{kegg}->{result}->{href});
%{$tmp[1]} = (title => "KEGG map", desc => "结果文件夹", href => $presets->{kegg}->{result}->{maphref});
$body->advList(array => \@tmp);

$body->text(text => $presets->{kegg}->{anno_text});
$table = ();
$table = file2array(file => $presets->{kegg}->{result}->{anno_table}, cols=>"1-8");
%hf = $body->table(array => $table, caption => "注释结果表格", file_link => $presets->{kegg}->{result}->{anno_table},footer => 1);
$hf{footer}->simList(list => $presets->{kegg}->{result}->{anno_footer});

$body->text(text => $presets->{kegg}->{path_text});
$table = ();
$table = file2array(file => $presets->{kegg}->{result}->{path_table}, cols=>"1-5");
%hf = $body->table(array => $table, caption => "KEGG 注释通路表格", file_link => $presets->{kegg}->{result}->{path_table},footer => 1);
$hf{footer}->simList(list => $presets->{kegg}->{result}->{path_footer});
@tmp = ();
%{$tmp[0]} = (title => "KEGG Pahtway 注释统计", desc => "结果HTML", href => $presets->{kegg}->{result}->{anno_html});
$body->advList(array => \@tmp);
my @kegg_lst=("LevelA","LevelB","Pathway");
@tmp = ();
@tmp = ("SG.KEGG_%NAME%_FreqBar.png");
@tmp_text = ("Level %NAME% 注释(频率)统计柱状图");
$tmp_images = dir2array(samples => \@kegg_lst, path => $presets->{kegg}->{result}->{href}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, footer => 1);
$hf{footer}->simList(list => $presets->{kegg}->{result}->{anno_image_text});

$body->simList(text => $presets->{kegg}->{map_text}, list => $presets->{kegg}->{map_list});
$body->image(image => $presets->{kegg}->{result}->{example}, text => $presets->{kegg}->{result}->{example_text});
@tmp = ();
%{$tmp[0]} = (title => "map", desc => "结果文件夹", href => $presets->{kegg}->{result}->{maphref});
$body->advList(array => \@tmp);

#$body->text(text => $presets->{krona}->{result}->{krona_text});
$body->simList(text => $presets->{krona}->{result}->{krona_text},list => $presets->{kegg}->{result}->{krona_list});
@tmp = ();
%{$tmp[0]} = (title => "Krona", desc => "HTML", href => $presets->{kegg}->{result}->{krona_html});
#%{$tmp[1]} = (title => "Krona", desc => "PDF", href => $presets->{kegg}->{result}->{krona_pdf});
$body->advList(array => \@tmp);
$body->image(image => $presets->{kegg}->{result}->{krona_png});
$body->simList(list => $presets->{kegg}->{result}->{krona_footer});
###########################################################
$body->title(text => "EggNOG", level => 2);
$body->text(text => $presets->{eggnog}->{intro});
$table = ();
$table = file2array(file => $presets->{eggnog}->{infos});
$body->table(array => $table);
$body->text(text => $presets->{eggnog}->{text});
@tmp = ();
%{$tmp[0]} = (title => "EggNOG注释", desc => "结果文件夹", href => $presets->{eggnog}->{result}->{href});
$body->advList(array => \@tmp);

$body->text(text => $presets->{eggnog}->{anno_text});
$table = ();
$table = file2array(file => $presets->{eggnog}->{result}->{eggnog_table},cols => "1,5,7-9");
%hf = $body->table(array => $table, caption => "EggNOG注释结果表格", file_link => $presets->{eggnog}->{result}->{eggnog_table}, footer => 1);
$hf{footer}->simList(list => $presets->{eggnog}->{result}->{eggnog_footer});
$table = ();
$table = file2array(file => $presets->{eggnog}->{result}->{cog_table},cols => "1,2,4");
%hf = $body->table(array => $table, caption => "COG注释结果表格", file_link => $presets->{eggnog}->{result}->{cog_table},footer => 1);
$hf{footer}->simList(list => $presets->{eggnog}->{result}->{cog_footer});

my @eggnog_lst=("COG","EggNOG");
@tmp = ();
@tmp = ("SG.COG_%NAME%_FreqBar.png");
@tmp_text = ("Level %NAME% 注释(频率)统计柱状图");
$tmp_images = dir2array(samples => \@eggnog_lst, path => $presets->{eggnog}->{result}->{href}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, footer => 1);
$hf{footer}->simList(list => $presets->{eggnog}->{result}->{cog_image_footer});

######################################################################################
$body->title(text => "CAZy", level => 2);
$body->text(text => $presets->{cazy}->{intro});
$table = ();
$table = file2array(file => $presets->{cazy}->{infos});
$body->table(array => $table);
$body->text(text => $presets->{cazy}->{text});
@tmp = ();
%{$tmp[0]} = (title => "CAZy注释", desc => "结果文件夹", href => $presets->{cazy}->{result}->{href});
$body->advList(array => \@tmp);

$body->text(text => $presets->{cazy}->{anno_text});
$table = ();
$table = file2array(file => $presets->{cazy}->{result}->{cazy_table}, cols => "1,2,5,6,7");
%hf = $body->table(array => $table, caption => "CAZy注释结果表格", file_link => $presets->{cazy}->{result}->{cazy_table}, footer => 1);
$hf{footer}->simList(list => $presets->{cazy}->{result}->{cazy_footer});
$table = ();
$table = file2array(file => $presets->{cazy}->{result}->{class_table});
%hf = $body->table(array => $table, caption => "CAZy-Class注释结果表格", file_link => $presets->{cazy}->{result}->{class_table},footer => 1);
$hf{footer}->simList(list => $presets->{cazy}->{result}->{class_footer});
$table = ();
$table = file2array(file => $presets->{cazy}->{result}->{family_table});
%hf = $body->table(array => $table, caption => "CAZy-Family注释结果表格", file_link => $presets->{cazy}->{result}->{family_table},footer => 1);
$hf{footer}->simList(list => $presets->{cazy}->{result}->{family_footer});

my @cazy_lst=("Level1","Level2");
@tmp = ();
@tmp = ("SG.CAZy_%NAME%_FreqBar.png");
@tmp_text = ("%NAME% 注释(频率)统计柱状图");
$tmp_images = dir2array(samples => \@cazy_lst, path => $presets->{cazy}->{result}->{href}, pattern => \@tmp, text => \@tmp_text);
%hf = $body->images(array => $tmp_images, god => $report, footer => 1);
$hf{footer}->simList(list => $presets->{cazy}->{result}->{cazy_image_footer});

######################################################################################
# Function Differential Analysis


######################################################################################
# Function TO Species

######################################################################################
# Generate Report
$report->write_report(PDF => 1);

