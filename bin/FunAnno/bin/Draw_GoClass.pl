#!/usr/bin/perl -w
use strict;

=head1 Procedure:
	draw_GoClass.pl

=head1 Parameters:
	-i Input file.
	-p Output file path and prefix.
	-t Type of picture you want to draw : diff / nodiff .

=head1 Example:
	perl draw_GoClass.pl -i xx -p xx -t diff
	perl draw_GoClass.pl -i xx -p xx -t nodiff

=cut

use Getopt::Long;
use List::Util qw(max);
my ($input,$output,$type);
GetOptions(
	"i=s"=>\$input,
	"p=s"=>\$output,
	"t=s"=>\$type,
);
die `pod2text $0` unless ($input&&$type&&$output);

if (( $type ne "diff" )&&( $type ne "nodiff" )){
	die "Type of picture you want to draw must be diff or nodiff.\n";
}

use SVG;
use File::Basename qw /basename/;
my $row=0;
my $max=0;
my %hash;
my %max_term = (
	biological_process => 30,
	cellular_component => 25,
	molecular_function => 20
);
open IN,"<$input";
<IN>;
while ( my $line = <IN> ){
	chomp $line;
	my @array = split /\t/,$line;
	next if(keys %{$hash{$array[0]}} >= $max_term{$array[0]});
	if ( $type eq "diff" ){
		$hash{$array[0]}{$array[1]} = "$array[2]_$array[3]";
		$max = $array[2] if($array[2] > $max);
		$max = $array[3] if($array[3] > $max);
	}else{
		$hash{$array[0]}{$array[1]} = "$array[2]";
		$max = $array[2] if($array[2] > $max);
	}
	$row++;
}
close IN;

my $file = basename($input);
my ($title, $level) = $file =~ /(\S+)\.(Level\d)\.xls/;
$title = "$level GO terms of $title";

my ($width,$height);
my $svg = SVG->new(width=>2100, height=>1000);
if ( $type eq "diff" ){
	$width =1600/($row*3+1);
}else{
	$width =1600/($row*2+1);
}
my $font_size=$width*2/3;

my $divisor;
my $max_1;
if($max <= 5){
	$divisor = $max;
}
elsif($max > 5 && $max <= 20){
	$divisor = 3;
}
elsif($max > 20 && $max <= 100){
	$divisor = 4;
}
else{
	$divisor = 5;
}
my $remainder = $max % $divisor;
if(0 == $remainder){
	$max_1= $max;
}
else{
	$max_1= $max + ( $divisor - $remainder);
}
$height = 500/$max_1;
#print "$max\t$max_1\t$height\n";

$svg->line(x1 => 300, y1 => 95, x2 => 300, y2 => 605, stroke=>'black',"stroke-width"=>2);
$svg->line(x1 => 300, y1 => 605, x2 => 1900, y2 => 605, stroke=>'black',"stroke-width"=>2);
$svg->line(x1 => 1900, y1 => 95, x2 => 1900, y2 => 605, stroke=>'black',"stroke-width"=>2);
$svg->line(x1 => 300, y1 => 95, x2 => 1900, y2 => 95, stroke=>'black',"stroke-width"=>2);
$svg->line(x1 => 125, y1 => 900, x2 => 1725, y2 => 900, stroke=>'black',"stroke-width"=>2);

my $yi=600;
for my $ii ( 0 .. $divisor ){
	my $each_height = 500 / $divisor;
	$yi=600-$ii*$each_height;
	my $out_i = $ii*$max_1/$divisor;
	$svg->line(x1 => 295, y1 => $yi, x2 => 300, y2 => $yi, stroke=>'black',"stroke-width"=>2);
	$svg->text(x => 280, y => $yi, width => $width, height => 50, "font-family"=>"Arial", "text-anchor"=>"end","font-size"=> "18", "-cdata" => "$out_i");
}

$svg->text(x => 210, y => 350, width => 30, height => 50, "font-family"=>"Arial", "text-anchor"=>"middle","font-size"=> "25", "-cdata" => "Num of Genes","transform"=>"rotate(-90, 210, 350)");
my $title_x = (2100 - length($title)) / 2;
$svg->text(x => $title_x, y => 50, width => 50, height => 30, "font-family"=>"Arial", "text-anchor"=>"middle","font-size"=> "25", "-cdata" => "$title");

if ( $type eq "diff" ){
	$svg->rect(x => 1920, y => 100, width => 15, height => 15, fill => "red");
	$svg->text(x => 1940, y => 112, width => 15, height => 15, "font-family"=>"Arial", "text-anchor"=>"start","font-size"=> "18", "-cdata" => "up");
	$svg->rect(x => 1920, y => 130, width => 15, height => 15, fill => "green");
	$svg->text(x => 1940, y => 142, width => 15, height => 15, "font-family"=>"Arial", "text-anchor"=>"start","font-size"=> "18", "-cdata" => "down");
}

my $x = 300;
my $locus=125;
my $xx;
my @col=("green","red","blue");
for my $k ( sort keys %hash ){
	$xx = $x-175;
	$svg->line(x1 => 300, y1 => 605, x2 => 125, y2 => 900, stroke=>'black',"stroke-width"=>2);
	my $mark_1=$x;
	my $Col      = shift @col;
	if ( $type eq "nodiff" ){
		for my $i ( sort {$hash{$k}{$b} <=> $hash{$k}{$a}} keys %{$hash{$k}} ){
			$x+=$width;
			my ($H,$h,$x1);
			$h = $hash{$k}{$i}*$height;
			$x1 = $x+$width/2;
			$H = 600-$h;
			$svg->rect(x => $x, y => $H, width => $width, height => $h, fill => "$Col");
			$svg->line(x1 => $x1, y1 => 605, x2 => $x1, y2 => 610, stroke=>'black',"stroke-width"=>2);
			$i =~ s/(.{35}).*/$1.../;
			$svg->text(x => $x1, y => 620, width => $width, height => 50, "font-family"=>"Arial", "text-anchor"=>"end","font-size"=> "15", "-cdata" => "$i","transform"=>"rotate(-60, $x1, 620)");
			$x+=$width;
		}
	}else{
		for my $i ( sort {max(split /\_/,$hash{$k}{$b}) <=> max(split /\_/,$hash{$k}{$a})} keys %{$hash{$k}} ){
			$x+=$width;
			my @array = split /\_/,$hash{$k}{$i};
			my ($x1,$H1,$h1,$H2,$h2);
			$H1 = 600-$array[0]*$height;
			$h1 = $array[0]*$height;
			$H2 = 600-$array[1]*$height;
			$h2 = $array[1]*$height;
			$x1 = $x+$width/2;
			$svg->rect(x => $x, y => $H1, width => $width, height => $h1, fill => "red");
			$x+=$width;
			$x1 = $x+$width/2;
			$svg->line(x1 => $x, y1 => 605, x2 => $x, y2 => 610, stroke=>'black',"stroke-width"=>2);
			$i =~ s/(.{35}).*/$1.../;
			$svg->text(x => $x, y => 610, width => $width, height => 50, "font-family"=>"Arial", "text-anchor"=>"end","font-size"=> "14", "-cdata" => "$i","transform"=>"rotate(-60, $x1, 620)");
			$svg->rect(x => $x, y => $H2, width => $width, height => $h2, fill => "green");
			$x+=$width;
		}
	}
	my $x1 = $x+$width;
	my $x2 = $x1-175;
	$locus = $locus+($x1-$mark_1)/2;
	$svg->text(x => $locus, y => 950, width => 30, height => 50, "font-family"=>"Arial", "text-anchor"=>"middle","font-size"=>22,"-cdata" => "$k");
	if ( $x1 eq "1900" ){
		$svg->line(x1 => 1900, y1 => 605, x2 => $x2, y2 => 900, stroke=>'black',"stroke-width"=>2);
	}else{
		my $x3=$x1-$width/2;
		$svg->line(x1 => $x3, y1 => 605, x2 => $x2, y2 => 900, stroke=>'black',"stroke-width"=>2);
	}
	$locus = $locus+($x1-$mark_1)/2;
}

my $out = $svg->xmlify;
open OUT,">$output\.svg";
print OUT $out;
close OUT;

`convert -density 300 $output\.svg $output\.png`;
