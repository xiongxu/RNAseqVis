#!/usr/bin/perl

# Copyright (c)   xiongxu 2011-2012
# Writer:         xuxiong <xuxiong19880610@163.com>
# Program Date:   2011.12.25
# Modifier:       xuxiong <xuxiong19880610@163.com>
# Last Modified:  2012.11.27

use strict;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use GD::Simple;
use warnings;

#Before writing your programmeyou must write the detailed time /discriptions /parameter and it's explanation,Meanwhile,annotation your programme in English if possible.

if (scalar(@ARGV)<4) {
	print "\n  When you run the program ,the sam|bam file needs to be sorted,the bam file can use the samtools(sort) to sort,the commadline:\"samtools sort infile.bam infile_sort &\",if the file format is sam,you can use the commmandline:\"sort -k 3,3 -k 4,4n infile.sam > infile_sort.sam &\" \n\n";
	print "Example1:\n  perl $0 /data2/verticillium_dahliae_vdls/gene/verticillium_dahliae_vdls.17_1_transcripts.gff /data2/human/seq_20120222_filter_reads/verticillium_dahliae_vdls/RNA_0224/sequence_3_clean.fq_0224/accepted_hits.sam 0228 5000 0 1\n";
	print "Example2:\n  perl $0 /data2/human/hg18.gff /data1/test/xux/test_software/human_1026/SpliceMap/A3_L1_1.fq/good_hits.sam 0227 5000 5000 1 /data2/human/kgXref.txt\n";
	print "Example3:\n  perl $0 /data2/arabidopsis/TAIR10/GFF/TAIR10_GFF3_genes.gff /data1/test/xux/arabidopsis_seq_20120614_seq_20120618/arabidopsis_0711/Ninanjie_mRNA_W_small_clean.fq_0710/accepted_hits.sam.uniq  Ninanjie_mRNA_W_small 10 10 0 /data2/arabidopsis/TAIR10/function_annotation/gene_aliases.20101027\n";
	print "Example4:\n  perl $0 /home/xuxiong/workdir/HSA_for_cufflinks.gff /home/xuxiong/workdir/lvliyaELL_L1_I008_tophat_mis2/accepted_hits.bam outfile 1000 1000 1\n";
	print "\nUsage:\n  perl $0 gff sam|bam outfile promoter_length downstream_length merge_dentical_reads[0|1] kgXref ......\n\n";
	exit;
}

###############Time_start##########
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print STDERR "\nStart Time :[$Time_Start]\n\n";
###################################
my $current_dir = `pwd`;chomp($current_dir);

my $gff = $ARGV[0];
my $infile = $ARGV[1];
my $outfile = $ARGV[2];
my $promoter=$ARGV[3] ;
my $downstream=$ARGV[4] ;
my $merge_dentical_reads= $ARGV[5];
my $kgXrefinfile = $ARGV[6] if (defined $ARGV[6] ) ;

my %kgXref=();
&load_kgXref($kgXrefinfile,\%kgXref) if (defined $ARGV[6] ) ;
my %gene=();
&load_gff($gff,\%gene,\%kgXref);
#print Dumper %gene;
my %cluster_gene=();
&cluster_gene(\%gene,\%cluster_gene,$promoter,$downstream);
my %know_junction=();
&read_known_junction(\%gene,\%know_junction);
&output_clustered_gene(\%cluster_gene,$outfile);
##############test#####
#undef(%gene);
#my @array_test;
#for (my $i=1;$i<=50;$i++) {
#&draw(\@array_test,$cluster_gene{"Chr1"}->[$i],$cluster_gene{"Chr1"}->[$i]->[0]->[0]->[0]->[6] eq "+" ? $cluster_gene{"Chr1"}->[$i]->[0]->[0]->[0]->[3]-$promoter : $cluster_gene{"Chr1"}->[$i]->[0]->[0]->[0]->[3],max(map {$_->[0]->[0]->[4]} @{$cluster_gene{"Chr1"}->[$i]}),"Chr1_cluster$i.png")
#}
#exit;
#######################

my %distribution=();
my $direction_bool=0;
chdir($current_dir);
if ($infile=~/\.sam/) {
	&load_sorted_sam($infile,$outfile,\%cluster_gene,\%distribution,$direction_bool,$promoter,$downstream);
}
elsif($infile=~/\.bam$/){
	$infile="samtools view $infile|";
	&load_sorted_sam($infile,$outfile,\%cluster_gene,\%distribution,$direction_bool,$promoter,$downstream);
}

sub load_kgXref{
	my ($infile,$hash_ref)=@_;
	open (IN ,"<$infile") || die "Can't open $!";
	while (<IN>) {
		chomp;
		my @line=split(/\t+/,$_);
		my $geneID=shift(@line);
		if ($geneID=~/([\w\.]+)/) {
			$hash_ref->{$1}=join("\t",@line);
		}
		else{
			$hash_ref->{$geneID}=join("\t",@line);
		}
	}
	close IN ;
	print STDERR "Done load kgXref\n";
}

sub load_gff{
	my ($GFF_INFILE,$hash,$kgXref)=@_;
	my @gene_content=();
	my $RNA_index=0;
	my $gene_start=0;
	my $CHROMOSOME="";
	open (IN,"$GFF_INFILE") || die $!;
	while (<IN>) {
		chomp;
		next if (/^\#|^$/) ;
		my @line=split(/\s+/,$_);
		if ($line[2]=~/^gene|pseudogene|transposable_element_gene$/ ) {
			if (@gene_content){
				foreach my $mRNA (@gene_content[1..$#gene_content]) {
					my $mRNAID="";
					if ($mRNA->[0]->[-1]=~/^ID=([\w\.]+)\;/) {
						$mRNAID=$1 ;
						push @{$mRNA->[0]},$kgXref->{$mRNAID} if (exists($kgXref->{$mRNAID})) ;
					}
					@{$mRNA}[1..$#{$mRNA}]=sort {$a->[3]<=>$b->[3] or $a->[4]<=>$b->[4]} sort {$a->[2] cmp $b->[2]} @{$mRNA}[1..$#{$mRNA}];
				}
#				&build_index(\@gene_content,$promoter,10000,$hash,$CHROMOSOME);
				push @{$hash->{$CHROMOSOME}},[@gene_content];
			}
			$RNA_index=0;
			($CHROMOSOME,$gene_start)=@line[0,3];
			@gene_content=();
			push @{$gene_content[$RNA_index]},[@line];
		}
		elsif($line[2]=~/^mRNA|miRNA|mRNA_TE_gene|ncRNA|pseudogenic_transcript|rRNA|snoRNA|snRNA|tRNA|transcript$/){
			$RNA_index++;
			push @{$gene_content[$RNA_index]},[@line];
		}
		elsif($line[2]=~/^chromosome|region$/i){
			next;
		}
		else{
			push @{$gene_content[$RNA_index]},[@line];
		}
		if (eof(IN)) {
			if (@gene_content){
				foreach my $mRNA (@gene_content[1..$#gene_content]) {
					@{$mRNA}[1..$#{$mRNA}]=sort {$a->[3]<=>$b->[3]} sort {$a->[2] cmp $b->[2]}  @{$mRNA}[1..$#{$mRNA}];
				}
#				&build_index(\@gene_content,$promoter,10000,$hash,$CHROMOSOME);
				push @{$hash->{$CHROMOSOME}},[@gene_content];
			}
		}
	}
	close(IN);
	&sort_gene($hash);
	print STDERR "Done load gff\n";
}

sub sort_gene{
	my $hash_ref=shift;
	foreach my $TEMP_CHRO (keys %{$hash_ref}) {
		@{$hash_ref->{$TEMP_CHRO}}= sort {$a->[0]->[0]->[3]<=>$b->[0]->[0]->[3] or $a->[0]->[0]->[4]<=>$b->[0]->[0]->[4] } @{$hash_ref->{$TEMP_CHRO}};
	}
}

sub build_index {
	my ($GENE_CONTENT,$PROMOTER_LENGTH,$SIZE,$HASH,$CHROM)=@_;
	my ($blockNum_start,$blockNum_end);
	if($GENE_CONTENT->[0]->[0]->[6]eq"+"){
		$blockNum_start=int(($GENE_CONTENT->[0]->[0]->[3]-$PROMOTER_LENGTH)/$SIZE);
		$blockNum_end=int(($GENE_CONTENT->[0]->[0]->[4])/$SIZE);
	}
	else{
		$blockNum_start=int($GENE_CONTENT->[0]->[0]->[3]/$SIZE);
		$blockNum_end=int(($GENE_CONTENT->[0]->[0]->[4]+$PROMOTER_LENGTH)/$SIZE);
	}
	map {push(@{$HASH->{$CHROM}{$_}},[@{$GENE_CONTENT}])}($blockNum_start..$blockNum_end);
}

sub region_combination{
	my ($start1,$end1,$start2,$end2)=@_;
	if ($end2<$start1) {
		return;
	}else{
		return(min($start1,$start2),max($end1,$end2));
	}
}

sub cluster_gene{
	my ($gene_ref,$gene_cluster_ref,$PROMOTER_LENGTH,$DOWNSTREAM_LENGTH)=@_;
	foreach my $TEMP_CHRO (sort keys %{$gene_ref}) {
		my $cluster_index=0;
		my ($Start,$End,$current_Start,$current_End)=(0,0,0,0);
		foreach my $temp_gene (@{$gene_ref->{$TEMP_CHRO}}) {
			if ($temp_gene->[0]->[0]->[6] eq "+") {
				($Start,$End)=($temp_gene->[0]->[0]->[3]-$PROMOTER_LENGTH,$temp_gene->[0]->[0]->[4]+$DOWNSTREAM_LENGTH);
				$Start=1 if ($Start<0) ;
			}else{
				($Start,$End)=($temp_gene->[0]->[0]->[3]-$DOWNSTREAM_LENGTH,$temp_gene->[0]->[0]->[4]+$PROMOTER_LENGTH);
				$Start=1 if ($Start<0) ;
			}
			my @current_region=&region_combination($Start,$End,$current_Start,$current_End);
			if (@current_region) {
				($current_Start,$current_End)=@current_region;
			}else{
				($current_Start,$current_End)=($Start,$End);
				$cluster_index++;
			}
			push @{$gene_cluster_ref->{$TEMP_CHRO}->[$cluster_index-1]},[@{$temp_gene}];
		}
	}
	print STDERR "Done cluster gene\n";
}

sub read_known_junction {
	my ($gene_ref,$know_junc_ref)=@_;
	foreach my $TEMP_CHRO (sort keys %{$gene_ref}) {
		foreach my $temp_gene (@{$gene_ref->{$TEMP_CHRO}}) {
			foreach my $temp_RNA (@{$temp_gene}[1..$#{$temp_gene}]) {
				my $junction_start=0;
				foreach my $temp_exon (@{$temp_RNA}[1..$#{$temp_RNA}]) {
	#				print join("\t",@{$temp_exon}),"\n";
					if ($temp_RNA->[0]->[2]=~/^mRNA$/) {
						next if ($temp_exon->[2]=~/^exon|start_codon|stop_codon$/) ;
					}
					if ($junction_start!=0) {
						$know_junc_ref->{$TEMP_CHRO}{$junction_start}{$temp_exon->[3]-1-$junction_start}++;
					}
					$junction_start= $temp_exon->[4];
				}
			}
		}
	}
#	foreach my $TEMP_CHRO (sort keys(%{$know_junc_ref})) {
#		foreach my $junction_start (sort {$a<=>$b} keys(%{$know_junc_ref->{$TEMP_CHRO}})) {
#			foreach my $junction_length (sort {$a<=>$b} keys(%{$know_junc_ref->{$TEMP_CHRO}{$junction_start}})) {
#				print $TEMP_CHRO,"\t",$junction_start,"\t",$junction_length,"\n";
#			}
#		}
#	}
	print STDERR "Done read_known_junction\n";
}

sub output_clustered_gene {
	my ($cluster_ref,$OUTFILE)=@_;
	open(OUT,">".$OUTFILE."cluster.gff") || die $!;
	foreach my $TEMP_CHRO (sort keys %{$cluster_ref}) {
		for(my $cluster_index=0;$cluster_index<@{$cluster_ref->{$TEMP_CHRO}};$cluster_index++) {
			print OUT ">$TEMP_CHRO\_cluster$cluster_index\n";
			foreach my $gene (@{$cluster_ref->{$TEMP_CHRO}->[$cluster_index]}) {
				foreach my $RNA (@{$gene}) {
					foreach my $temp_line (@{$RNA}) {
						print OUT join("\t",@{$temp_line}),"\n";
					}
				}
			}
		}
	}
	close(OUT);
	print STDERR "Done output cluster\n";
}

sub load_sorted_sam	{
	my ($SAM_FILE,$OUTFILE,$CLUSTER_GENE,$distribution,$directional,$promoter_length,$downstream_length)=@_;
	open(IN,$SAM_FILE) || die $!;
	open(OUT2,">".$OUTFILE.".intergenic") || die $!;
	my %intergenic_peak=();
	my @gene_reads=();
	my @intergenic_reads=();
	my $current_chr=(sort keys(%{$CLUSTER_GENE}))[0];
	my %reads_distribution=();
	my $cluster_index=0;
	mkdir("$current_dir/$current_chr");
	chdir("$current_dir/$current_chr");
	my $cluster_start=min(map {$_->[0]->[0]->[3]} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]});
	my $cluster_end=max(map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]});
	my $cluster_end_strand= (map {$_->[0]->[0]->[6]} grep {$_->[0]->[0]->[4]==$cluster_end} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]})[0];
	$cluster_start= $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[6] eq "+" ? $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$promoter_length : $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$downstream_length;
	$cluster_end= $cluster_end_strand eq "-" ? $cluster_end+$promoter_length : $cluster_end+$downstream_length;
	my $tag_chromosome_end=0;
	while (<IN>) {
		chomp;
		next if (/^\@|^\#/) ;
		my @line=split(/\t+/,$_);
		next if ($line[3]==0 or $line[2]=~/Zv/i) ;
		my $read_strand = $line[1] & (0x0010) ? "-" : "+";
#		next if ($line[5]!~/^[\dMNISP=XHD]+$/) ;
		my $read_end;
		if ($line[5]=~/^(\d+)M$/) {
			$read_end=$line[3]+$1-1;
		}
		elsif($line[5]=~/^(\d+)M(\d+)N(\d+)M$/){
			next if (min($1,$3)<8) ;
			$read_end=$line[3]+$1+$2+$3-1;
		}
		else{
			next;
		}
		if ($line[2] ne $current_chr) {
			mkdir("$current_dir/$line[2]");
			chdir("$current_dir/$line[2]");
			print STDERR "Done SAM $current_chr\n";
			@gene_reads=();
			@intergenic_reads=();
			$cluster_index=0;
#			$cluster_end=(sort {$a->[0]->[0]->[4]<=>$b->[0]->[0]->[4]} map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$line[2]}->[$cluster_index]})[-1]->[0]->[0]->[4];
			$cluster_end=max(map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$line[2]}->[$cluster_index]});
			$cluster_end_strand= (map {$_->[0]->[0]->[6]} grep {$_->[0]->[0]->[4]==$cluster_end} @{$CLUSTER_GENE->{$line[2]}->[$cluster_index]})[0];
			$cluster_start= $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[6] eq "+" ? $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$promoter_length : $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$downstream_length;
			$cluster_end= $cluster_end_strand eq "-" ? $cluster_end+$promoter_length : $cluster_end+$downstream_length;
		}
		if ($line[3]>$cluster_end ) {
			if ($cluster_index+1>$#{$CLUSTER_GENE->{$current_chr}}) {
				if ($tag_chromosome_end==0) {
					&draw(\@gene_reads,$CLUSTER_GENE->{$current_chr}->[$cluster_index],$cluster_start,$cluster_end,"$current_chr\_cluster$cluster_index.png",$current_chr) if (scalar(@gene_reads)>10) ;
				}
				$distribution->{"intergenic"}++;
				$tag_chromosome_end=1;
				next;
			}else{
				$tag_chromosome_end=0;
			}
			&draw(\@gene_reads,$CLUSTER_GENE->{$current_chr}->[$cluster_index],$cluster_start,$cluster_end,"$current_chr\_cluster$cluster_index.png",$current_chr) if (scalar(@gene_reads)>10) ;
			@gene_reads=();
			$cluster_index++;
			
			$cluster_end=max(map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]});
			$cluster_end_strand= (map {$_->[0]->[0]->[6]} grep {$_->[0]->[0]->[4]==$cluster_end} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]})[0];
			$cluster_start= $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[6] eq "+" ? $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$promoter_length : $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$downstream_length;
			$cluster_end= $cluster_end_strand eq "-" ? $cluster_end+$promoter_length : $cluster_end+$downstream_length;
			while ($line[3]>$cluster_end) {
				$cluster_index++;
				if ($cluster_index>$#{$CLUSTER_GENE->{$current_chr}}){
					$cluster_index--;
					$tag_chromosome_end=1;
					last;
				}
				$cluster_end=max(map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]});
				$cluster_end_strand= (map {$_->[0]->[0]->[6]} grep {$_->[0]->[0]->[4]==$cluster_end} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]})[0];
				$cluster_start= $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[6] eq "+" ? $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$promoter_length : $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$downstream_length;
				$cluster_end= $cluster_end_strand eq "-" ? $cluster_end+$promoter_length : $cluster_end+$downstream_length;
			}
		}
		
		if ($line[3]+sum(split(/[M|N]/,$line[5]))-1<$cluster_start) {
			$distribution->{"intergenic"}++;
			push @intergenic_reads,[($line[3],$read_strand,$line[5])];
		}
		else{
			if (@intergenic_reads) {
				&call_peak($current_chr,\@intergenic_reads,\%intergenic_peak,\*OUT2,5);
				@intergenic_reads=();
			}
			push @gene_reads,[($line[3],$read_strand,$line[5],$line[9])];
		}
		$current_chr=$line[2];
	}
	close(OUT2);
	close(IN);
}

sub call_peak {
	my ($chromosome,$array_ref,$peak_info,$OUTFILE_HANDEL,$Peak_Threshold)=@_;
	my @peak_array;
	map {push(@{$_},$_->[0]+sum(split(/[M|N]/,$_->[2]))-1)} @{$array_ref};
#	foreach my $reads (@{$array_ref}) {
#		my @matched = split(/[M|N]/,$reads->[-1]);
#		$reads->[-1]=$reads->[0]+sum(@matched)-1;
#	}
	@{$array_ref}= sort {$a->[0]<=>$b->[0] or $a->[-1]<=>$b->[-1]} @{$array_ref};
	my ($peak_index,$current_Start,$current_End)=(0,0,0);
	while (@{$array_ref}) {
		my @array_content=@{shift(@{$array_ref})};
		my @current_peak=&region_combination($array_content[0],$array_content[-1],$current_Start,$current_End);
		if (@current_peak) {
			($current_Start,$current_End)=@current_peak;
		}else{
			($current_Start,$current_End)=($array_content[0],$array_content[-1]);
			$peak_index++;
		}
		push @{$peak_array[$peak_index-1]},[@array_content];
	}
	for (my $i=0;$i<@peak_array ;$i++) {
		if (scalar(@{$peak_array[$i]})<=$Peak_Threshold) {
			splice(@peak_array,$i,1);
			$i--;
		}else{
			push @{$peak_info->{$chromosome}},[(min(map {$_->[0]} @{$peak_array[$i]}),max(map {$_->[-1]} @{$peak_array[$i]}),scalar(@{$peak_array[$i]}))];
			print $OUTFILE_HANDEL "$chromosome\tPeak$i\t",join("\t",(min(map {$_->[0]} @{$peak_array[$i]}),max(map {$_->[-1]} @{$peak_array[$i]}),scalar(@{$peak_array[$i]}))),"\n";
			print $OUTFILE_HANDEL join("\n",map {join("\t",@{$_})} @{$peak_array[$i]}),"\n";
		}
	}
	undef(@peak_array);
}

sub draw {
	my ($reads_ref,$cluster,$Cluster_start,$Cluster_end,$BASE_NAME,$chromosome)=@_;
	my $left=300;
	my $top=300;
	my $scale=1;
	my $scale_height=5000;
	my $scale_width=8000;
	my $Width=$scale_width/$scale+2*$left;

	my $cluster_Length=abs($Cluster_end-$Cluster_start);
	my $ratio=$scale_width/$cluster_Length;

	my $RNA_line_height=0.02;
	my $read_height=0.002;
	my $wig_height=0.1;
	my (@junction,@no_junction);
	my %junction_pos=();
	&differentiate_junction($reads_ref,\@junction,\@no_junction,\%junction_pos);
	my $novel_junction_flag=0;
	if (%junction_pos) {
#		print Dumper %junction_pos;
		foreach my $junction_start (sort {$a<=>$b} keys(%junction_pos)) {
			if (!exists($know_junction{$chromosome}{$junction_start})) {
				$novel_junction_flag=1;
#				print STDERR $chromosome,"\t",$junction_start,"\n";
				last;
			}else{
				foreach my $junction_length (sort {$a<=>$b} keys(%{$junction_pos{$junction_start}})) {
					if (!exists($know_junction{$chromosome}{$junction_start}{$junction_length})) {
						$novel_junction_flag=1;
#						print STDERR $chromosome,"\t",$junction_start,"\t",$junction_length,"\n";
						last;
					}
				}
			}
		}
	}
	return 0 if ($novel_junction_flag==0);

	my @wig_ref=();
	&tidy_sam_to_wig($reads_ref,\@wig_ref,$merge_dentical_reads);
	my $max_wig=max(map {$_->[2]} @wig_ref)+0.15;
	
	my $RNA_count=sum(map {$#{$_}} @{$cluster});

	my $Height=($RNA_count*$RNA_line_height+0.06+($max_wig+0.15)*$wig_height)*$scale_height/$scale+2*$top;
	
	my $img = GD::Simple->new($Width,$Height);

	$img->font(gdSmallFont);
	$img->bgcolor('blue');
	$img->fgcolor('blue');
	$img->fontsize(40);
	my $RNA_index=0;
	my $temp_end=10e10;
	foreach my $temp_gene (@{$cluster}) {
#		foreach my $temp_RNA (@{$temp_gene}[1..$#{$temp_gene}]) {
		foreach my $temp_RNA (@{$temp_gene}) {
#			next if ($temp_RNA->[0]->[2]=~/^mRNA|miRNA|mRNA_TE_gene|ncRNA|pseudogenic_transcript|rRNA|snoRNA|snRNA|tRNA|transcript$/) ;
			next if ($temp_RNA->[0]->[2]=~/^gene|pseudogene|transposable_element_gene$/) ;
			$RNA_index-- if ($temp_RNA->[0]->[3] > $temp_end+100*$scale/$ratio) ;
			$img->penSize(5,5);
			######draw RNA information###########
			$img->fgcolor('blue');
			$img->moveTo($left+($temp_RNA->[0]->[3]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.009)*$scale_height/$scale);
			$img->angle(0);
			$img->font('Times:italic');
			$img->string(join(" ",@{$temp_RNA->[0]}[0,3,4,6,-1]));
			######draw mRNA|ncRNA length#########
			$img->fgcolor('blue');
			$img->moveTo($left+($temp_RNA->[0]->[3]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.015)*$scale_height/$scale);
			$img->lineTo($left+($temp_RNA->[0]->[4]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.015)*$scale_height/$scale);
#			$img->line($left+($temp_RNA->[0]->[3]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.015)*$scale_height/$scale,$left+($temp_RNA->[0]->[4]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.015)*$scale_height/$scale);
#			print join("\t",@{$temp_RNA->[0]}),"\n";
			######draw each CDS|three_prime_UTR|five_prime_UTR|ncExon########
			foreach my $temp_exon (@{$temp_RNA}[1..$#{$temp_RNA}]) {
#				print join("\t",@{$temp_exon}),"\n";
				$img->penSize(1,1);
				if ($temp_gene->[0]->[0]->[2]=~/^mRNA$/) {
					next if ($temp_exon->[2]=~/^exon|start_codon|stop_codon$/) ;
				}
				if ($temp_exon->[2]=~/^three_prime_UTR|five_prime_UTR$/) {
					$img->bgcolor('fuchsia');
					$img->fgcolor('fuchsia');
				}else{
					$img->bgcolor('aqua');
					$img->fgcolor('aqua');
				}
				my $exon_x1=($temp_exon->[3]-$Cluster_start)*$ratio/$scale;
				my $exon_x2=($temp_exon->[4]-$Cluster_start)*$ratio/$scale;
				$img->rectangle($left+$exon_x1,$top+($RNA_line_height*$RNA_index+0.01)*$scale_height/$scale,$left+$exon_x2,$top+($RNA_line_height*$RNA_index+0.02)*$scale_height/$scale);
			}
			#####mark the transcript direction#####
			if ($temp_RNA->[0]->[6] eq "+") {
				$img->fgcolor('red');
				$img->penSize(3,3);
				$img->moveTo($left+($temp_RNA->[0]->[4]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.01)*$scale_height/$scale);
				$img->angle(45);
				$img->lineTo($left+($temp_RNA->[0]->[4]-$Cluster_start)*$ratio/$scale+50,$top+($RNA_line_height*$RNA_index+0.015)*$scale_height/$scale);
				$img->moveTo($left+($temp_RNA->[0]->[4]-$Cluster_start)*$ratio/$scale+50,$top+($RNA_line_height*$RNA_index+0.015)*$scale_height/$scale);
				$img->angle(135);
				$img->lineTo($left+($temp_RNA->[0]->[4]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.02)*$scale_height/$scale);
			}else{
				$img->fgcolor('blue');
				$img->penSize(3,3);
				$img->moveTo($left+($temp_RNA->[0]->[3]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.01)*$scale_height/$scale);
				$img->angle(135);
				$img->lineTo($left+($temp_RNA->[0]->[3]-$Cluster_start)*$ratio/$scale-50,$top+($RNA_line_height*$RNA_index+0.015)*$scale_height/$scale);
				$img->moveTo($left+($temp_RNA->[0]->[3]-$Cluster_start)*$ratio/$scale-50,$top+($RNA_line_height*$RNA_index+0.015)*$scale_height/$scale);
				$img->angle(45);
				$img->lineTo($left+($temp_RNA->[0]->[3]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.02)*$scale_height/$scale);
			}
			$RNA_index++;
			$temp_end = $temp_RNA->[-1]->[4];
		}
	}
	#############draw_RNA_frame##################
	$img->bgcolor(undef);
	$img->fgcolor('black');
	$img->penSize(1,1);
	$img->rectangle($left,$top,$left+$scale_width/$scale,$top+$RNA_index*$RNA_line_height*$scale_height/$scale);
	#############draw x-axis#####################
	$img->penSize(5,5);
	$img->fgcolor('black');
	my $locate=$RNA_index*$RNA_line_height+0.02;
	for (my $i=$Cluster_start;$i<=$Cluster_end ;$i+=int($cluster_Length*0.2)) {
		$img->moveTo($left+($i-$Cluster_start)*$ratio/$scale,$top+$locate*$scale_height/$scale);
		$img->lineTo($left+($i-$Cluster_start)*$ratio/$scale,$top+($locate-0.01)*$scale_height/$scale);
		$img->moveTo($left+($i-$Cluster_start)*$ratio/$scale-10,$top+($locate+0.015)*$scale_height/$scale);
		$img->angle(0);
#		$img->fontsize(40);
		$img->font('Times:italic');
		$img->string($i);
	}
	$img->moveTo($left,$top+$locate*$scale_height/$scale);
	$img->lineTo($left+$cluster_Length*$ratio/$scale,$top+$locate*$scale_height/$scale);
	$locate+=0.029;
	############draw junction####################
	$locate+=0.055;
	&draw_junction($left,$top,$scale_height,$scale_width,$scale,\%junction_pos,$locate,$ratio,$Cluster_start,$img);
	############draw wig#########################
	&draw_wig($wig_height,$left,$top,$scale_height,$scale_width,$scale,\@wig_ref,$locate,$ratio,$Cluster_start,$img,$max_wig);
	$locate+=$max_wig*$wig_height;

	open (OUT,">$BASE_NAME") or die $!;
	binmode OUT;
	print OUT $img->png;
	close OUT;
#	print STDERR  "Done draw $BASE_NAME\n";
}

sub differentiate_junction {
	my ($array_ref,$return_ref1,$return_ref2,$return_ref3)=@_;
	for (my $i=0;$i<@{$array_ref} ;$i++) {
		if ($array_ref->[$i]->[2]=~/N/) {
			push @{$return_ref1},$array_ref->[$i];
			my @tmp_size=split(/[M|N]/,$array_ref->[$i]->[2]);
			$return_ref3->{$array_ref->[$i]->[0]-1+$tmp_size[0]}{$tmp_size[1]}++;
		}else{
			push @{$return_ref2},$array_ref->[$i];
		}
	}
}

sub draw_junction{
	my ($Left,$Top,$Scale_height,$Scale_width,$Scale,$hash_ref,$Locate,$cluster_length_ratio,$Cluster_start,$IMG)=@_;
	$IMG->angle(0);
	foreach my $start (sort {$a<=>$b} keys(%{$hash_ref})) {
		foreach my $length (sort {$a<=>$b} keys(%{$hash_ref->{$start}})) {
			$IMG->moveTo($Left+($start+($length+1)/2-$Cluster_start)*$cluster_length_ratio/$Scale,$Top+$Locate*$Scale_height/$Scale);
			$IMG->bgcolor('purple');
			$IMG->arc($length*$cluster_length_ratio/$Scale,0.05*$Scale_height/$Scale,180,360,gdNoFill);
#			$IMG->moveTo($Left+($start+($length+1)/2-$Cluster_start)*$cluster_length_ratio/$Scale,$Top+($Locate-0.025)*$Scale_height/$Scale);
#			$IMG->string($hash_ref->{$start}{$length});
		}
	}
}

sub draw_wig {
	my ($Wig_height,$Left,$Top,$Scale_height,$Scale_width,$Scale,$array_ref,$Locate,$cluster_length_ratio,$Cluster_start,$IMG,$max)=@_;
	$IMG->penSize(1,1);
	$IMG->angle(0);
	$IMG->bgcolor(undef);
	$IMG->fgcolor('black');
	$IMG->rectangle($Left,$Top+$Locate*$Scale_height/$Scale,$Left+$Scale_width/$Scale,$Top+($Locate+$max)*$Scale_height/$Scale);
	my $poly = new GD::Polygon;
	$IMG->moveTo($Left+($Cluster_start-$Cluster_start)*$cluster_length_ratio/$Scale,$Top+$Locate*$Scale_height/$Scale);
	$poly->addPt($Left+($array_ref->[0]->[0]-$Cluster_start)*$cluster_length_ratio/$Scale,$Top+$Locate*$Scale_height/$Scale);
	foreach my $hit (@{$array_ref}) {
		$poly->addPt($Left+($hit->[0]-$Cluster_start)*$cluster_length_ratio/$Scale,$Top+($hit->[2]*$Wig_height+$Locate)*$Scale_height/$Scale);
		$poly->addPt($Left+($hit->[1]-$Cluster_start)*$cluster_length_ratio/$Scale,$Top+($hit->[2]*$Wig_height+$Locate)*$Scale_height/$Scale);
	}
	$poly->addPt($Left+($array_ref->[-1]->[1]-$Cluster_start)*$cluster_length_ratio/$Scale,$Top+$Locate*$Scale_height/$Scale);
	$IMG->fgcolor(undef);
	$IMG->bgcolor('blue');
	$IMG->polygon($poly);
}

sub color_gradient_ten () {#datanow,data1,data2,color1,color2#######颜色随值变化的函数
	my @svg_x=@_;
	my @out_color;
	if ($svg_x[0] >=$svg_x[2]) {
		@out_color=@{$svg_x[4]};
	}
	elsif ($svg_x[0] <=$svg_x[1]) {
		@out_color=@{$svg_x[3]};
	}
	else {
		my $new_red=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($svg_x[4]->[0]-$svg_x[3]->[0])+$svg_x[3]->[0]);
		my $new_gre=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($svg_x[4]->[1]-$svg_x[3]->[1])+$svg_x[3]->[1]);
		my $new_blu=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($svg_x[4]->[2]-$svg_x[3]->[2])+$svg_x[3]->[2]);
		@out_color=($new_red,$new_gre,$new_blu);
	}
	return @out_color;
}

sub reads_split_0 {
	my ($start,$string,$reads_start,$matched_length)=@_;
	@{$matched_length}=split(/[M|N]/,$string);
	undef(@{$reads_start});
	push @{$reads_start},$start;
	push @{$reads_start},map {$start+sum(@{$matched_length}[0..$_])} (0..$#{$matched_length}-1) if (scalar(@{$matched_length}>1));
}

sub tidy_sam_to_wig{
	my ($hash_ref,$arr_ref,$Merge_identical_reads)=@_;
	if ($Merge_identical_reads) {
	#################combine the identical reads ##################
		my %hash_ref_tmp;
		foreach my $s (@{$hash_ref}) {
			$hash_ref_tmp{join("\t",@{$s})}++;
		}
		@{$hash_ref}=();
		foreach my $t (sort {(split/\s+/,$a)[0]<=>(split/\s+/,$b)[0]} keys %hash_ref_tmp) {
			push @{$hash_ref},[(split(/\s+/,$t),$hash_ref_tmp{$t})];
		}
	}
	else{
		@{$hash_ref}=sort {$a->[0]<=>$b->[0]} @{$hash_ref};
	}
	################################################################
	my %start=();
	my %end=();
	for (my $k=0;$k<@{$hash_ref} ;$k++) {
		my @tmp_size=split(/[M|N]/,$hash_ref->[$k]->[2]);
		if ($Merge_identical_reads) {
			$start{$hash_ref->[$k]->[0]}+=$hash_ref->[$k]->[-1];
			$end{$hash_ref->[$k]->[0]+$tmp_size[0]-1}+=$hash_ref->[$k]->[-1];
			if (scalar(@tmp_size)>1) {
				$start{$hash_ref->[$k]->[0]+$tmp_size[0]+$tmp_size[1]}+=$hash_ref->[$k]->[-1];
				$end{$hash_ref->[$k]->[0]+$tmp_size[0]+$tmp_size[1]+$tmp_size[2]-1}+=$hash_ref->[$k]->[-1];
			}
		}else{
			$start{$hash_ref->[$k]->[0]}++;
			$end{$hash_ref->[$k]->[0]+$tmp_size[0]-1}++;
			if (scalar(@tmp_size)>1) {
				$start{$hash_ref->[$k]->[0]+$tmp_size[0]+$tmp_size[1]}++;
				$end{$hash_ref->[$k]->[0]+$tmp_size[0]+$tmp_size[1]+$tmp_size[2]-1}++;
			}
		}
	}
	my $count=0;
	my $prevkey=0;
	my %all_keys=();
	map {$all_keys{$_}++} keys(%start);
	map {$all_keys{$_}++} keys(%end);
	foreach my $pos (sort {$a<=>$b} keys(%all_keys)) {
		push @{$arr_ref},[($prevkey, $pos, $count)] if ($prevkey!=0);
		$prevkey=$pos;
		$count+=$start{$pos} if (exists($start{$pos}));
		$count-=$end{$pos} if (exists($end{$pos}));
	}
	for (my $k=1;$k<@{$arr_ref};$k++) {
		if ($arr_ref->[$k]->[0]==$arr_ref->[$k-1]->[1] && $arr_ref->[$k]->[2] == $arr_ref->[$k-1]->[2]) {
			splice(@{$arr_ref},$k-1,2,[($arr_ref->[$k-1]->[0],$arr_ref->[$k]->[1],$arr_ref->[$k]->[2])]);
			$k--;
		}
		if ($arr_ref->[$k]->[2]>0) {
			$arr_ref->[$k]->[2]=log($arr_ref->[$k]->[2])/log(10);
			$arr_ref->[$k]->[2]=0.1 if ($arr_ref->[$k]->[2]==0) ;
		}
	}
	return;
}

sub reads_split{
	my ($start,$string,$reads_start,$reads_center,$matched_length,$reads_count)=@_;
	my @array=split(/[M|N]/,$string);
#	print "@array\n";
	my @matched_index=grep {$_%2==0} (0..$#array);
	@{$matched_length}=map {$array[$_]} @matched_index;
	my $read_length=sum(@{$matched_length});
	@{$reads_count}=map {$_/$read_length} @{$matched_length};
	undef(@{$reads_start});
	push @{$reads_start},$start;
	push @{$reads_start},map {$start+sum(@array[0..$_])} grep {$_%2==1} (0..$#array) if (scalar(@array>1));
	@{$reads_center}=map {$reads_start->[$_]-1+int(($matched_length->[$_]+1)/2)} (0..$#{$reads_start});
}

sub AbsolutePath{		#获取指定目录或文件的决定路径
	my ($type,$input) = @_;
	my $return;
	if ($type eq 'dir'){
		my $pwd = `pwd`;
		chomp $pwd;
		chdir($input);
		$return = `pwd`;
		chomp $return;
		chdir($pwd);
	}
	elsif($type eq 'file'){
		my $pwd = `pwd`;
		chomp $pwd;
		my $dir=dirname($input);
		my $file=basename($input);
		chdir($dir);
		$return = `pwd`;
		chomp $return;
		$return .="\/".$file;
		chdir($pwd);
	}
	return $return;
}

###############Time_end###########
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print STDERR "\nEnd Time :[$Time_End]\n\n";

###############Sub_format_datetime
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
