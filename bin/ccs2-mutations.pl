#!/usr/bin/perl -w

################################################################################
# Using PacBio SMRT sequencing to track Pol V replication activity
# Copyright (C) 2018 New England Biolabs, Inc.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

use strict;
use Getopt::Long;

my $o_realign = 1;

GetOptions(
    "realign!" => \$o_realign,
    );

### command-line arguments
if( @ARGV == 0 )
{
    print "usage: $0 [options] input.bam reference.fasta\n\n";
    print "options\n";
    print "  --realign\t(default '$o_realign')\n";
    exit;
}

my $bamfile = shift @ARGV;
my $reffile = shift @ARGV;

### load reference
my $refs = load_references($reffile);

### compute homopolymer regions
my $hpmap = build_hpmap($refs);

# data structures for parsing CIGAR string
my %digit = (
    "0" => 1,
    "1" => 1,
    "2" => 1,
    "3" => 1,
    "4" => 1,
    "5" => 1,
    "6" => 1,
    "7" => 1,
    "8" => 1,
    "9" => 1 );

my %operator = (
    "M" => 1,
    "I" => 1,
    "D" => 1,
    "=" => 1,
    "X" => 1 );

print "Movie,ZMW,Ref,Pos,Length,Type,QV,AtEnd,RefBP,AltBP,IndelType,InHP,Bases,HPLen,HPChar\n";

open(BAM,"samtools view $bamfile |") || die "Can't open '$bamfile'";

while( my $line = <BAM> )
{
    my $s1 = ""; # aligned reference
    my $s2 = ""; # aligned read
    my $bq = ""; # base quality

    ### parse fields
    my ( $qname,
	 $flag,
	 $rname,
	 $pos,
	 $mapq,
	 $cigar,
	 $rnext,
	 $pnext,
	 $tlen,
	 $read,
	 $qual,
	 @tag ) = split(/\t/,$line);

    ### fwd/rev strand
    my $strand = ( $flag & 0x10 ) ? 1 : 0;

    ### skip unmapped reads
    next if( $flag & 0x4 );

    ### skip secondary alignment
    next if( $flag & 0x100 );

    ### skip supplementary alignment
    next if( $flag & 0x800 );

    # ----- process CIGAR ------------------------------------------------------

    my $rseq = $$refs{$rname};

    my $p1  = $pos - 1; # start of aligned reference
    my $p2  = 0;        # start of aligned read
    my $num = 0;        # length of the block to operate
    my $len = 0;        # read length
    
    my @cigar = split(//,$cigar);
    
    while( @cigar )
    {
	my $c = shift @cigar;
	
	if( exists $digit{$c} )
	{
	    $num .= $c;
	    next;
	}
	
	if( $c eq "M" || $c eq "=" || $c eq "X" )
	{
	    $s1 .= substr($rseq,$p1,$num);
	    $s2 .= substr($read,$p2,$num);
	    $bq .= substr($qual,$p2,$num);

	    $p1 += $num;
	    $p2 += $num;
	    
	    $len += $num;
	}
	elsif( $c eq "D" )
	{
	    $s1 .= substr($rseq,$p1,$num);
	    $s2 .= ("-" x $num);
	    $bq .= (" " x $num);
	    
	    $p1 += $num;
	    
	    $len += $num;
	}
	elsif( $c eq "I" )
	{
	    $s1 .= ("-" x $num);
	    $s2 .= substr($read,$p2,$num);
	    $bq .= substr($qual,$p2,$num);

	    $p2 += $num;
	}
	elsif( $c eq "S" )
	{
	    $p2 += $num;
	}
	else
	{
	    print STDERR "error: unknown operator in CIGAR string '$c'\n\n";
	    exit;
	}
	
	$num = "";
    }
    
    my $ref_start = $pos;
    my $ccs_start = 1;
    
    my ($movie,$zmw,$other) = split(/\//,$qname);
    
    mutations($rname,$movie,$zmw,$s1,$s2,$bq,$ref_start-1,$ccs_start-1,$strand);
}

################################################################################
#                             SUBROUTINES                                      #
################################################################################

sub mutations {
    my ($ref_name,$movie,$zmw,$rseq,$read,$qual,$ref_start,$ccs_start,$strand) = @_;
    
    my $len1 = length($rseq);
    my $len2 = length($read);
    
    if( $len1 != $len2 )
    {
	print "error :: sequence length does not match\n\n\n";
	exit;
    }
    
    ### convert strings to arrays
    my @rseq = split(//,uc($rseq));
    my @read = split(//,uc($read));
    
    ### convert quality values
    my @qual = map { ord($_) - 33 } split(//,$qual);

    ### update deletion qual values
    fill_in_deletion_qual(\@read,\@qual,$strand);

    ### init position in the reference sequence
    my $rp = $ref_start - 1;
    
    ### init position in the sequencing read
    my $sp = $ccs_start - 1;
    
    ### init indel position in the reference sequence
    my $pindel = $rp;
    
    ### init auxilary indel data
    my @bases = ();
    
    for( my $i = 0; $i < $len1; $i++ )
    {
	my $rb = $rseq[$i];
	my $sb = $read[$i];
	
	### update reference position
	$rp++ if( $rb ne "-" );
	
	### update sequencing read position
	$sp++ if( $sb ne "-" );
	
	### update indel position
	$pindel = $rp if( $rb ne "-" && $sb ne "-" );

	if( $rb ne "-" && $sb ne "-" )
	{
	    ### SNP
	    
	    ### print call info
	    printf( "%s,%s,%s,%i,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",
		    $movie,
		    $zmw,
		    $ref_name,
		    $rp + 1,
		    1,
		    "SNP",
		    qual_value(\@rseq,\@read,\@qual,$strand,$i,\@bases),
		    ( $i == 0 || $i + 1 == $len1 ? "True" : "False" ),
		    $rb,
		    $sb,
		    "NA",
		    "NA",
		    "NA",
		    "NA",
		    "NA" );

	    @bases = ();
	}
	else
	{
	    ### INDEL
	    if( $rb ne "-" && $sb eq "-" )
	    {
		### DELETION
		
		### indel bases
		push @bases, $rb;
		
		### [end of alignment] or [end of deletion]
		if( $i + 1 == $len1 || $read[$i+1] ne "-" )
		{
		    ### homopolymer info (not currently used in downstream analysis)
		    my $HPLen = $$hpmap{$ref_name}[$pindel+1];
		    my $InHP = ( $HPLen > 1 ) ? "True" : "False";
		    my $HPChar = substr($$refs{$ref_name},$pindel+1,1);
		    
		    ### print call info
		    printf( "%s,%s,%s,%i,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",
			    $movie,
			    $zmw,
			    $ref_name,
			    $pindel + 1,
			    scalar(@bases),
			    "INDEL",
			    qual_value(\@rseq,\@read,\@qual,$strand,$i,\@bases),
			    "False",
			    "NA",
			    "NA",
			    "Deletion",
			    $InHP,
			    join("",@bases),
			    $HPLen,
			    $HPChar );

		    @bases = ();
		}
	    }
	    elsif( $rb eq "-" && $sb ne "-" )
	    {
		### INSERTION
		
		### indel bases
		push @bases, $sb;
		
		### [end of alignment] or [end of insertion]
		if( $i + 1 == $len1 || $rseq[$i+1] ne "-" )
		{
		    ### homopolymer info (not currently used in downstream analysis)
		    my $HPLen = $$hpmap{$ref_name}[$pindel+1];
		    my $InHP = ( $HPLen > 1 ) ? "True" : "False";
		    my $HPChar = substr($$refs{$ref_name},$pindel+1,1);
		    
		    ### print call info
		    printf( "%s,%s,%s,%i,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",
			    $movie,
			    $zmw,
			    $ref_name,
			    $pindel + 1,
			    scalar(@bases),
			    "INDEL",
			    qual_value(\@rseq,\@read,\@qual,$strand,$i,\@bases),
			    "False",
			    "NA",
			    "NA",
			    "Insertion",
			    $InHP,
			    join("",@bases),
			    $HPLen,
			    $HPChar );

		    @bases = ();
		}
	    }
	    else
	    {
		print STDERR "error :: rb='$rb' sb='$sb'\n\n";
		exit;
	    }
	}
    }
}

sub load_references {
    my ($file) = @_;
    
    my %refs = ();
    my $name = "";
    
    open(FA,$file) || die "Can't open '$file': $!";
    
    while( my $line = <FA> )
    {
	chomp($line);
	
	if( substr($line,0,1) eq ">" )
	{
	    $name = substr($line,1);
	}
	else
	{
	    $refs{$name} .= $line;
	}
    }
    
    close(FA);
    
    return \%refs;
}

sub build_hpmap {
    my ($refs) = @_;
    
    my %hpmap = ();

    for my $name ( keys %$refs )
    {
	my $seq = $$refs{$name};
	
	for( my $i = 0; $i < length($seq); $i++ )
	{
	    my $len = 1;
	    my $k = $i + 1;
	    
	    while( $k < length($seq) && substr($seq,$i,1) eq substr($seq,$k,1) )
	    {
		$len++;
		$k++;
	    }

	    push @{$hpmap{$name}}, $len;
	}
    }

    return \%hpmap;
}

sub fill_in_deletion_qual {
    my ($read,$qual,$strand) = @_;
    
    if( $strand == 0 )
    {
	### forward strand
	### AAGTTAAACAAAATTATTTCTAGACCAAT >>> sequencing direction
	### AAGTTAAAC-------------GACCAAT
	###               <-------*
	
	for( my $i = @$qual-2; $i >= 0; $i-- )
	{
	    if( $$read[$i] eq "-" )
	    {
		$$qual[$i] = $$qual[$i+1];
	    }
	}
    }
    else
    {
	### reverse strand
	### AAGTTAAACAAAATTATTTCTAGACCAAT <<< sequencing direction
	### AAGTTAAAC-------------GACCAAT
	###         *------>

	for( my $i = 1; $i < @$qual; $i++ )
	{
	    if( $$read[$i] eq "-" )
	    {
		$$qual[$i] = $$qual[$i-1];
	    }
	}
    }
}

sub qual_value {
    my ($rseq,$read,$qual,$strand,$i,$bases) = @_;
    
    if( $$rseq[$i] ne "-" && $$read[$i] ne "-" )
    {
	### SNP
	# return $i;
	return $$qual[$i];
    }
    else
    {
	### indel length
	my $len = scalar(@$bases);

	if( $$rseq[$i] ne "-" && $$read[$i] eq "-" )
	{
	    ### Deletion

	    ### 1234567890123456
	    ###          i
	    ### ATGCTGATGTAGATAA ref
	    ### ATGCTG----AGATAA ccs

	    if( $strand == 0 )
	    {
		# return $i+1;
		return $$qual[$i+1];
	    }
	    else
	    {
		my $offset = ( $o_realign ? qual_position_offset($rseq,$read,$i,$len) : 0 );
		my $qpos = ( $offset == 0 ? $i-$len : $i+$offset );

		# return $qpos;
		return $$qual[$qpos];
	    }
	}
	elsif( $$rseq[$i] eq "-" && $$read[$i] ne "-" )
	{
	    ### Insertion
	    ### 1234567890123456
	    ###          i
	    ### ATGCTG----AGATAA ref
	    ### ATGCTGATGTAGATAA ccs

	    my $hp = 1;

	    for( my $i = 0; $i < @$bases-1; $i++ )
	    {
		if( $$bases[$i] ne $$bases[$i+1] )
		{
		    $hp = 0;
		    last;
		}
	    }

	    if( $hp )
	    {
		### homopolymer indel

		if( $strand == 0 )
		{
		    ### first base
		    # return $i-$len+1;
		    return $$qual[$i-$len+1];
		}
		else
		{
		    ### last base
		    my $offset = ( $o_realign ? qual_position_offset($read,$rseq,$i,$len) : 0 );
		    my $qpos = $i + $offset;

		    # return $qpos;
		    return $$qual[$qpos];
		}
	    }
	    else
	    {
		if( $strand == 0 )
		{
		    ### take the lowest qual value for non-homopolymer indel
		    my $min_qv = 93;
		    
		    for( my $k = 0; $k < @$bases; $k++ )
		    {
			if( $$qual[$i-$k] < $min_qv )
			{
			    $min_qv = $$qual[$i-$k];
			}
		    }
		    
		    # return $i;
		    return $min_qv;
		}
		else
		{
		    my $offset = ( $o_realign ? qual_position_offset($read,$rseq,$i,$len) : 0 );
		    my $qpos = $i + $offset;
		    
		    ### take the lowest qual value for non-homopolymer indel
		    my $min_qv = 93;
		    
		    for( my $k = 0; $k < @$bases; $k++ )
		    {
			if( $$qual[$qpos-$k] < $min_qv )
			{
			    $min_qv = $$qual[$qpos-$k];
			}
		    }
		    
		    # return $qpos;
		    return $min_qv;
		}
	    }
	}
    }
}


sub qual_position_offset {
    my ($seq1,$seq2,$i,$len) = @_;
    
    ### seq1 is expected to be a non-gapped sequence
    ### seq2 is expected to be a gapped sequence
    ### use seq1/seq2 accordingly for deletions/insertions
    
    ### 1234567890123456
    ###          i
    ### ATGCTGATGTAGATAA seq1
    ### ATGCTG----AGATAA seq2
    
    ### swap as many bases as possible without introducing mismatches
    ### (end-to-end realignment)
    my $k = $i;

    while( $$seq1[$k-$len+1] eq $$seq2[$k+1] )
    {
	$k++;
    }

    ### check indel before and after swap
    my $iold = join("",@{$seq1}[($i-$len+1)..$i]);
    my $inew = join("",@{$seq1}[($k-$len+1)..$k]);
    
    ### only proceed when indels before and after swap match each other
    ### effectively, we apply realign to homopolymers and tandem repeats only
    ### everything else is left untouched
    
    if( $inew eq $iold )
    {
	return ($k-$i);
    }
    
    return 0;
}
