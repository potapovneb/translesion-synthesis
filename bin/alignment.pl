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

### command-line arguments
if( @ARGV == 0 )
{
    print "usage: $0 input.bam reference.fasta\n\n";
    exit;
}

my $bamfile = shift @ARGV;
my $reffile = shift @ARGV;

### load reference
my $refs = load_references($reffile);

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

print "Movie,ZMW,Reference,Flag,MAPQ,AlnStart,AlnEnd,AlnLength,AlnCoverage,Correct,Substitution,Deletion,Insertion,Clipped,CigarD,CigarI\n";

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

    my $ncor = 0;
    my $nsub = 0;
    my $ndel = 0;
    my $nins = 0;
    my $ncli = 0;
    my $acov = 0;

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
	    $ncli += $num;
	}
	else
	{
	    print STDERR "error: unknown operator in CIGAR string '$c'\n\n";
	    exit;
	}
	
	$num = "";
    }
    
    for( my $i = 0; $i < length($s1); $i++ )
    {
	my $rb = substr($s1,$i,1);
	my $cb = substr($s2,$i,1);
	
	if( $rb ne "-" && $cb ne "-" )
	{
	    if( $rb eq $cb )
	    {
		$ncor++;
	    }
	    else
	    {
		$nsub++;
	    }
	}
	elsif( $rb ne "-" && $cb eq "-" )
	{
	    $ndel++;
	}
	elsif( $rb eq "-" && $cb ne "-" )
	{
	    $nins++;
	}

	if( $rb ne "-" )
	{
	    $acov++;
	}
    }

    my ($movie,$zmw,$other) = split(/\//,$qname);

    print join(",",
	       $movie,
	       $zmw,
	       $rname,
	       $flag,
	       $mapq,
 	       $pos,
	       $pos - 1 + $acov,
	       length($s1),
	       $acov,
	       $ncor,
	       $nsub,
	       $ndel,
	       $nins,
	       $ncli,
	       counta($cigar,"D"),
	       counta($cigar,"I"),
	), "\n";
}

################################################################################
#                             SUBROUTINES                                      #
################################################################################

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

sub counta {
    my ($str,$chr) = @_;

    my $num = 0;

    for( my $i = 0; $i < length($str); $i++ )
    {
	if( substr($str,$i,1) eq $chr )
	{
	    $num++;
	}
    }

    return $num;
}
