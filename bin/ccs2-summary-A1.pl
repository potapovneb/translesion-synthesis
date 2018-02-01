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

### command-line options
my $o_qv  = 0;
my $o_np  = 1;
my $o_mapq = 254;
my $o_rlen_min_a = -1;
my $o_rlen_max_a = -1;
my $o_rlen_min_f = -1;
my $o_rlen_max_f = -1;
my $o_zm  = -1;
my $o_sna = 1e6;
my $o_snc = 1e6;
my $o_sng = 1e6;
my $o_snt = 20.0;
my $o_lb  = -1;
my $o_ub  = -1;
my $o_log = "";
my $o_wl = "";
my $o_expand_deletion = 1;
my $o_mutation_dir = "mutation-bwa-pacbio";
my $o_testfile = "";

GetOptions(
    "qv=f"  => \$o_qv,
    "np=f"  => \$o_np,
    "mapq=f"  => \$o_mapq,
    "rlen-min-a=f" => \$o_rlen_min_a,
    "rlen-max-a=f" => \$o_rlen_max_a,
    "rlen-min-f=f" => \$o_rlen_min_f,
    "rlen-max-f=f" => \$o_rlen_max_f,
    "zm=f"  => \$o_zm,
    "sna=f" => \$o_sna,
    "snc=f" => \$o_snc,
    "sng=f" => \$o_sng,
    "snt=f" => \$o_snt,
    "lb=f"  => \$o_lb,
    "ub=f"  => \$o_ub,
    "log"   => \$o_log,
    "wl=s"  => \$o_wl,
    "expand-deletion!" => \$o_expand_deletion,
    "mutation-dir=s" => \$o_mutation_dir,
    "testfile=s" => \$o_testfile,
    );

### command-line arguments
if( @ARGV == 0 )
{
    print "usage: $0 [options] reference.fasta sampleDir1 SampleDir2 ...\n";
    print "\n";
    print "options:\n";
    print " --qv\t\tminimum quality value ($o_qv)\n";
    print " --np\t\tminimum number of passes ($o_np)\n";
    print " --zm\t\tanalyse a specific ZMW ($o_zm)\n";
    print " --sna\t\tmaximum SNR-A ($o_sna)\n";
    print " --snc\t\tmaximum SNR-C ($o_snc)\n";
    print " --sng\t\tmaximum SNR-G ($o_sng)\n";
    print " --snt\t\tmaximum SNR-T ($o_snt)\n";
    print " --lb\t\tskip first N bases ($o_lb)\n";
    print " --ub\t\tskip last N bases ($o_ub)\n";
    print " --log\t\tprint log info ($o_log)\n";
    print " --wl\t\twhitelist file ($o_wl)\n";
    print " --expand-deletion\t($o_expand_deletion)\n";
    print "\n";
    exit;
}

my $reference = shift @ARGV;
my @sampleDirs = @ARGV;

### load reference data
my $refs = load_references($reference);
my $hmap = homopolymer($refs);

### init counts
my @type = qw(AA AC AT AG CA CC CT CG TA TC TT TG GA GC GT GG Deletion Insertion D1 DX I1 IX);
my %count = map { $_ => 0 } @type;

### global storage for filtered data
my %DATA = ();

my $testlist = ( $o_testfile ne "" ) ? load_blacklist($o_testfile) : {};

### process samples one-by-one
for my $sampleDir ( @sampleDirs )
{
    my $variants_file  = sprintf( "%s/%s/variants.csv.7z",       $sampleDir, $o_mutation_dir );
    my $zmws_file      = sprintf( "%s/%s/zmws.csv.7z",           $sampleDir, $o_mutation_dir );
    my $alns_file      = sprintf( "%s/%s/aln.csv.7z",            $sampleDir, $o_mutation_dir );
    my $chimeric_file  = sprintf( "%s/chimeric/chimeric.csv.7z", $sampleDir );
    my $whitelist_file = sprintf( "%s/%s/whitelist.csv.7z",      $sampleDir, $o_mutation_dir );

    if( $o_log )
    {
	print STDERR $variants_file, "\n";
	print STDERR $chimeric_file, "\n";
	print STDERR $whitelist_file, "\n";
	print STDERR $zmws_file, "\n";
	print STDERR $alns_file, "\n";
    }

    ### load pacbio read stats
    my $zmws = load_zmw_data($zmws_file);
    
    ### load alignment stats
    my $alns = load_zmw_data($alns_file);

    ### load chimeric reads
    my $blacklist = load_blacklist($chimeric_file);

    ### load whitelist
    my $whitelist = ( $o_wl ne "" ? load_blacklist($whitelist_file) : {} );


    ### filter variants
    filter($variants_file, $zmws, $alns, $blacklist, $whitelist);
}

### print column headers
printout_data(\%DATA);

sub filter {
    my ($variants_file,$zmws,$alns,$blacklist,$whitelist) = @_;

    my @head = ();
    
    open(VARIANTS,"7za e -so $variants_file |") || die "Can't open '$variants_file': $!";
    
    while( my $line = <VARIANTS> )
    {
	chomp($line);
	
	### parse line
	my @tokens = split(/,/,$line);

	my %entry = ();
	
	if( @head == 0 )
	{
	    @head = @tokens;
	    next;
	}
	else
	{
	    for( my $i = 0; $i < @tokens; $i++ )
	    {
		$entry{$head[$i]} = $tokens[$i];
	    }
	}
	
	if( ! exists $$refs{$entry{"Ref"}} )
	{
	    print STDERR "error :: missing reference '$entry{'Ref'}'\n\n";
	    exit;
	}
	
	### filter
	my $refLength = length($$refs{$entry{"Ref"}});
	
	### frequently used vars
	my $movie = $entry{"Movie"};
	my $zmw = $entry{"ZMW"};

	### filter by number of passes
	next if( $$zmws{$movie}{$zmw}{"NP"} < $o_np );
	
	### filter by mapping quality
	next if( $o_mapq != -1 && $$alns{$movie}{$zmw}{"MAPQ"} < $o_mapq );

	### filter by read length
	next if( $o_rlen_min_a != -1 && $$zmws{$movie}{$zmw}{"ReadLength"} < $refLength - $o_rlen_min_a );
	next if( $o_rlen_max_a != -1 && $$zmws{$movie}{$zmw}{"ReadLength"} > $refLength + $o_rlen_max_a );
	next if( $o_rlen_min_f != -1 && $$zmws{$movie}{$zmw}{"ReadLength"} < $refLength * $o_rlen_min_f );
	next if( $o_rlen_max_f != -1 && $$zmws{$movie}{$zmw}{"ReadLength"} > $refLength * $o_rlen_max_f );

	### make sure read covers the entire region of interest (in other words, avoid truncated reads)
	next if( $o_lb != -1 && $$alns{$movie}{$zmw}{"AlnStart"} > $o_lb );
	next if( $o_ub != -1 && $$alns{$movie}{$zmw}{"AlnEnd"} < $refLength - $o_ub );

	### filter by SNR
	next if( $$zmws{$movie}{$zmw}{"SnrA"} > $o_sna );
	next if( $$zmws{$movie}{$zmw}{"SnrC"} > $o_snc );
	next if( $$zmws{$movie}{$zmw}{"SnrG"} > $o_sng );
	next if( $$zmws{$movie}{$zmw}{"SnrT"} > $o_snt );
	
	### exclude blacklisted reads (sample-specific blacklist)
	next if( %$blacklist && exists $$blacklist{$movie}{$zmw} );
	
	### exclude NOT whitelisted reads (sample-specific whitelist)
	next if( %$whitelist && ! exists $$whitelist{$movie}{$zmw} );
	
	### exclude NOT whitelisted reads (global whitelist)
	next if( %$testlist && ! exists $$testlist{$movie}{$zmw} );

	### ----- custom block (begin) -----------------------------------------
	
	process_entry(\%entry,\%DATA);

	### ----- custom block (end) -------------------------------------------
	
	### print out data lines (optional)
	if( $o_log )
	{
	    print STDERR $line, "\n";
	}
    }
    
    close(VARIANTS);
}

sub load_references {
    my ($file) = @_;
    
    my %refs = ();
    
    my $name = "";
    my $seq = "";
    
    open(FA,$file) || die "Can't open '$file': $!";
    
    while( my $line = <FA> )
    {
	chomp($line);
	
	if( substr($line,0,1) eq ">" )
	{
	    if( $name ne "" )
	    {
		$refs{$name} = uc($seq);
	    }
	    
	    $name = substr($line,1);
	    $seq = "";
	}
	else
	{
	    $seq .= $line;
	}
    }
    
    if( $name ne "" )
    {
	$refs{$name} = uc($seq);
    }
    
    close(FA);
    
    return \%refs;
}

sub load_blacklist {
    my ($file) = @_;

    my %data = ();

    if( substr($file,-3) eq ".7z" )
    {
	open(SR,"7za e -so $file |") || die "Can't open '$file': $!";
    }
    elsif( substr($file,-3) eq ".gz" )
    {
	open(SR,"gzip -cd $file |") || die "Can't open '$file': $!";
    }
    elsif( substr($file,-3) eq ".bz2" )
    {
	open(SR,"bzip2 -cd $file |") || die "Can't open '$file': $!";
    }
    else
    {
	open(SR,$file) || die "Can't open '$file': $!";
    }

    while( my $line = <SR> )
    {
	chomp($line);

	my ($movie,$zmw) = split(/,/,$line);

	$data{$movie}{$zmw} = 1;
    }

    close(SR);

    return \%data;
}

sub load_zmw_data {
    my ($file) = @_;
    
    my %zmws = ();
    my @head = ();

    open(IN,"7za x -so $file |") || die "Can't open '$file': $!";
    
    while( my $line = <IN> )
    {
	chomp($line);
	
	### parse line
	my @tokens = split(/,/,$line);
	
	if( @head == 0 )
	{
	    ### extract column names
	    @head = @tokens;
	}
	else
	{
	    ### store data fields
	    my %entry = ();
	    
	    for( my $i = 0; $i < @head; $i++ )
	    {
		$entry{$head[$i]} = $tokens[$i];
	    }

	    $zmws{$entry{"Movie"}}{$entry{"ZMW"}} = \%entry;
	}
    }

    close(IN);

    return \%zmws;
}

sub complement {
    my ($seq) = @_;

    my %bases = (
	"A" => "T",
	"C" => "G",
	"G" => "C",
	"T" => "A",
	"a" => "t",
	"c" => "g",
	"g" => "c",
	"t" => "a",
	"-" => "-",
	);

    my $complement = join("",map { $bases{$_} } split(//,$seq));

    return $complement;
}

sub process_entry {
    my ($entry,$data) = @_;
    
    if( $$entry{"Type"} ne "SNP" || $$entry{"RefBP"} ne $$entry{"AltBP"} )
    {
	print STDERR join( ",",
			   $$entry{"Movie"},
			   $$entry{"ZMW"},
			   $$entry{"Ref"},
			   $$entry{"Pos"},
			   $$entry{"Type"},
			   $$entry{"QV"},
			   $$entry{"RefBP"},
			   $$entry{"AltBP"},
			   $$entry{"IndelType"},
			   $$entry{"Bases"}), "\n";
    }

    my $qt = ( $$entry{"QV"} < $o_qv ) ? "LQ" : "HQ";
    
    if( $$entry{"Type"} eq "SNP" )
    {
	my $rb = uc($$entry{"RefBP"});
	my $cb = uc($$entry{"AltBP"});
	
	if( $rb eq $cb )
	{
	    $$data{$$entry{"Pos"}}{$qt}{"Correct"}++;
	}
	else
	{
	    $$data{$$entry{"Pos"}}{$qt}{"Substitution"}++;
	    $$data{$$entry{"Pos"}}{$qt}{"SubSpectrum"}{$rb.$cb}++;
	}
    }
    else
    {
	if( $$entry{"IndelType"} eq "Insertion" )
	{
	    ### count indel occurrences at a given position not indel length
	    $$data{$$entry{"Pos"}}{$qt}{"Insertion"}++;
	    $$data{$$entry{"Pos"}}{$qt}{"InsSpectrum"}{$$entry{"Bases"}}++;
	}
	else
	{
	    if( $o_expand_deletion )
	    {
		$$data{$$entry{"Pos"}}{$qt}{"DelSpectrum"}{$$entry{"Bases"}}++;
		
		### deletion is attached to a particular position but can span several positions
		for( my $i = 1; $i <= $$entry{"Length"}; $i++ )
		{
		    ### update deletion count for all spanned positions
		    $$data{$$entry{"Pos"}+$i}{$qt}{"Deletion"}++;

		    ### single-base deletions
		    if( $$entry{"Length"} == 1 )
		    {
			$$data{$$entry{"Pos"}+$i}{$qt}{"D1"}++;
		    }
		}
	    }
	    else
	    {
		### count indel occurrences at a given position not indel length
		$$data{$$entry{"Pos"}}{$qt}{"Deletion"}++;
		$$data{$$entry{"Pos"}}{$qt}{"DelSpectrum"}{$$entry{"Bases"}}++;
	    }
	}
    }
}

sub printout_data {
    my ($data) = @_;
    
    ### print resulting headers
    print join( ",",
		"Position",
		"RefBP",
		"HPLen",
		"Correct (HQ)",
		"Substitution (HQ)",
		"Deletion (HQ)",
		"Insertion (HQ)",
		"D1 (HQ)",
		"DX (HQ)",
		"I1 (HQ)",
		"IX (HQ)",
		"SubSpectrum (HQ)",
		"DelSpectrum (HQ)",
		"InsSpectrum (HQ)",
		"Correct (LQ)",
		"Substitution (LQ)",
		"Deletion (LQ)",
		"Insertion (LQ)",
		"D1 (LQ)",
		"DX (LQ)",
		"I1 (LQ)",
		"IX (LQ)",
		"SubSpectrum (LQ)",
		"DelSpectrum (LQ)",
		"InsSpectrum (LQ)",
	), "\n";

    my @rname = keys %$refs;

    for my $pos ( sort { $a <=> $b } keys %$data )
    {
	my @values = ( $pos,
		       substr($$refs{$rname[0]},$pos-1,1),
		       $$hmap{$rname[0]}[$pos-1] );
	
	for my $qt ( "HQ", "LQ" )
	{
	    ### counts
	    for my $type ( "Correct", "Substitution", "Deletion", "Insertion", "D1", "DX", "I1", "IX" )
	    {
		my $count = 0;
		
		if( exists $$data{$pos}{$qt}{$type} )
		{
		    $count = $$data{$pos}{$qt}{$type};
		}
		
		push @values, $count;
	    }
	    
	    ### spectrum
	    for my $spectrum ( "SubSpectrum", "DelSpectrum", "InsSpectrum" )
	    {
		my $spec = "";
		
		if( exists $$data{$pos}{$qt}{$spectrum} )
		{
		    ### spectrum
		    my $list = $$data{$pos}{$qt}{$spectrum};
		    
		    ### sort by frequency
		    my @sorted = sort { $$list{$b} <=> $$list{$a} } keys %$list;
		    
		    ### build resulting string
		    $spec = join(";", map { sprintf("%s=%s",$_,$$list{$_}) } @sorted);
		}
		
		push @values, $spec;
	    }
	}
	
	print join(",",@values), "\n";
    }
}

sub homopolymer {
    my ($refs) = @_;
    
    my %data = ();
    
    for my $name ( keys %$refs )
    {
	my @seq = split(//,$$refs{$name});
	
	my $count = 1;
	
	for( my $i = 0; $i < @seq; $i++ )
	{
	    if( $i == 0 || $seq[$i-1] ne $seq[$i] )
	    {
		$count = 1;
	    }
	    else
	    {
		$count++;
	    }
	    
	    $data{$name}[$i] = $count;
	}
	
	for( my $i = @seq - 1; $i > 0; $i-- )
	{
	    if( $seq[$i-1] eq $seq[$i] )
	    {
		$data{$name}[$i-1] = $data{$name}[$i];
	    }
	}
    }

    return \%data;
}
