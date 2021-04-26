#!/usr/bin/env perl
use strict;

my @HEADFLD = qw/chr start end type commonSome commonAll dbsnpRef dbsnpAlt allRef allAlt minRef minAlt id/;

while (my $line = <STDIN>)
{
    chomp $line;
    my @f = split(/\t/,$line);

    # get max MAF
    my @allMafs = uniqVals($f[9]);
    my $maxMaf = maxval(@allMafs);

    # get common alleles
    my @allRef = uniqVals($f[10]);
    my @allAlt = uniqVals($f[11]);
    if (0)
    {
	print STDERR sprintf("refAlleles has more than one:\t%s\t%s\n",$f[10],makeList(@allRef))
	    if (scalar(@allRef)>1);
	print STDERR sprintf("altAlleles has more than one:\t%s\t%s\n",$f[11],makeList(@allAlt))
	    if (scalar(@allAlt)>1);
    }

    my @dbsnpRef = uniqVals($f[4]);
    my @dbsnpAlt = uniqVals($f[6]);
    my $usedbsnp = ( (scalar(@dbsnpAlt)<scalar(@allAlt)) || (scalar(@dbsnpRef)<scalar(@allRef)) );
    my @minRef = ($usedbsnp) ? @dbsnpRef : @allRef;
    my @minAlt = ($usedbsnp) ? @dbsnpAlt : @allAlt;

    # Notes
    my $notes = $f[14];
    my $commonSome = ($notes =~ /commonSome/);
    my $commonAll = ($notes =~ /commonAll/);
    
    my @newf = (@f[0],$f[1]+1,$f[2],@f[13], # Coords are 1-based inclusive for bcftools annotate
		$commonSome || "0", $commonAll || "0",
		makeList(@minRef), makeList(@minAlt),
		makeList(@dbsnpRef), makeList(@dbsnpAlt), makeList(@allRef),makeList(@allAlt),
		@f[3],$maxMaf);

    print join("\t",@newf)."\n";
}



# - - - Funcs

sub maxval
{
    my (@vals) = @_;

    my $maxv = -inf;
    foreach my $val (@vals)
    {
	$maxv = $val if ($val>$maxv);
    }
    return $maxv;
}

sub makeList
{
    my (@flds) = @_;
    return join(",",@flds);
}

sub uniqVals
{
    my ($fld) = @_;

    my @allAlleles = split(/,/,$fld);
    my $outAllelesH = {};
    foreach my $allele (@allAlleles) { $outAllelesH->{$allele}++ if ($allele); }
    my @outAlleles = keys(%$outAllelesH);
    
#    print STDERR join("\t",$fld,@outAlleles)."\n";

    return @outAlleles;
}
