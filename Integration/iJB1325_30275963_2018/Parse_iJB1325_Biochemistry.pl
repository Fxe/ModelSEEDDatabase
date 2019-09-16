#!/usr/bin/env perl
#use PlantSEED::Formulas;
use warnings;
use strict;
my @temp=();

my $Biochemistry = "iJB1325_ATCC1015";
my $Compounds = $Biochemistry."_Compound_Table.txt";
my $Reactions = $Biochemistry."_Reaction_Table.txt";

open(FH, "< $Compounds");
my $header=1;
my %Original_Compounds=();
my $Default_Cpt="";
while(<FH>){
    chomp;
    if($header){$header--;next}
    @temp=split(/\t/,$_);

    ###############################################################
    #Redundant code for finding unicode hexademical code points
    foreach my $temp (@temp){
	my $found=0;
	foreach my $char (split(//,$temp)){
	    if($char =~ /[^\x00-\x7F]+/){
		$found=1;
		my $n = ord $char;
		open my $BYTE, '>:utf8', \ my $bytes;
		print {$BYTE} $char;
#		printf "%s\t%s\t%x\t%b\t%x %x\t %b %b\n", $char, $n, $n, $n, (unpack('CC', $bytes)) x 2;
	    }
	}
#	print $temp,"\n" if $found;
    }
    ###############################################################

    my $cpd_cpt = $temp[0];
    
    #Clean up identifier
    $temp[0] =~ s/^M_+//;

    #Clean up name
    $temp[1] =~ s/\x{ce}\x{b1}/alpha/;
    $temp[1] =~ s/\x{ce}\x{b2}/beta/;
    $temp[1] =~ s/\x{ce}\x{b3}/gamma/;
    $temp[1] =~ s/\x{e4}/e/;
    $temp[1] =~ s/\x{e2}\x{88}\x{92}/-/;

    #Remove Compartment
    my $cpd = $temp[0];
    my $cpt = "c";

    if($cpd =~ s/([a-z])$//){
	$cpt = $1;
    }

    #Define KEGG
    $temp[5] = "" if !$temp[5];

    if($Default_Cpt eq ""){
	$Default_Cpt = $temp[2];
    }

    $Original_Compounds{$cpd_cpt}={'ID'=>$cpd,
				   'NAMES'=>$temp[1],
				   'COMPARTMENT'=>$temp[2],
				   'KEGG'=>$temp[5]};
}
close(FH);

open(FH, "< $Reactions");
$header=1;
my %Original_Reactions=();
my %Cpds_in_Rxns=();
while(<FH>){
    chomp;
    if($header){$header--;next}

    @temp=split(/\t/,$_);

    my $rxn = $temp[0];

    #Clean up identifier
    $rxn =~ s/^R_+//;

    #Clean up unicode character in name
    $temp[1] =~ s/\x{a0}//;

    #Go through reactants
    my ($rev,$reactants,$products)=@temp[2..5];

    #Skipping all boundary exchange reactions
    next if !$products;

    my @reactants = split(/;/,$reactants);
    my @eqn=();
    my %cpts = (); #got to double-check
    my $cpt_count = 0;
    foreach my $entry (@reactants){
	$entry =~ s/\[([-eE\d.]+)\]$//;
	my $coeff = $1;

	my $cpt = $Default_Cpt;
	my $cpd = $entry;
	if(!exists($Original_Compounds{$entry})){
	    print "Warning: cannot find ".$entry."\n";
	}else{
	    $cpt = $Original_Compounds{$entry}{'COMPARTMENT'};
	    $cpd = $Original_Compounds{$entry}{'ID'};
	}

	$Cpds_in_Rxns{$cpd}=1;

	if(defined($coeff) && $coeff != 1){
	    $coeff="(".$coeff.")";
	}else{
	    $coeff=undef($coeff);
	}

	if(!exists($cpts{$cpt})){
	    $cpts{$cpt}=$cpt_count;
	    $cpt_count++;
	}

	$cpt = "[".$cpts{$cpt}."]";
	
	my $rct_str = $cpd.$cpt;
	if(defined($coeff)){
	    $rct_str = $coeff." ".$rct_str;
	}

	push(@eqn,$rct_str);
    }

    my $reversibility = "<=>";
    if($rev ne "False"){
	$reversibility = "=>";
    }
    push(@eqn,$reversibility);

    my @products = split(/;/,$products);
    foreach my $entry (@products){
	$entry =~ s/\[([-eE\d.]+)\]$//;
	my $coeff = $1;

	my $cpt = $Default_Cpt;
	my $cpd = $entry;
	if(!exists($Original_Compounds{$entry})){
	    print "Warning: cannot find ".$entry."\n";
	}else{
	    $cpt = $Original_Compounds{$entry}{'COMPARTMENT'};
	    $cpd = $Original_Compounds{$entry}{'ID'};
	}

	$Cpds_in_Rxns{$cpd}=1;

	if(defined($coeff) && $coeff != 1){
	    $coeff="(".$coeff.")";
	}else{
	    $coeff=undef($coeff);
	}

	if(!exists($cpts{$cpt})){
	    $cpts{$cpt}=$cpt_count;
	    $cpt_count++;
	}

	$cpt = "[".$cpts{$cpt}."]";

	my $pdt_str = $cpd.$cpt;
	if(defined($coeff)){
	    $pdt_str = $coeff." ".$pdt_str;
	}
	push(@eqn,$pdt_str);
    }
    my $eqn_str = join(" ",@eqn);

    $Original_Reactions{$rxn}={'ID'=>$rxn,
			       'NAMES'=>$temp[1],
			       'EQUATION'=>$eqn_str};
}
close(FH);

my $filestub = $Compounds;
$filestub =~ s/_Compound_Table\.txt$//;

open(OUT, "> ".$filestub."_Compounds.tbl");
my @Headers=("ID","NAMES","KEGG","COMPARTMENT");
print OUT join("\t",@Headers),"\n";
foreach my $id (sort keys %Original_Compounds){
    next if !exists($Cpds_in_Rxns{$Original_Compounds{$id}{'ID'}});

    foreach my $h (@Headers){
	print OUT $Original_Compounds{$id}{$h};
	print OUT "\t" unless $h eq $Headers[$#Headers];
    }
    print OUT "\n";
}
close(OUT);

open(OUT, "> ".$filestub."_Reactions.tbl");
@Headers=("ID","NAMES","EQUATION");
print OUT join("\t",@Headers),"\n";
foreach my $id (sort keys %Original_Reactions){
    foreach my $h (@Headers){
	print OUT $Original_Reactions{$id}{$h};
	print OUT "\t" unless $h eq $Headers[$#Headers];
    }
    print OUT "\n";
}
close(OUT);
