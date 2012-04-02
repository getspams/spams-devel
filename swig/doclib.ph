#!/usr/bin/perl -w
use strict;
$main::mlab_dir = "../../src_release";

my %undocumented = ("mult",1,"im2col_sliding",1);

%main::conv_names = (
    "Sort", "sort",
    "CalcAAt", "calcAAt",
    "CalcXAt", "calcXAt",
    "CalcXY", "calcXY",
    "CalcXYt", "calcXYt",
    "CalcXtY", "calcXtY",
    "Bayer", "bayer",
    "ConjGrad", "conjGrad",
    "InvSym", "invSym",
    "Normalize", "normalize",
    "SparseProject", "sparseProject",
    "Lasso", "lasso",
    "FistaFlat", "fistaFlat",
    "ProximalFlat", "proximalFlat",
    "TrainDL", "trainDL",
    "TrainDL_Memory", "trainDL_Memory",
    "CD", "cd",
    "LassoMask", "lassoMask",
    "LassoWeighted", "lassoWeighted",
    "OMP", "omp",
    "OMPMask", "ompMask",
    "SOMP", "somp",
    "FistaGraph", "fistaGraph",
    "FistaTree", "fistaTree",
    "ProximalGraph", "proximalGraph",
    "ProximalTree", "proximalTree",
);

%main::tomlab = ();
while( my($k,$v) = each(%main::conv_names)) {
    $main::tomlab{$v} = $k;
}

sub mlab_name {
    my($name) = @_;
    my $mlab_name = "";
    if(defined($main::tomlab{$name})) {
	$mlab_name = $main::tomlab{$name};
    } else {
	print STDERR "Cannot convert name <$name>\n";
	$mlab_name = $name;
    }
    $mlab_name;
}

sub newname {
    my($mlab_name) = @_;
    my $name = "";
    if(defined($main::conv_names{$mlab_name})) {
	$name = $main::conv_names{$mlab_name};
    } else {
	print STDERR "Cannot convert name <$mlab_name>\n";
	$name = $mlab_name;
    }
    $name;
}
sub read_spams {
    my($file,$find_func,$indx,$progs,$spams) = @_;
    open(IN,"<$file") || die "$file open err $!\n";
    my $i = -1;
    my $prog = "";
    my $found = 0;
    my $name = "";
    while(<IN>) {
	chomp;
	$i++;
	($prog,$found) = &$find_func("$_",$found);
	if($found) {
	    if("$prog") {
		if(! ($prog =~ /^_/) && ! defined($undocumented{$prog})) {
		    $name = $prog;
		} else {$found = 0;}
	    }
	    if($found > 0) {
		push(@$indx,$i);
		$$progs{$i} = $name;
		$found = 0;
	    }
	}
	push(@$spams,$_);
    }
    close(IN);
}

sub get_doc {
    my($f,$r_mode,$mlab_prog,$myprog,$doc) = @_;
    if(! open(IN,"<$f") ) {
	print "ERR $f open err $!\n";
	return 0;
    }
    my $stat = 0;
    my $tmp = [()];
    my $key = "";
    my $prefix = $r_mode ? "spams." : "";
    while(<IN>) {
	chomp;
	if(/mex([A-Z][_A-z]+)/) {
	    my $s = $1;
	    (defined($main::conv_names{$s}) ) || die "Inconnu : $s\n";
	    my $s1 = $main::conv_names{$s};
	    s/mex$s/$prefix$s1/g;
	}
	if(! $stat) {
	    (s/^%\s*Usage\s*:\s*//) || next;
	    $stat = 1;
	    if($r_mode) {s/\s*=\s*/ <- /;}
	    push(@$tmp,$_);
	    $key = 'Usage';
	    next;
	}
	if(s/^%\s([^\s:]+)\s*:\s*//) {
	    my $x = $1;
	    my $i = $#$tmp;
#	    # remove last empty lines
#	    while($i >= 0) {
#		($$tmp[$i] =~ /^\s*$/) || last;
#		$i--;
#	    }
#	    $#$tmp = $i;
	    $$doc{$key} = $tmp;
	    $tmp = [($_)];
	    $key = $x;
	    if ($x eq "Author") {
		$tmp = [()];
		push(@$tmp,"Julien MAIRAL (spams, matlab interface and documentation);");
		push(@$tmp,"Jean-Paul CHIEZE (R interface)");
		push(@$tmp,"");
		$$doc{$x} = $tmp;
		last;
	    }
	} else {
	    s/^%\s?//;
	    if(/^\s*param:\s*struct/) {
		$$doc{$key} = $tmp;
		$tmp = [()];
		$key = "Param";
		next;
	    }
	    s/param\.lambda([^\w])/param.lambda1$1/;
	    s/(param\.[^\s:]+)\s*:/$1/;
	    if($key eq "Param") {
		s/^\s*param\.([\w]+)\s*,\s*param\.([\w]+)\s*/    $1, $2: /;
		s/^\s*param\.([^\s]+)\s/    $1: /;
	    }
	    s/param\.//g;
	    if($r_mode) {
		s/Matlab\s+function\s+pcg/R function solve/;
		s/Matlab\s+expression\s+XAt[^\s\;]+/R expression/;
		s/Matlab/R/;
	    } else {
		s/Matlab\s+function\s+pcg/python function solve/;
		s/Matlab\s+expression\s+XAt[^\s\;]+/python expression/;
		s/Matlab/python/;
	    }
	    if($key eq "Usage") {
		s/\s*=\s*/ <- /;
	    }
	    push(@$tmp,$_);
	}
    }
    close(IN);
    1;
}

sub get_modifs {
    my($r_mode,$f,$modifs) = @_;
    my $inblock = 0;
    my ($tmp,$key,$op);
    my $expr = "";
    open(IN,"<$f") || return;
    while(<IN>) {
	chomp;
	(/^\s*$/) && next;
	(/^\s*\#/) && next;
	if(s/^\[([PR])\]//) {  # this line is only for R or python
	    my $x = ($1 eq "P") ? 0 : 1;
	    ($x == $r_mode) || next;
	}
	if(! $r_mode) { s/<-/=/;}
	if($inblock) {
	    if(/^end/) {
		$inblock = 0;
		$$modifs{$key} = { 'op' => $op, 'data' => $tmp};
		if("$expr") {
		    my $x = $$modifs{$key};
		    $$x{'subst'} = $expr;
		    $expr = "";
		}
		next;
	    }
	    push(@$tmp,$_);
	} else {
	    (/^begin\s+([^\s]+)\s+([^\s]+)$/) || next;
	    $op = $1;
	    $key = $2;
	    $tmp = [()];
	    $inblock = 1;
	}
    }
    close(IN);
}

# try to split Description into short description and detail
sub split_description {
    my($doc) = @_;
    my $tmp = $$doc{'Description'};
    my $det = [()];
    ($#$tmp < 3) && return;
    for(my $i = 0;$i <= $#$tmp;$i++) {
	my $s = $$tmp[$i];
	if(($s =~ /^\s*$/) || ($s =~ /\.$/)) {
	    my $j = $i;
	    $i++;
	    while($i <= $#$tmp) {
		push(@$det,$$tmp[$i++]);
	    }
	    $$doc{'detail'} = $det;
	    $#$tmp = $j;
	    last;
	}
		
    }
}

sub apply_modifs {
    my($doc,$format,$modifs) = @_;
    my($op,$tmp);
    while(my ($key,$x) = each(%$modifs)) {
	$op = $$x{'op'};
	$tmp = $$x{'data'};
	if($op eq "repl") {
	    $$doc{$key} = $tmp;
	
	} else {
	    (defined($$format{$key})) || next;
	    my $lst = (defined($$doc{$key})) ? $$doc{$key} : [()];
	    if ( $op eq "addfirst") {
		die "Addfirst not implemented\n";
	    } elsif ( $op eq "addlast") {
		push(@$lst,@$tmp);
		$$doc{$key} = $lst;
	    } elsif ( $op eq "subst") {
		my $e = $$tmp[0];
		$e =~ s/^\s+//;
		for(my$i = 0;$i <= $#$lst;$i++) {
		    my $s = $$lst[$i];
		    eval("\$s =~ $e");
		    $$lst[$i] = $s;
		}
	    } else {
		die "Unknown op $op\n";
	    }
	}
    }
}

# IN: $r_mode = 0/1 for python/R
# $mlab_prog = prog name in matlab
# $myprog = prog name for python or R
# $doc : empty doc hash table
# $format : hash table describing format of the different parts of doc
# Out : $doc of the function
sub prepare_doc {
    my($r_mode,$mlab_prog,$myprog,$doc,$format) = @_;
    my $f = "$main::mlab_dir/mex$mlab_prog.m";
    my $fref = "./refman/$myprog.in";
    my %modifs = ();
    get_doc($f,$r_mode,$mlab_prog,$myprog,$doc) || return;
    split_description($doc);
    get_modifs($r_mode,$fref,\%modifs);
    # apply modifs
    apply_modifs($doc,$format,\%modifs);
}

    
1;
