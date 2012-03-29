#!/usr/bin/perl -w
use strict;
$main::mlab_dir = "../../src_release";

my %undocumented = ("mult",1,"im2col_sliding",1);


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
    while(<IN>) {
	chomp;
	s/mex$mlab_prog/$myprog/g;
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
	    # remove last empty lines
	    while($i >= 0) {
		($$tmp[$i] =~ /^\s*$/) || last;
		$i--;
	    }
	    $#$tmp = $i;
	    $$doc{$key} = $tmp;
	    $tmp = [($_)];
	    $key = $x;
	    if ($x eq "Author") {
		push(@$tmp,"Julien MAIRAL, 2010 (spams, matlab interface and documentation)");
		push(@$tmp,"        Jean-Paul CHIEZE, 2011 (R interface)");
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
		s/^(\s+)param\.([^\s]+)\s/$1$2: /;
	    } else {
		s/param\.//g;
	    }
	    s/Matlab\s+function\s+pcg/R function solve/;
	    s/Matlab\s+expression\s+XAt[^\s\;]+/R expression/;
	    s/Matlab/R/;
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
    my($f,$modifs) = @_;
    my $inblock = 0;
    my ($tmp,$key,$op);
    open(IN,"<$f") || return;
    while(<IN>) {
	chomp;
	(/^\s*$/) && next;
	(/^\s*\#/) && next;
	if($inblock) {
	    if(/^end/) {
		$inblock = 0;
		$$modifs{$key} = { 'op' => $op, 'data' => $tmp};
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

# try to split Description into short description end detail
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
	    } elsif ( $op eq "addlast") {
		push(@$lst,@$tmp);
		$$doc{$key} = $lst;
	    } else {
		print "Unknown op $op\n";
	    }
	}
    }
}

sub prepare_doc {
    my($mlab_prog,$myprog,$doc,$format) = @_;
    my $f = "$main::mlab_dir/mex$mlab_prog.m";
    my $fref = "./refman/$myprog.in";
    my %modifs = ();
    get_doc($f,1,$mlab_prog,$myprog,$doc) || return;
    get_modifs($fref,\%modifs);
    # apply modifs
    apply_modifs($doc,$format,\%modifs);
    split_description($doc);
}

    
1;
