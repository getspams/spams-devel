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

# in : function name
# out : matlab name (mex...)
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

# in : matlab name
# out : function name 
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
# read .R or .py file containing function definitions
# in : $file = file path
#     $find_func : sub to find the function definition
#     $spams (array) = lines of input file
# out:   $progdefs (hash)  key = progname, val = (i1,i2) indexes of first and last line of a function def
#        $progs (hash)  = function name by function line index
sub read_spams {
    my($file,$find_func,$progdefs,$progs,$spams) = @_;
    open(IN,"<$file") || die "$file open err $!\n";
    my $i = -1;
    my $prog = "";
    my $found = 0;
    my $name = "";
    my $i1 = 0;
    while(<IN>) {
	chomp;
	$i++;
	($prog,$found) = &$find_func("$_",$found);
	if($found) {
	    if("$prog") {
		if(! ($prog =~ /^_/) && ! defined($undocumented{$prog})) {
		    $name = $prog;
		    $i1 = $i;
		} else {$found = 0;}
	    }
	    if($found > 0) {
		$$progdefs{$name} = [($i1,$i)];
		$$progs{$i} = $name;
		$found = 0;
	    }
	}
	push(@$spams,$_);
    }
    close(IN);
}

# read a .m doc file
# in : $f = file path
#      $r_mode (bool) = true if R
#      $mlag_prog = name of matlab function
#      $myprog = name of non matlab function
# out : $doc (hash) = arrays of lines by doc section 
sub get_doc {
    my($f,$r_mode,$mlab_prog,$myprog,$doc) = @_;
    if(! open(IN,"<$f") ) {
	print "ERR $f open err $!\n";
	exit 1;
	return 0;
    }
    my $stat = 0;
    my $tmp = [()];
    my $key = "";
    my $prefix = $r_mode ? "spams." : "";
    my $lang = $r_mode ? "R" : "python";
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
		$$tmp[0] =~ s/$/ (spams, matlab interface and documentation)/;
		push(@$tmp,"      Jean-Paul CHIEZE 2011-2012 ($lang interface)");
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
	    s/Matlab\s+function\s+pcg/$lang function solve/;
	    s/Matlab\s+expression\s+XAt[^\s\;]+/$lang expression/;
	    s/Matlab/$lang/;
	    if($key eq "Usage") {
		s/\s*=\s*/ <- /;
	    }
	    push(@$tmp,$_);
	}
    }
    close(IN);
    1;
}

sub get_def {
    my ($spams,$progdefs,$myprog,$idt) = @_;
    my $x;
    my @def = ();
    if(! defined($$progdefs{$myprog})) {
	print STDERR "WARNING! No def for <$myprog>\n";
	return ();
    }
    my $ix = $$progdefs{$myprog};
    for(my $i = $$ix[0];$i <= $$ix[1];$i++) {
	$x = $$spams[$i];
	$x =~ s/^\s+//;$x =~ s/\s+$//;
	push(@def,$x);
    }
    $def[0] =~ s/^[^\(]+\(//;
    $def[$#def] =~ s/\)[^\(]*$/\)/;
    $x = join('',@def);
    my @tmp = split(/\s*,\s*/,$x);
    my $lgr = $idt;
    $#def = -1;
    my @line = ();
    foreach $x (@tmp) {
	my $n = length($x);
	if(($lgr + $n) > 90) {
	    push(@def,join(",",@line) . ",");
	    @line = ($x);
	    $lgr = $idt + $n + 1;
	} else {
	    push(@line,$x);
	    $lgr += $n + 1;
	}
    }
    if($#line >= 0) {
	push(@def,join(",",@line) . ")");
    } else {
	$def[$#def] =~ s/,$/)/;
    }
    @def;
}

# reads files modifying original doc
# in : $r_mode = true if R
#      $f = data file
# $spams : array of python or R file
# $progdefs (hash)  key = progname, val = (i1,i2) indexes of first and last line of a function def
# special replacement of %DEF% in Usage section
# 
# out : $modifs (hash) = modifications by doc section
sub get_modifs {
    my($r_mode,$f,$myprog,$modifs,$spams,$progdefs) = @_;
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
	    if($key eq 'Usage' && ! $r_mode) {
		s/<-/=/;
	    }
	    if(/%DEF%/) {
		my $name = $r_mode ? "spams.$myprog" : $myprog;
		s/%DEF%.*$/spams.$myprog\(/;
		my $n = length($_);
		my $idt = " " x $n;
		my @def = get_def($spams,$progdefs,$myprog,$n);
		$_ .= shift(@def);
		push(@$tmp,$_);
		foreach $_ (@def) {
		    push(@$tmp,"$idt$_");
		}
	    } else {
		push(@$tmp,$_);
	    }
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
# modify $doc according to $modifs
# 
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
# $spams : array of python or R file
# $progdefs (hash)  key = progname, val = (i1,i2) indexes of first and last line of a function def
# Out : $doc of the function
sub prepare_doc {
    my($r_mode,$mlab_prog,$myprog,$doc,$format,$spams,$progdefs) = @_;
    my $f = "$main::mlab_dir/mex$mlab_prog.m";
    my $fref = "./refman/$myprog.in";
    my %modifs = ();
    get_doc($f,$r_mode,$mlab_prog,$myprog,$doc) || return;
    split_description($doc);

    get_modifs($r_mode,$fref,$myprog,\%modifs,$spams,$progdefs);
    # apply modifs
    apply_modifs($doc,$format,\%modifs);
}

############## latex output ###############
@main::texkeys = ('Name','Usage','Description','Inputs', 'Output','Author','Note','Examples');

$main::tex_docformat = {
    'Name' => {},
    'Description' => {},
    'Usage' => {},
    'Inputs' => {},
    'Output' => {'indent' => 1},
    'Author' => {'tag' => 'Authors'},
    'Note' => {'tag' => 'Note', 'optional' => 1},
    'Examples' => {'tag' => 'Examples', 'optional' => 1},
};

# in : $indx : array of last line of function def
sub make_tex_doc {
    my($r_mode,$dir,$indx,$progs,$spams,$progdefs) = @_;
    foreach my $i (@$indx) {
	my $myprog = $$progs{$i};
	my $mlab_prog = mlab_name($myprog);
	my %doc = ();
	prepare_doc($r_mode,$mlab_prog,$myprog,\%doc,$main::tex_docformat,$spams,$progdefs);
	print "++ $myprog\n";
	write_tex_man("$dir/functions/$myprog.in",$myprog,\%doc,$main::tex_docformat);
    }
    modif_tex_src($r_mode,$dir);
}
sub merge_doc {
    my($doc,$tagin,$tagout) = @_;

    if(defined($$doc{$tagin})) {
	my $docp = $$doc{$tagin};
	my $doci = $$doc{$tagout};
	push(@$doci,@$docp);
	undef($$doc{$tagin});
    }
}  
sub write_tex_man {
    my ($f,$prog,$doc,$rdformat) = @_;
    my($rdf,$i,$key,$tmp);
    merge_doc($doc,'Param','Inputs');
    merge_doc($doc,'detail','Description');
    open(OUT,">$f") || die "$f create err $!\n";
    print OUT "#\n";
    foreach $key (@main::texkeys) {
	(defined($$rdformat{$key})) || next;
	$rdf = $$rdformat{$key};
	if(! defined($$doc{$key})) {
	    if(! defined($$rdf{'optional'})) {
		print STDERR "!! $prog : $key MISSING.\n";
	    }
	    next;
	}
	$tmp = $$doc{$key};
	my $tag = (defined($$rdf{'tag'})) ? $$rdf{'tag'} : $key;
	print OUT "# $tag: ";
	if(defined($$rdf{'prog'})) {
	    my $func = $$rdf{'prog'};
	    my @res = &$func($tmp);
	    print OUT join("\n",@res), "\n}\n";
	} else {
	    print OUT join("\n#    ",@$tmp), "\n#\n";
	}
    }
    
    close(OUT);
}

sub modif_tex_src {
    my ($r_mode,$dir) = @_;
    open(IN,"<../../doc/doc_spams.tex") || die "Cannot read ../../doc/doc_spams.tex: $!\n";
    open(OUT,">$dir/doc_spams.tex") || die "Cannot create $dir/doc_spams.tex\n";
    while(<IN>) {
	chomp;
	if(/mex([A-Z][_A-z]+)/) {
	    my $s = $1;
	    my $x = $s;
	    $x =~ s/\\_/_/;
	    $s =~ s/\\/\\\\/;
	    if (! defined($main::conv_names{$x}) ) {
		print STDERR "Unkown $s\n";
	    } else {
		my $s1 = $main::conv_names{$x};
		$x = $s1;
		$s1 =~ s/_/\\_/;
		s/mex$s/spams.$s1/g;
	    }
	    if(/verbatiminput/) {
	      s|../src_release/spams\.(.*)\.m|functions/$1.in|;
	      s/\\_/_/;
	      if (! -r "$dir/functions/$x.in") {
		  $_ = "\\verbatiminput{functions/missing.in}";
	      }
	    }
	}
	print OUT "$_\n";
    }
    close(IN);
    close(OUT);
}
1;
