#! /usr/bin/perl
 
use Getopt::Long;
 
my $cutoff=5;
my $blacklist='';
my $restrict='';  #the restrict file
my $whitelist='';
GetOptions ('cutoff=s' => \$cutoff,
            'blacklist=s' => \$blacklist, 
            'restrict=s' => \$restrict, 
            'whitelist=s'=> \$whitelist );
 
# get the unique blacklisted functions
%seen = ();
if ( $blacklist ){
    open(BLACKLIST, "< $blacklist") or die "Couldn't open $blacklist for reading: $!\n";
    while(<BLACKLIST>){
	chomp;
	foreach $item (split(/[ ,;]+/)) {
	    $seen{$item}++;
	}
    }
}
@blist = keys %seen;
       
%seen = ();
if( $whitelist ){
    open(WHITELIST, "< whitelist") or die "Couldn't open $whitelist for reading: $!\n";
    while(<WHITELIST>){
	chomp;
	foreach $item (split(/[ ,;]+/)) {
	    $seen{$item}++;
	}
    }
}
@wlist = keys %seen;
 
 
@rest=();
unless($restrict =~ /^$/ ){
    open( RESTRICT, "< $restrict") or die "Couldn't opendir $restrict for reading: $!\n" ;
    while( <RESTRICT> ){
	chomp;
	push( @rest, $_ );
    }
}
 
%calltree = ();
%allfun = ();
 
LINE: while (<>) {
    if (/^sample\.interval=/) {
	# <TODO>
    # do something about the timings
	s/sample\.interval=//;
	$sample = $_ / 1e6;
    # </TODO>
    } else {
	chomp;
    
    ##
    ## whitelist, only keep function that are matching the whitelist
    ##     
	if( $whitelist){
      s/"//g;
      $out='';
      foreach $word ( split(/ / ) ){  
	  foreach $white (@wlist){ 
          $pat = '^'.$white.'$' ;
          $out .= ' "'.$word.'"' if $word =~ /$pat/ ;
	  }
      }
      $_ = $out;
      
	}
    
 
    ##
    ## blacklist, simply remove blacklisted functions from the list
    ##     
	if($blacklist){
	    foreach $black (@blist){
		s/ *"$black"//g ;
	    }
	}
    
    ##
    ## handle the "restrict" code here, apply each rule to decide
    ## which line to keep and possibly what part of the line
    ##
	$keep=0;
	if ($restrict) { 
	    foreach $bl (@rest){
		($regex,$before,$after) = split(/,/, $bl ) ; 
		next unless /"$regex"/ ;
		$keep=1;
        # deal with before
        
		unless ( $after =~ /^\*$/ ){ 
		    $rx = ( '"[^"]*" ' x $after).'"'.$regex.'"' ;
          s/.*($rx)/\1/ ;      
        }
        
		unless ( $before =~ /^\*$/ ){ 
		    $rx = '"'.$regex.'"'.( ' "[^"]*"' x $before) ;
          s/($rx).*/\1/ ;
		}
      
	    } 
	    next LINE if $keep == 0;
	}
    
	next LINE unless $_; 
    ## 
    ## count with what is left
    ##
	@line = reverse split(/ /) ;
	$caller = shift(@line);
	$allfun{$caller}++;
	while( $called = shift(@line) ){
	    $allfun{$called}++;
	    $calltree{$caller}{$called} ++ ;
	    $caller = $called;
	}
    
    }
}
 
##
## generate the dot code to make the graph
##
print "digraph {\n" ;
print 'graph [ rankdir = "LR"]; '."\n";
foreach $fun (keys %allfun){ 
    $_ = $fun;
  s/"$//;
    $value = $allfun{$fun} ;
    $seconds = $value * $sample ;
    $shape='rect';
  ## restrict functions are on the watchlist
    if ($restrict) {
	foreach $bl (@rest){
	    my ($regex,$before,$after) = split(/,/, $bl ) ; 
	    $regex='"'.$regex.'"';
	    $shape='ellipse,color=orange,style=filled' if $fun =~ $regex; 
	}
    }
    print "$fun [shape=$shape,fontsize=6,label=$_\\n$seconds seconds\", width=1.5] \n" if $value > ($cutoff-1) ;
}
  
foreach $caller (keys %calltree){
    for $called ( keys %{ $calltree{$caller} } ) {
	$value = $calltree{$caller}{$called} ;
	$seconds= $value * $sample ;    
	print " $caller -> $called [label=\"" . $seconds. "s\" ,fontsize=6]\n" if ( $value > $cutoff );
    }
}
print "}\n"
