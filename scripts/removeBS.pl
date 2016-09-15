# !/usr/bin/perl
# Takes a tree in parenthesis format and changes it so that bootstrap values are removed

#print("$ARGV[0]\n");
$intree=$ARGV[0];

$outtree=$intree."_noBS";
my %h;
$newtree="";
$inseq=0;
$inBS=0;
$inlength=0;
$start=1;
open ( INTREE,"<$intree")|| die "Erreur de lecture $intree, Erreur: $!";
open ( OUTTREE,">$outtree")|| die "Erreur de lecture $outtree, Erreur: $!";

while (<INTREE>) {
  chomp;
  if (/\(/) {
    $temp=$_;
   
    @list=split(//,$temp); #splits the tree in characters
    foreach $ch (@list) {
   
      if (($ch eq ' ')) #a GC content is associated to each sequence name (after a space character), or a length if there is no bootstrap value
	{
	  #print("seqname $seqname");
	  $newtree=$newtree.$seqname;
	  $seqname="";
	  $inseq=0;
	  $inlength=0;
	  $inBS=1;
	  $start=0
	}
      elsif(($ch eq ';')){
	$newtree=$newtree.$ch."\n";
	#$newtree="";
        $inseq=0;
        $inBS=0;
        $inlength=0;
	}
      elsif(($ch eq ')')){
	$newtree=$newtree.$ch;
	$inBS=1;
      }
      elsif($ch eq ':') {
	#print("seqname $seqname");
	$newtree=$newtree.$seqname;
	$inBS=0;
	$seqname="";
	$inseq=0;
	$newtree=$newtree.":";
	$inlength=1;
      }
      elsif (($ch=~/[0-9]/)&&($inlength!=1)&&($inseq!=1)&&($start!=1))
	{
	  $inBS=1;
	}
      elsif ((!($ch=~/\w/))&&($inseq!=1)&&($inBS!=1))
	{
	  $newtree=$newtree.$ch;
	 # $inlength=0;
	}
      elsif ($ch=~/[A-Z]/)
	{
	  $seqname=$seqname.$ch;
	  $inseq=1;
	  $inlength=0;
	}
      elsif ($inseq==1)
	{
	  $seqname=$seqname.$ch;
	  $inlength=0;
	}
      elsif ($inBS==1)
	{  $inlength=0;}
	#elsif ($ch eq ';') {
	#print OUTTREE "$newtree;\n";
	#$newtree="";	
	#}
      else {
	$newtree=$newtree.$ch;
	   }	
    }
  }
}
if ($newtree =~/\(/ ) {
  print OUTTREE "$newtree";
}

close INTREE;
close OUTTREE;
