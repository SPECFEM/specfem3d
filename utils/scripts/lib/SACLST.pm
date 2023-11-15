package SACLST;

# SACLST.pl
use strict;
use warnings;
use TrimWhitespace;

BEGIN {
    use Exporter();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    $VERSION     = 1.00;

    @ISA         = qw(Exporter);
    @EXPORT      = qw(&sac_files_get_values &sac_files_prune &sac_files_sort &sac_files_match &sac_files_skip);
    %EXPORT_TAGS = ( );     # eg: TAG => [ qw!name1 name2! ],
    @EXPORT_OK   = qw( );
}

our @EXPORT_OK;

my $saclst      = "/opt/seismo-util/bin/saclst";

# Subroutine: sac_files_get_values
#      Agruments: keyword - sac heade keyword
#                 files   - a list of sacfiles
#      Purpose: to retun a list of values from a list of sacfiles
#               Multiple values can be obtained at the same time by
#               using multiple keywords in a single string and they
#               are then returned as a string in the same order
sub sac_files_get_values {
    my($keyword, @files) = @_;
    my (@values, $line, $name, $value);
    if($keyword =~ m/^\s*$/) {
	return();
    }
    open(SACLST, "$saclst $keyword f @files|") or die "Error opening saclst";
    while($line = <SACLST>) {
	chomp($line);
	($name, $value) = split(/\s+/, $line, 2);
	$value =~ s/\s+/ /g;
	$value = TrimWhitespace($value);
	push(@values, $value);
    }
    close(SACLST);
    return(@values);
}


# Subroutine: sac_files_prune
#      Argument: keyword     - a sac header keyword to prune on
#                small_value - the smallest value desired
#                large_value - the largest value desired
#                files       - a list of sacfiles
#      Purpose: to remove sacfiles outside a range of values
#      Bugs: Only works on numerical values
sub sac_files_prune {
    my($keyword, $small_value, $large_value, @files) = @_;
    my (@files2, $name, $value, $line);
    open(SACLST, "$saclst $keyword f @files|") or die "Error opening saclst";
    while($line = <SACLST>) {
	($name, $value) = split(/\s+/, $line);
	if (not defined $value) {die("Check saclst : value of $keyword\n");}
	if($value > $small_value && $value < $large_value) {
	    push(@files2, $name);
	}
    }
    close(SACLST);
    return(@files2);
}


# Subroutine sac_files_sort
#      Arguments: keyword - a sac header keyword to sort on
#                 updown  - direction; < 0 is reverse direction
#                 files   - list of sacfiles
#      Purpose: to sort a list of sacfile depending on a particular value
#      Bugs: Onlt works on numerical values
sub sac_files_sort {
    my($keyword, $updown, @files) = @_;
    my (@files2, %values, $line, $name, $value);
    open(SACLST, "$saclst $keyword f @files|") or die "Error opening saclst";
    while($line = <SACLST>) {
	($name, $value) = split(/\s+/, $line);
	$values{$name} = $value;
    }
    close(SACLST);
    @files2 = sort { $values{$a} <=> $values{$b} } keys %values;
    if($updown < 0) {
	@files2 = reverse(@files2);
    }
    return(@files2);
}

# subroutine sac_files_match
#      Arguments : keywords - a sac header keyword to match
#                  values - value of the header to match
#                  files - sac files to match
#      Purpose : find in the sac files those that match keyword has
#                the desired value
#                maybe helpful in taking certain component, station,
#                network etc. from a large set of sacfiles
#      Bugs:   Only works on k headers

sub sac_files_match {
  my($keywords,$values,@files) = @_;
  my(@keywords) = split(" ",$keywords);
  my(@values) = split(" ",$values);
  if ($#keywords != $#values ) {die("The number of keywords does not match that
  of the values\n");}
  if (@keywords == 0) {die("No keywords get passed to sub\n");}
  my(@outfiles) = ();
  my(@file_values) = ();
  my($i,$file);

  foreach $file (@files) {
    (undef,@file_values) = split(" ",`$saclst @keywords f $file`);
    for ($i = 0;$i<@values; $i++) {
      if ($values[$i] ne $file_values[$i]) {last;}
    }
    if ($i == @values) {push(@outfiles,$file);}
  }
  return(@outfiles);

}

# Subroutine sac_files_skip
#      Arguments: n in    - number (integer) ...
#                 files - list of sacfiles
#      Purpose: takes a list of sac file names and 
#      Returns a list with only every $n-th item
#      Basically List Reduction
sub sac_files_skip {
    my($nn, @files) = @_;
    my($i, $n,$in,@ret_files);
    ($n,$in) = split(" ",$nn);
    if (not $in)  {$in = 0;}
    if ($in >= $n or $in < 0) {die("only in = 0,n-1\n");}
    if($n <= 1) {
	return(@files);
    }
    for($i = $in-1; $i < @files; $i += $n) {
	push(@ret_files, $files[$i]);
    }
    return(@ret_files);
}

END { }

1;
