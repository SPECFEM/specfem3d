package FILE_TOOLS;
# Qinya Liu, Caltech, May 2007

use warnings;
use Exporter;
use lib '/opt/seismo-util/lib/perl';

@ISA = ('Exporter');
@EXPORT = ('list2files','files2list','file_in_files','files_minus',
          );

# list2files

sub list2files {
  my($list) = @_;
  if (@_ != 1) {die("list2files list\n");}
  my(@outfile) = ();
  open(LIST,"$list") or die("Error open file $list\n");
  @list = <LIST>;
  foreach $list (@list) {
    ($list) = split(" ",$list);
    if ($list !~/^\s*$/) {push(@outfile,$list);}
  }
  close(LIST);
  return(@outfile);
}

# files2list

sub files2list{
  my($list,@files) = @_;
  if (@_ < 2) {die("files2list list files\n");}
  open(LIST,">$list") or die("Error opening file $list for writing\n");
  foreach $file (@files) {
    print LIST "$file\n";
  }
  close(LIST);
}

#file_in_files
sub file_in_files {
  my($onefile,@files) = @_;
  my($file);
  foreach $file (@files) {
    if ($file eq $onefile) {return 1;}
  }
  return 0;
}


# files_minus
sub files_minus {
  my($files,$mfiles) = @_;
  my(@files,@mfiles,@outfiles,$file);
  (@files) = split(" ",$files); (@mfiles) = split(" ",$mfiles);
  foreach $file (@files) {
	if (not file_in_files($file,@mfiles)) {push(@outfiles,$file);}
  }
  return(@outfiles);
}

# ls2files -- most of its functionality can be replaced by 
#             glob("pattern"), but pattern
#             needs to include at least * or ?.
sub ls2files {
  if (@_ != 1) {die("usage: ls2files(pattern)\n");}
  my($pattern) = @_;
  my(@files,$file,@outfiles);
  @files = `ls -1d $pattern 2> /dev/null`;
  foreach $file (@files) {
    chomp($file); if (-f $file) {push(@outfiles,$file);}}
  return(@outfiles);
}


1;
