#! /usr/bin/perl

# shrink elements in an existing OpenDX mesh file to make it easier
# to visualize them

# David Michea, University of Pau, France, January 2007

use Data::Dumper;

if (scalar(@ARGV) < 2)
{     print "usage : $0 file2split.dx filesplit.dx [spacing ratio (try 1)]\n";
      exit;
}
$filein = $ARGV[0];
$fileout = $ARGV[1];
$distance = $ARGV[2] || 0.3;
$numnode=0;
$nb_corners=8;

# read datas
open (DF, "$filein");
while(my $line=<DF>)
{     chomp $line;
      if ($line =~ m/^(\s+\d+\.?\d*\s*)+$/)
      {     if ($line =~ m/\./)
            {     my @temp=split(/\s+/, $line);
                  shift @temp;

                  $nodes[$numnode] = \@temp;
                  $numnode++;
            }
            elsif ($line =~ m/(\s+\d+\s*){8}/)
            {     my @temp=split(/\s+/,$line);
                  shift @temp;
                  $hex[$numhex] = \@temp;
                  $numhex++;
            }
            else
            {     # couleur
            }
      }
}
close (DF);
# dump datas into structure
my @cubes;
foreach $cube (@hex)
{     my @corners;
      foreach $corner (@{$cube})
      {     my @coords = @{$nodes[$corner]};
            push(@corners, \@coords);
      }
      push (@cubes, \@corners);
}
# split the brick
eclate($distance,@cubes);
# write dx file
dump_dx(@cubes);
exit;

sub eclate
{     my ($dist, @cubes) = @_;
      my ($mid_x,$mid_y,$mid_z) = mid_pos(@cubes);
      foreach $cube(@cubes)
      {     ($mx, $my, $mz) = moy_pos($cube);
            $dx = ($mx-$mid_x)*$dist;
            $dy = ($my-$mid_y)*$dist;
            $dz = ($mz-$mid_z)*$dist;
            redefine($cube,$dx,$dy,$dz);
      }
}

sub mid_pos
{     my (@cubes) = @_;
      my ($mx,$my,$mz)=(1000,1000,1000);
      my ($Mx,$My,$Mz)=(-1000,-1000,-1000);
      my $compt=0;
      foreach $cube(@cubes)
      {     foreach my $corner (@{$cube})
            {     my ($tx,$ty,$tz) = @$corner;
                  $mx=$tx if ($tx < $mx);
                  $my=$ty if ($ty < $my);
                  $mz=$tz if ($tz < $mz);
                  $Mx=$tx if ($tx > $Mx);
                  $My=$ty if ($ty > $My);
                  $Mz=$tz if ($tz > $Mz);
            }
      }
      return (($Mx+$mx)/2,($My+$my)/2,($Mz+$mz)/2);
}

sub moy_pos
{     my $cube = shift;
      my ($x,$y,$z)=(0,0,0);
      my $compt=0;
      foreach my $corner (@{$cube})
      {     my ($tx,$ty,$tz) = @$corner;
            $x+=$tx;
            $y+=$ty;
            $z+=$tz;
            $compt++;
      }
      ($x,$y,$z) = map ($_/$compt, ($x,$y,$z));
      return ($x,$y,$z)
}

sub redefine
{     my ($cube,$dx,$dy,$dz) = @_;
      foreach my $corner (@{$cube})
      {     $corner->[0]+=$dx;
            $corner->[1]+=$dy;
            $corner->[2]+=$dz;
      }
}

sub dump_dx
{     my (@cubes) = @_;
      open (OUT, ">$fileout");
      $nb_nodes = $numhex * $nb_corners;
      print OUT " object 1 class array type float rank 1 shape 3 items         $nb_nodes  data follows\n";
      foreach my $uncube (@cubes)
      {     foreach my $uncoin (@$uncube)
            {     print OUT join(" ", @$uncoin),"\n";
            }
      }
      print OUT " object 2 class array type int rank 1 shape 8 items         $numhex  data follows\n";
      for ($i=0;$i<$nb_nodes;$i++)
      {     print OUT $i," ";
            print OUT "\n" unless (($i+1)%$nb_corners);
      }
      print OUT "\n\n",' attribute "element type" string "cubes"',"\n",'attribute "ref" string "positions"',"\n",'object 3 class array type float rank 0 items         ',$numhex,"  data follows\n";
      for ($i=1;$i<$numhex+1;$i++) {print OUT $i,"\n";}
      print OUT <<EOF
 attribute "dep" string "connections"
 object "irregular positions irregular connections" class field
 component "positions" value 1
 component "connections" value 2
 component "data" value 3
 end
EOF
;
      close (OUT);
}

