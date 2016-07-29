#!/usr/bin/bash

#bash make_me.sh

#get points
export SCRIPT_PATH=../../scripts

perl -e 'my @upper = (1,1,1);
my @lower = (0,0,0);
for(my $i = 0; $i < 30000; $i++)
{
  my $x = rand() * ($upper[0] - $lower[0]) + $lower[0];
  my $y = rand() * ($upper[1] - $lower[1]) + $lower[1];
  my $z = rand() * ($upper[2] - $lower[2]) + $lower[2];
  print "$x $y $z\n";
}' > pts
perl $SCRIPT_PATH/points2minmax.pl pts > mm

bash run_me.sh pts 0 0.1

rm pts
