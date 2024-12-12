#!/usr/bin/perl -w

$offset = 10;

while($inline=<>){
    @lines = split(/[\s ]+/,$inline);
    $start = $lines[9];
    $end = $lines[10];
    if($lines[9] > $lines[10]){$start = $lines[10]; $end= $lines[9];}

    @nuc = split(/_/,$lines[2]);
    if($lines[2] =~ /_Complete/){$start2 = 201; $end2 = 201 + $nuc[1] - $nuc[0] + 1;}
    if($lines[2] =~ /_Ins/ || $lines[2] =~ /_Del/){$start2 = 201; $end2 = 201;}

    $ctl = 0;
    if($start < $start2 -$offset && $end >= $start2 + $offset){$ctl = 1;}
    if($start < $end2 -$offset && $end >= $end2 + $offset){$ctl = 1;;}

    if($ctl == 1){print $inline;}
}

