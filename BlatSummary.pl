#!/usr/bin/perl -w

%size = ();
$lib = `cat library_info.txt`;
@libs = split(/\n/,$lib);
for($i = 0; $i < @libs; $i++){
    @lines = split(/[\s ]+/,$libs[$i]);
    $size{$lines[0]} = $lines[1];
}

$lencutoff = 0.8;

%count = ();

while($inline=<>){
    @lines = split(/[\s ]+/,$inline);
    if($lines[3] >= 90 && $lines[4] >= $size{$lines[0]} * $lencutoff){
	$handle = $lines[0].";".$lines[2];
	if(!defined($count{$handle})){$count{$handle} = 0;} 
	$count{$handle} = $count{$handle} +1;
    }
}

@id = keys %count;

for($i=0; $i < @id; $i++){
    $handle = $id[$i];
    $handle =~ s/;/\t/;
print $handle, "\t", $count{$id[$i]}, "\n";}


