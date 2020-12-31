#
$column1 = $ARGV[0];
$column2 = $ARGV[1];

if ($#ARGV+1 < 1) {
    print "\n\tUsage: copyColumn.pl <fromColumn> <toColumn> <input >output\n";
    print "\n\tNotes: The column indexing is 0-based.\n";
    exit -1;
}

foreach $line (<STDIN>) {
    chomp($line);
    @entries = split(",", $line);
    #print @entries;
    #$size = @entries;
    #print $size;

    if ($entries[$column1] ne "") {
	$entries[$column2] = $entries[$column1];
	print(join(",",@entries), "\n");
    }
    else {
	print(join(",",@entries), "\n");
    }

}


