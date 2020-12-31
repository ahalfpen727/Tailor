#!/usr/bin/envd perl

$commentChar = $ARGV[0];
# print "commentChar = $commentChar\n";

if ($#ARGV+1 < 1) {
    print "\n\tUsage: stripComments.pl char <input >output\n";
    print "\n\tNotes: Output contains all text to the left of the first comment char.\n";
    exit -1;
}

sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };
# sub trim {
#    return $_[0] =~ s/^\s+|\s+$//rg;
# }

foreach $line (<STDIN>) {
  chomp($line);
  $line=trim($line);
  # print "line = $line\n";

  $firstChar = substr($line, 0, 1);
  $sevenChars = substr($line, 0, 7);

  if ($firstChar eq $commentChar) {
    # print "firstChar = $firstChar\n";
    # do nothing!
  }
  elsif ($sevenChars eq "export ") {
    # do nothing!
  }
  elsif ($line =~ /\(/){
    # do nothing!
  }
  elsif ($line =~ /\[/){
    # do nothing!
  }
  elsif ($line !~ /=/){
    # do nothing!
  }
  else {
    @entries = split($commentChar, $line);
    $length=@entries;
    # print "length = $length\n";

    print "$entries[0]\n";

    #print @entries;
    #$size = @entries;
    #print $size;
  }
}


