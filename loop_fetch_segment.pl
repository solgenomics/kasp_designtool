
#boucle_fetch_segment.pl

#Search ID in a multifasta  and fetch sequence fragment by launching python fetch_segment
# outputs are in individual files

#!/usr/bin/perl -w
#syntax:  perl boucle_fetch_segment.pl input_file

#Open file
#use strict;
#use warnings;

my $fichlu=$ARGV[0];

print "fichlu value is $ARGV[0] and $fichlu\n";

open(IDENT, "$fichlu") || die ("File doesn't exist");

# Corps du script

foreach my $line (<IDENT>) {

    print "line = $line\n";

    chop($line);

    #if ( $line eq "") {continue}
    #if (/^(\s)*$/) { print "ligne vide\n"; continue;}

    my @tab=split('\t', $line);

#Attribute variables

    my $chromo = $tab[0];
    my $gene = $tab[1];
    my $pi = $tab[2];
    my $pf = $tab[3];

    #print("chromo = $chromo \n");
    #print("gene = $gene \n");
    #print("pi = $pi \n");
    #print("pf = $pf \n");

# add var to launch python fetch_segment script
#loop on cassava genome

       my $command= " /home/gjb99/GBS/tools/KaSpaR/fetch_segment -id " . $chromo . " -pi " . $pi . " -pf " . $pf . " -db /home/gjb99/GBS/manihot_ref6/Mesculenta/Manihot_v6.1/assembly/Mesculenta_305_v6_18K/Mesculenta_305_v6_18K.fa -new_header ".$gene." -out " . $gene . "\.fasta"."\n";

#/Users/guillaume/Desktop/manihot_ref6/Mesculenta/Manihot_v6.1/assembly/Mesculenta_305_v6.fa
		print $command. "\n";

 		open(FASTA,">gene.fasta");
		print FASTA $gene."\n".$command;
		close(FASTA);

#launch in bash command

   		system($commande) == 0 or die "la commande a échoué\n";

}


close(IDENT);
