use strict;
use warnings;
use Statistics::Basic qw(:all);

my %ge; &get_target_gene (\%ge);
my %fi; &fix_filename ('All_Cell_Lines',\%fi);
my $fno=0;

open (my $fh2,">all_cell_lines_result.txt");
#process each file
foreach my $tf (sort {$a cmp $b} keys %fi) {
	foreach my $f (sort {$a cmp $b} keys %{$fi{$tf}}) {
		$fno++;
		open (my $fh1,"All_Cell_Lines/$f") or die $!;
		my %rs; my %pk1; my %pk2;
		while (<$fh1>) {
			my ($ch,$ll,$rr,$dd,$ss)=split /\t/;
			next unless $ch && $ll && $rr && $ss;
			if ($ll =~ /\D/ || $rr =~ /\D/ || $ll >$rr) {
				print "wrong column 3 or 4\n"; 
				next;
			}
			$pk1{$_}=1;
			foreach my $lr (keys %{$ge{$ch}}) {
				my ($l,$r)=split /\t/,$lr;
				next if $l >$rr || $r <$ll;
				push @{$rs{$ge{$ch}{$lr}}},$ss;
				$pk2{$_}=1;
			}
		}
		close $fh1;

		print $fh2 ">$tf\t$fno\t$f\n";
		print $fh2 "all_peaks\t",get_stat(\%pk1);
		print $fh2 "overlapped_peaks\t",get_stat(\%pk2);
		my $gno=scalar keys %rs;
		print $fh2 "gene_number\t$gno\n";

		foreach my $ge (sort {$a cmp $b} keys %rs) {
			my @ss=@{$rs{$ge}};
			my $mean=mean(\@ss);
			my $sd=stddev(\@ss);
			print $fh2 "$ge\t",$#ss+1,"\t$mean\t$sd\n";
		}
	print $fh2 "\n";
	}
}
close $fh2;

#get target gene
sub get_target_gene {
	my ($ge)=@_;
	open (hand1,"promoter_list_052015.txt") or die $!;
	while (<hand1>) {
		$_ =~ s/\s+$//;
		my @a1=split (/\s+/, $_);
		if ($a1[2] =~ /\D/ || $a1[3] =~ /\D/) {
			print "start or end coordinates contains non-number\n";
		}
		my $l = $a1[2];
		my $r = $a1[3];
		$a1[1] =~ tr/a-z/A-Z/;
		${$ge}{"chr$a1[1]"}{"$l\t$r"}="$a1[0]\t$a1[1]\t$a1[2]\t$a1[3]";
	}
	close hand1;
}

#get statistics
sub get_stat {
	my ($pk)=@_;
	my $pka=scalar keys %{$pk};
	my @all;
	foreach my $p (keys %{$pk}) {
		my ($ch,$ll,$rr,$dd,$ss)=split /\t/,$p;
		push @all,$ss;
	}

	my $mean=mean(\@all);
	my $median=median(\@all);
	my $stddev=stddev(\@all);
	return "$pka\t$mean\t$median\t$stddev\n";
}

####log2
sub log2                {
my $n = shift;
return log($n)/log(2);  }    

####group file
sub fix_filename {
	my ($fna,$fa)=@_;
	opendir (dir1, $fna) or die $!;
	while (my $f=readdir(dir1)) {
		if ($f !~ /\.narrowPeak/) {
		next;
		}
		my $ff=$f;
		$ff=~ s/wgEncode//;
		$ff=~ s/Tfbs/_/;
		$ff=~ s/Gm12878/_/;
		$ff=~ s/(UniPk|StdPk)\./\./;
		$ff=~ s/(IggrabPk|IggmusPk)\./\./;
		$ff=~ s/(Iggmus|Iggrab)\./\./;
		$ff=~ s/(Pcr1x|Pcr2x|iknucla)//;
		$ff=~ s/(Znf\d)\w+/$1/;
		$ff=~ s/an{0,1}b\d+\./\./;
		$ff=~ s/(200401194|01388V0416101|sc6059|sc13268V0416101|sc81325V0422111)//;
		$ff=~ s/(12771|39875|sc101553V0422111)\./\./;
		$ff=~ s/(sc22827|sc81188V0422111|sc150V0422111|sc584|sc15914c20)\./\./;
		$ff=~ s/(a301218a|V0416101|sc631V0416101|sc502V0422111|sc137065)\./\./;
		$ff=~ s/(sc6327|sc17834V0422111|sc81335V0422111|sc71910V0422111)\./\./;
		$ff=~ s/(a300|V0416101|sc74442V0422111|sc281|sc25388V0416102|sc30189)\./\./;
		
		$ff =~ s/Hepg2|A549|Imr90|Mcf7|U2os/_/;
		
		$ff=~ s/Etoh02\./\./;
		$ff=~ s/(sc240V0416102|V0416102|sc5916)//;		
		$ff=~ s/(sc5916|sc240V0416102|sc81325V0422111|V0422111)\./\./;

		$ff=~ s/(anb100279|V0416101|a300|Forskln|sc150V0416101|sc636V0416101)\./\./;
		$ff=~ s/(ab68301|sc30189|sc5916V0416101|39875|sc101058V0416101|sc6553V0416101)\./\./;
		$ff=~ s/(sc5916|sc631|sc101058|sc6553|sc6554|sc6296|sc8987|sc6558|sc81192V0422111)\./\./;
		$ff=~ s/(200401194|ab50322|sc477|ab85725|sc271530V0422111|sc81335V0422111|sc582|V0416102)\./\./;
		$ff=~ s/(m8194|ab9263|V0422111|Ucd|UcdPk|sc281|sc101184|Pravast|Insln)\./\./;
		$ff=~ s/(sc101184)\./\./;
		
		$ff=~ s/(nb10060411|sc30189|ab85725)//;
		
		$ff=~ s/(sc269UcdUniPk|Estro|Serumstim|Serumstvd|Veh|Ucd|UcdPk)\./\./;
		$ff=~ s/(sc269)\./\./;

		$ff=~ s/(Ucd|UcdPk)\./\./;

		$ff =~ s/_/-/g;

		if ($ff =~ /-(\w+)\.narrowPeak/) {
			$ff=$1;
		}
		${$fa}{$ff}{$f}=1;
	}
	closedir dir1;
	foreach my $ff (sort {$a cmp $b} keys %{$fa}) {
		foreach my $f (sort {$a cmp $b} keys %{${$fa}{$ff}}) {
			print "$ff\t$f\n";
		}
	}
}
