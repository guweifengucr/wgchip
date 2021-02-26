use strict;
use warnings;
use Parallel::ForkManager;
use Statistics::Basic qw(:all);

my %ge; &get_target_gene (\%ge); my %rs; my %av; my $id;
open (my $fh,"all_cell_lines_result.txt") or die $!;
while (<$fh>) {
	$_ =~ s/\s+$//;
	next if !/\w/;
	my @a1=split (/\t/, $_);
	if (/^>/) {
		$a1[0] =~ s/^>//;
		$id=$a1[0]."|$a1[1]";
		foreach my $ge (keys %ge) {
			$rs{$ge}{$id}=0;
		}
		next;
	}
	if (/^all_peaks/) {
		$a1[2] =int($a1[2]+0.5);
		$a1[4] =int($a1[4]+0.5);
		$av{$id}=$a1[2].'|'.$a1[4];
		next;
	}
	next if /^(overlapped_peaks|gene_number)/;
	$rs{$a1[0]}{$id}=$a1[5];
}
close $fh;

open ($fh,">all_cell_lines_result_reformated.txt");
my $na; my $va; my $lb; my $sd;
foreach (sort {$a cmp $b} keys %av) {
	my ($n,$l)=split /\|/;
	my ($v,$s)=split /\|/,$av{$_};
	$na.="\t$n";
	$lb.="\t$l";
	$va.="\t$v";
	$sd.="\t$s";
}
print $fh "label$lb\n";
print $fh "TF$na\n";
print $fh "avg_PK$va\n";
print $fh "stddev$sd\n";
foreach my $ge (sort {$a cmp $b} keys %rs) {
	print $fh $ge;
	foreach my $tf (sort {$a cmp $b} keys %{$rs{$ge}}) {
		my $va=$rs{$ge}{$tf}; $va =~ s/\,//;
		$va=int($va+0.5);
		print $fh "\t$va";
	}
	print $fh "\n";
}
close $fh;

#get target gene
sub get_target_gene {
	my ($ge)=@_;
	open (hand1,"promoter_list_052015.txt") or die $!;
	while (<hand1>) {
		$_ =~ s/\s+$//;
		my @a1=split /\s+/;
		${$ge}{$a1[0]}=1;
	}
	close hand1;
}
