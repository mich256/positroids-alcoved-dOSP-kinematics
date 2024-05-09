package cyclic_chain_polytopes;
use parent 'Exporter';
our @EXPORT_OK = qw(chains_to_polytope1 chains_to_polytope2);

sub chains_to_polytope1 {
	my @chains = @_;
	my @arr = map { map {[0]} 0..$n } 0..$#chains+2*$n;
	for (my $i = 0; $i < $n; $i++) {
		$arr[$i][$i+1] = 1;
		$arr[$i+$n][0] = 1;
		$arr[$i+$n][$i+1] = -1;
	}
	my $count = 2*$n;
	foreach $c (@chains) {
		foreach (@$c) {
			$arr[$count][0] = 1;
			$arr[$count][$_] = -1;
		}
		$count++;
	}
	return @arr;
}

sub right_circ_des {
	my @Y = @_;
	my $size = @Y;
	my $a = $Y[0];
	my $b = 0;
	my $cdesr = 0;
	for (my $i = 1; $i < $size; $i++) {
		if ($i == $size-1 & @Y[$i] > @Y[0]) {
			$cdesr++;
			last;
		}
		$b = $Y[$i];
		if ($a > $b) {
			$cdesr++;
		}
		$a = $b;
	}
	return $cdesr;
}

sub chains_to_polytope2 {
	my @chains = @_;
	my @arr = map { map {[0]} 0..$n } 0..$#chains+2*$n;
	for (my $i = 0; $i < $n; $i++) {
		$arr[$i][$i+1] = 1;
		$arr[$i+$n][0] = 1;
		$arr[$i+$n][$i+1] = -1;
	}
	my $count = 2*$n;
	foreach $c (@chains) {
		$arr[$count][0] = right_circ_des(@$c);
		foreach (@$c) {
			$arr[$count][$_] = -1;
		}
		$count++;
	}
	return @arr;
}

	
