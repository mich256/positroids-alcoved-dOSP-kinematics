# circular shift an array. Thanks stack overflow https://stackoverflow.com/questions/29491286/circular-shift-of-array

sub circ_shift {
  my ($list, $n) = @_;
  push @$list, splice @$list, 0, $n;
}

sub chain_to_array {
  my @chain = @_;
  my $size = @chain; # last index of chain, lengh of chain -1
  my @arr = map { map { [0] } 0..$n } 1..$size; # create a zero array with n+1 columns and size rows
  my $count = 1;
  $arr[0][$chain[0]] = 1;
  for (my $i = 1; $i < $size; $i++) {
    $arr[$count][$chain[$i-1]] = -1;
    $arr[$count][$chain[$i]] = 1;
    $count++;
  }
  $arr[$count][$chain[$size]] = -1;
  $arr[$count][0] = 1;
  return @arr;
}

sub chain_to_fan {
  my @chain = @_;
  my $size = @chain; # last index of chain, lengh of chain -1
  my @arr_cones = map { [] } 0..$size; # create an empty array
  for (my $i=0; $i < $size; $i++) {
    my @arr = chain_to_array(@chain);
    my $m = new Matrix<Int>(\@arr);
    my $c = new Cone(INEQUALITIES=>$m);
    $arr_cones[$i] = $c;
    circ_shift(@chain);
  }
  my $a = new Array<Int>(\@arr_cones);
  return check_fan_objects($a); 
}

sub intersect_arr_fans {
  my @arr_fans = @_;
  my $cc = $arr_fans[0];
  foreach (@arr_fans) {
    $cc = common_refinement($cc, $_);
  }
  return $cc;
}

sub cyclic_order_complex {
  my @max_chains = @_;
	my $size = scalar @max_chains;
  my @arr_fans = map { [] } 1..$size;
  my $count = 0;
  foreach (@max_chains) {
    $arr_fans[$count] = chain_to_fan($_);
    $count++;
  }
  return intersect_arr_fans(@arr_fans);
}
