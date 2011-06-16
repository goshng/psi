use strict;
use Parallel::ForkManager;
my $max_procs = 50;
#my $pdbsws = "/home/goshng/work/pdbsws"; 
my $output_dir = "/home/goshng/work/pdm/psi/trunk/mosix/output";
my $data_dir = "/home/goshng/work/pdbsws";
my $pdbsws = $data_dir;
my @names = qw( Fred Jim Lily Steve Jessica Bob Dave Christine Rico Sara );

# Procedure prototypes
sub run_psi ($);
sub notify_psi ($$);

# hash to resolve PID's back to child specific information
my $pm =  new Parallel::ForkManager($max_procs);

# Setup a callback for when a child finishes up so we can
# get it's exit code
$pm->run_on_finish(
  sub { my ($pid, $exit_code, $ident) = @_;
    #print "** $ident just got out of the pool ".
    #  "with PID $pid and exit code: $exit_code\n";
    if ($exit_code != 0)
      {
        notify_psi ($ident, "Fail $exit_code");
      }
  }
);
$pm->run_on_start(
  sub { my ($pid,$ident)=@_;
    #print "** $ident started, pid: $pid\n";
    notify_psi ($ident, "Start");
  }
);
$pm->run_on_wait(
  sub {
    #print "** Have to wait for one children ...\n"
  },
  300
  #3
);

sub notify_psi ($$)
{
  my $pdbid = shift;
  my $msg = shift;

  my $begin_file = $pdbid."mail";
  open BFILE, ">$begin_file";
  print BFILE "PDB FILE: $pdbid";
  close BFILE;
  my $command_line = "mail -s \"Mosixpdbsws$msg\" goshng\@gmail.com < $begin_file";
  system $command_line;
  unlink $begin_file;
}

sub run_psi ($)
{
  my $file = shift;
  my $r = 0;

  my $pdb_id = substr $file, 3, 4;
  my $chain_id = substr $file, 8, 1;
  #print "start to run psi\n";
  my $command_str = "GSL_RNG_SEED=`date +\%s` nice mosrun -j2-13 ../../bd/src/psi".
                    " --exe-all".
                    " --auto-data $pdb_id$chain_id".
                    " --directory $output_dir".
                    " --data-directory $data_dir".
                    " --burn-in 1".
                    " --sample-freq 2".
                    " --sample-size 1000".
                    " --gridpoint-p-begin -0.5".
                    " --gridpoint-p-end 0.5".
                    " --gridpoint-p-number 101".
                    " --gibbs-size 10".
                    " --gibbs-burn 2".
                    " --gibbs-freq 2".
                    " --limit-tau 1.1".
                    " --limit-psr 1.01\n";
  $r = system $command_str;
  #$r = system "cp 1";

  if ($? == -1) {
    #print "failed to execute: $!\n";
    $r = -1;
  }
  elsif ($? & 127) {
    #printf "child died with signal %d, %s coredump\n",
    #    ($? & 127),  ($? & 128) ? 'with' : 'without';
    $r = 127;
  }
  else {
    #printf "child exited with value %d\n", $? >> 8;
    $r = $? >> 8
  }

  #print $command_str;
  #$r = 1;
  return $r;
}

# Main functions

opendir PDBSWS, $pdbsws;

while (my $file = readdir (PDBSWS))
#foreach (sort readdir <PDBSWS>)
{
  next unless $file =~ /prd$/;
  my $pid = $pm->start($file) and next;
  # This code is the child process
  my $r = run_psi ($file); 
  $pm->finish($r); # pass an exit code to finish
}
closedir PDBSWS;

#print "Waiting for Children...\n";
$pm->wait_all_children;
#print "Everybody is out of the pool!\n";

