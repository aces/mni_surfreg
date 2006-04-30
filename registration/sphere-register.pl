#! xPERLx -w

use strict;
use MNI::Startup qw(optvars opttable);
use Getopt::Tabular;
use MNI::Spawn;
use MNI::DataDir;
use MNI::FileUtilities;

my $show_version = 0;

my @opt_table =
  (@DefaultArgs,
   [ "-version", "const", 1, \$show_version, "display version information and exit" ],
  );
GetOptions( \@opt_table, \@ARGV ) || exit 1;

if ( $show_version ) {
  print "sphere-register from mni_surfreg version xVERSIONx\n";
  $Verbose = 0; # disable printing elapsed CPU time at exit
  exit 0;
}

die "usage: $0 [options] source target smap_out\n" unless @ARGV==3;
my( $source, $target, $map_final ) = @ARGV;

die "Output map [$map_final] exists; use -clobber to overwrite\n"
  if ( -f $map_final && ! $Clobber );

# The programs used.
# Search for them first in xBINDIRx.
RegisterPrograms( [qw(create_tetra
		      initial-surface-map
		      surftracc
		      refine-surface-map
		      cp
		      wc)],
		  "xBINDIRx:" . $ENV{'PATH'} );

AddDefaultArgs( 'surftracc', ['-debug',1] ) if $Debug;

# Create temporary work directory, and generate a set of spherical
# model files.
#
$TmpDir .= "/" unless ( $TmpDir =~ m,.*/$, );
MNI::FileUtilities::check_output_path($TmpDir);

my $max_num_faces = get_num_faces( $source );
my $target_num_faces = get_num_faces( $target );

die "Source has $max_num_faces but target has $target_num_faces\n"
  unless ( $max_num_faces == $target_num_faces );

print "Input surfaces have $max_num_faces faces\n"
  if $Verbose;

my $num_faces = 1280;
my $map = generate_initial_map( $num_faces, $max_num_faces );
$map = register( $num_faces, $max_num_faces, $map );

while ( $num_faces < $max_num_faces ) {
    $num_faces = $num_faces * 4;

    $map = refine( $num_faces, $max_num_faces, $map );
    $map = register( $num_faces, $max_num_faces, $map );
}

Spawn(['cp', '-f', $map,$map_final]);


sub get_num_faces {
    my $vv_file = shift;
    my $out;
    Spawn([ 'wc', '-l', $vv_file ], stdout => \$out );

    return 81920 unless $Execute;


    my @out = split( ' ', $out );
    die "Unexpected output: $out\n" unless @out == 2;

    my %vertex_to_face_map = (
          642 => 1280,
         2562 => 5120,
        10242 => 20480,
        40962 => 81920,
       163842 => 327680
    );

    my $n_vert = $out[0];
    my $n_face = $vertex_to_face_map{$n_vert};

    die "Number of vertices [$n_vert] not handled\n"
      unless defined $n_face;

    return $n_face;
}

sub get_model_file {
    my $f = shift;

    my $filename = $TmpDir . "model_$f.obj";
    if ( ! -f $filename ) {
      Spawn([ 'create_tetra', $filename, 0,0,0, 1,1,1, $f ]);
    }
    return $filename;
}

sub generate_initial_map {
    my( $f1, $f2 ) = @_;
    my $model1 = get_model_file( $f1 );
    my $model2 = get_model_file( $f2 );
    my $fn = $TmpDir . "initial.sm";
    Spawn([ 'initial-surface-map', $model1, $model2, $fn ]);
    return $fn;
}

sub register {
    my( $f1, $f2, $map_in ) = @_;
    my $model1 = get_model_file( $f1 );
    my $model2 = get_model_file( $f2 );
    my $map_out = $TmpDir . "map_reg_" . $f1 . ".sm";
    Spawn([ 'surftracc',
	    $model2, $source,
	    $model2, $target,
	    $model1,
	    $map_in, $map_out ]);
    return $map_out;
}

sub refine {
    my( $f1, $f2, $map_in ) = @_;
    my $model1 = get_model_file( $f1 );
    my $model2 = get_model_file( $f2 );
    my $map_out = $TmpDir . "map_ref_" . $f1 . ".sm";
    Spawn([ 'refine-surface-map',
	    $map_in, $model1, $model2, $map_out ]);
    return $map_out;
}
