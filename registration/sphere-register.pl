#! xPERLx -w

use strict;
use MNI::Startup;
use Getopt::Tabular;
use MNI::Spawn;
use MNI::DataDir;
use MNI::FileUtilities qw(test_file check_output_dirs);

die "usage: $0 source target smap_out" unless @ARGV==3;
my( $source, $target, $map_final ) = @ARGV;

# The programs used.  
# Search for them first in xBINDIRx.
RegisterPrograms( [qw(create_tetra
		      initial-surface-map
		      surftracc
		      refine-surface-map
		      cp)],
		  "xBINDIRx:" . $ENV{'PATH'} );

AddDefaultArgs( 'surftracc', ['-debug',1] );

# Create temporary work directory, and generate a set of spherical
# model files.
# 
check_output_dirs($TmpDir);

my %model;
foreach my $f (1280, 5120, 20480, 81920) {
    $model{$f} = $TmpDir . "model_$f.obj";
    Spawn(['create_tetra',$model{$f},0,0,0,1,1,1,$f]);
}


my $map = generate_initial_map(1280);
$map = register(1280,$map);

$map = refine(5120,$map);
$map = register(5120,$map);

$map = refine(20480,$map);
$map = register(20480,$map);

$map = refine(81920,$map);
$map = register(81920,$map);

Spawn(['cp',$map,$map_final]);


sub generate_initial_map {
    my $f = shift;
    my $fn = $TmpDir . "initial.sm";
    Spawn(['initial-surface-map', $model{$f}, $model{81920}, $fn]);
    return $fn;
}

sub register {
    my( $f, $map_in ) = @_;
    my $map_out = $TmpDir . "map_reg_" . $f . ".sm";
    Spawn(['surftracc', 
	   $model{81920},$source, 
	   $model{81920},$target,
	   $model{$f},
	   $map_in, $map_out]);
    return $map_out;
}

sub refine {
    my( $f, $map_in ) = @_;
    my $map_out = $TmpDir . "map_ref_" . $f . ".sm";
    Spawn(['refine-surface-map', 
	   $map_in, $model{$f}, $model{81920}, $map_out]);
    return $map_out;
}

