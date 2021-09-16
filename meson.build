project('gevolution','cpp','c',
        default_options : ['cpp_std=c++17',
        'warning_level=3',
        'optimization=3'],
        version: '1.2')

cpp=meson.get_compiler('cpp')
mpi=dependency('mpi',language: 'cpp')
boost = dependency('boost', modules: ['mpi','serialization','unit_test_framework'])
fftw3 = dependency('fftw3')
hdf5 = dependency('hdf5')
gsl = dependency('gsl')
latfield = dependency('latfield')
mpirun = find_program('mpirun')
bash = find_program('bash')
pkg=import('pkgconfig')

conf_data = configuration_data({ 
    'PHINONLINEAR' :get_option('PHINONLINEAR'),
    'BENCHMARK' :get_option('BENCHMARK'),
    'EXACT_OUTPUT_REDSHIFTS' :get_option('EXACT_OUTPUT_REDSHIFTS'),
    'VELOCITY' :get_option('VELOCITY'),
    'COLORTERMINAL' :get_option('COLORTERMINAL'),
    'CHECK_B' :get_option('CHECK_B'),
    'HAVE_CLASS' :get_option('HAVE_CLASS'),
    'HAVE_HEALPIX' :get_option('HAVE_HEALPIX'),
    'GEVOLUTION_OLD_VERSION' :get_option('GEVOLUTION_OLD_VERSION')
    })

deps = [mpi,fftw3,hdf5,gsl,boost,latfield]

subdir('include')
subdir('src')


libgevolution = library('gevolution',
    gevolution_headers,
    gevolution_sources,
    config,
    include_directories: include,
    install: true,
    dependencies: deps)

install_headers(gevolution_headers,subdir: 'gevolution')
install_headers(config, subdir: 'gevolution')

pkg.generate(libgevolution,
    description: 'gevolution')

gitinfo_generator = files('git_version.sh')
version_template = files('version.template')

version = custom_target('git-info',
    build_by_default: true,
    input : version_template,
    output: 'version.h',
    command : [bash,gitinfo_generator,meson.source_root(),'@OUTDIR@','@INPUT@'],
    build_always_stale: true)

gevolution = executable('gevolution',
    main,version,
    link_with: libgevolution,
    dependencies: deps,
    include_directories: include)
