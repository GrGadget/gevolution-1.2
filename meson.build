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

include=include_directories('.')
conf_data = configuration_data({ 
    'PHINONLINEAR' :get_option('PHINONLINEAR'),
    'BENCHMARK' :get_option('BENCHMARK'),
    'EXACT_OUTPUT_REDSHIFTS' :get_option('EXACT_OUTPUT_REDSHIFTS'),
    'VELOCITY' :get_option('VELOCITY'),
    'COLORTERMINAL' :get_option('COLORTERMINAL'),
    'CHECK_B' :get_option('CHECK_B'),
    'HAVE_CLASS' :get_option('HAVE_CLASS'),
    'HAVE_HEALPIX' :get_option('HAVE_HEALPIX')
    })

deps = [mpi,fftw3,hdf5,gsl,boost,latfield]

gevolution_headers = files([
    'background.hpp',
    'class_tools.hpp',
    'debugger.hpp',
    'gevolution.hpp',
    'hibernation.hpp',
    'ic_basic.hpp',
    'ic_prevolution.hpp',
    'ic_read.hpp',
    'metadata.hpp',
    'newtonian_pm.hpp',
    'output.hpp',
    'parser.hpp',
    'particle_handler.hpp',
    'Particles_gevolution.hpp',
    'prng_engine.hpp',
    'radiation.hpp',
    'real_type.hpp',
    'tools.hpp',
    'velocity.hpp'
    ])
gevolution_sources = files([
    'background.cpp',
    'debugger.cpp',
    'gevolution.cpp',
    'hibernation.cpp',
    'ic_basic.cpp',
    'ic_prevolution.cpp',
    'ic_read.cpp',
    'newtonian_pm.cpp',
    'output.cpp',
    'parser.cpp',
    'Particles_gevolution.cpp',
    'radiation.cpp',
    'tools.cpp',
    'velocity.cpp'
    ])

config=configure_file(
    input: 'config.h.in', 
    output: 'config.h',
    configuration: conf_data)

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

main = files('main.cpp')
gevolution = executable('gevolution',
    main,version,
    link_with: libgevolution,
    dependencies: deps)
