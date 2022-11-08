db_test = executable('db_test',
    'db_test.cpp',
    link_with: libgevolution,
    dependencies: deps,
    include_directories: include)
    
test('Distributed database',mpirun,
   args: ['--oversubscribe','-np','4',db_test.full_path()],
   is_parallel: false)
