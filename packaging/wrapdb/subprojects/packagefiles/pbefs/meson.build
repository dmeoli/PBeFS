project(
  'pbefs',
  'cpp',
  version: '0.1.0',
  license: 'LGPL-3.0-only',
  meson_version: '>=1.1.0',
  default_options: ['cpp_std=c++20'],
)

# the library with the C++ threads; the runtimes of OpenMP, MPI and FastFlow
# are headers that a project takes with its own dependency on them. MSVC
# exports no symbol unless told to, and hence there the library is static
cpp = meson.get_compiler('cpp')
pbefs_type = cpp.get_argument_syntax() == 'msvc' ? 'static_library' : 'library'
pbefs_inc = include_directories('include')
threads_dep = dependency('threads')
pbefs_lib = build_target(
  'PBeFS',
  'src/Distribution.cpp',
  'src/Grid.cpp',
  'src/SearchUtils.cpp',
  'src/Solution.cpp',
  target_type: pbefs_type,
  include_directories: pbefs_inc,
  dependencies: threads_dep,
  install: true,
)
pbefs_dep = declare_dependency(
  link_with: pbefs_lib,
  include_directories: pbefs_inc,
  dependencies: threads_dep,
)
meson.override_dependency('pbefs', pbefs_dep)

install_subdir(
  'include/PBeFS',
  install_dir: get_option('includedir'),
)
import('pkgconfig').generate(
  pbefs_lib,
  name: 'pbefs',
  description: 'Parallel best-first search: A* and branch-and-bound',
  url: 'https://github.com/dmeoli/PBeFS',
)

test(
  'use-meson',
  executable(
    'use-meson',
    'packaging/test/main.cpp',
    dependencies: pbefs_dep,
  ),
)
