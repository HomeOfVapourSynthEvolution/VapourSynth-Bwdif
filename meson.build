project('Bwdif', 'cpp',
  default_options: ['buildtype=release', 'warning_level=2', 'b_lto=true', 'b_ndebug=if-release', 'cpp_std=c++17'],
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.51.0',
  version: '4.1'
)

cxx = meson.get_compiler('cpp')

gcc_syntax = cxx.get_argument_syntax() == 'gcc'

if get_option('buildtype') == 'release'
  add_project_arguments(gcc_syntax ? ['-fno-math-errno', '-fno-trapping-math'] : '/GS-', language: 'cpp')
endif

if gcc_syntax
  vapoursynth_dep = dependency('vapoursynth').partial_dependency(compile_args: true, includes: true)
  install_dir = vapoursynth_dep.get_variable(pkgconfig: 'libdir') / 'vapoursynth'
else
  vapoursynth_dep = dependency('', required: false)
  install_dir = get_option('libdir') / 'vapoursynth'
endif

sources = [
  'Bwdif/Bwdif.cpp',
  'Bwdif/Bwdif.h'
]

libs = []

if host_machine.cpu_family().startswith('x86')
  project_args = ['-DBWDIF_X86']
  if gcc_syntax
    project_args += ['-mfpmath=sse', '-msse2']
  endif
  add_project_arguments(project_args, language: 'cpp')

  sources += [
    'Bwdif/Bwdif_SSE2.cpp',
    'Bwdif/VCL2/instrset.h',
    'Bwdif/VCL2/instrset_detect.cpp',
    'Bwdif/VCL2/vector_convert.h',
    'Bwdif/VCL2/vectorclass.h',
    'Bwdif/VCL2/vectorf128.h',
    'Bwdif/VCL2/vectorf256.h',
    'Bwdif/VCL2/vectorf256e.h',
    'Bwdif/VCL2/vectorf512.h',
    'Bwdif/VCL2/vectorf512e.h',
    'Bwdif/VCL2/vectori128.h',
    'Bwdif/VCL2/vectori256.h',
    'Bwdif/VCL2/vectori256e.h',
    'Bwdif/VCL2/vectori512.h',
    'Bwdif/VCL2/vectori512e.h',
    'Bwdif/VCL2/vectori512s.h',
    'Bwdif/VCL2/vectori512se.h',
    'Bwdif/VCL2/vectormath_common.h',
    'Bwdif/VCL2/vectormath_exp.h',
    'Bwdif/VCL2/vectormath_hyp.h',
    'Bwdif/VCL2/vectormath_lib.h',
    'Bwdif/VCL2/vectormath_trig.h'
  ]

  libs += static_library('avx2', 'Bwdif/Bwdif_AVX2.cpp',
    cpp_args: gcc_syntax ? ['-mavx2', '-mfma'] : '/arch:AVX2',
    dependencies: vapoursynth_dep,
    gnu_symbol_visibility: 'hidden'
  )

  libs += static_library('avx512', 'Bwdif/Bwdif_AVX512.cpp',
    cpp_args: gcc_syntax ? ['-mavx512f', '-mavx512bw', '-mavx512dq', '-mavx512vl', '-mfma'] : '/arch:AVX512',
    dependencies: vapoursynth_dep,
    gnu_symbol_visibility: 'hidden'
  )
endif

shared_module('bwdif', sources,
  dependencies: vapoursynth_dep,
  link_with: libs,
  install: true,
  install_dir: install_dir,
  gnu_symbol_visibility: 'hidden'
)
