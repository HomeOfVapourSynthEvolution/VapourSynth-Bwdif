project('Bwdif', 'cpp',
  default_options: ['buildtype=release', 'b_ndebug=if-release', 'cpp_std=c++17'],
  meson_version: '>=0.48.0',
  version: '2'
)

sources = [
  'Bwdif/Bwdif.cpp',
  'Bwdif/Bwdif.h'
]

vapoursynth_dep = dependency('vapoursynth').partial_dependency(compile_args: true, includes: true)

libs = []

if host_machine.cpu_family().startswith('x86')
  add_project_arguments('-DBWDIF_X86', '-mfpmath=sse', '-msse2', language: 'cpp')

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
    dependencies: vapoursynth_dep,
    cpp_args: ['-mavx2', '-mfma'],
    gnu_symbol_visibility: 'hidden'
  )

  libs += static_library('avx512', 'Bwdif/Bwdif_AVX512.cpp',
    dependencies: vapoursynth_dep,
    cpp_args: ['-mavx512f', '-mavx512bw', '-mavx512dq', '-mavx512vl', '-mfma'],
    gnu_symbol_visibility: 'hidden'
  )
endif

shared_module('bwdif', sources,
  dependencies: vapoursynth_dep,
  link_with: libs,
  install: true,
  install_dir: join_paths(vapoursynth_dep.get_pkgconfig_variable('libdir'), 'vapoursynth'),
  gnu_symbol_visibility: 'hidden'
)
