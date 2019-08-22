project('optflow', 'c',
  default_options : ['default_library=static', 'buildtype=release', 'warning_level=0',
    'c_std=gnu99'])

cc = meson.get_compiler('c')
if cc.get_id() == 'gcc'
  add_project_arguments('-Wno-unused-result', '-Wno-implicit-int', language: 'c')
endif

os = host_machine.system()

if not ['msvc', 'intel-cl']. contains(cc.get_id()) and not (os == 'windows' and cc.get_id() == 'pgi')
  math_lib = cc.find_library('m')
else
  math_lib = []
endif


if not (cc.get_id() == 'msvc') and add_languages('fortran', required: false)
  fc = meson.get_compiler('fortran')

  blas = dependency('blas', required: false)
  if not blas.found()
    blas = fc.find_library('blas', required: false)
  endif
  if not blas.found()
    blas_proj = subproject('blas')
    blas = declare_dependency(link_with: blas_proj.get_variable('blas'))
  endif

  linpack = library('linpack', 'SOFTWARE/dsvdc.f',
    dependencies: blas)

  fleet_exe = executable('fleet', 'FLEET/fleet.c',
    dependencies: math_lib,
    link_with: linpack)
  test('Fleet', fleet_exe,
    args: ['new2rastertreet.', '2.', '9', '0.05', meson.source_root() / 'TESTDATA/TREE_DATA/TRANS', meson.source_root() / 'out', meson.source_root() / 'out'],
    timeout: 30,
    is_parallel: false)

  singh_exe = executable('singh', 'SINGH/singh.c',
    link_with: linpack,
    dependencies: math_lib,
    include_directories: 'SOFTWARE')
  test('Singh', singh_exe,
    args: ['new2binarytreet.', '20', meson.source_root() / 'TESTDATA/TREE_DATA/TRANS', meson.source_root() /'out', '-B', '150', '150', '-i', '25', '-n', '2', '-w', '2', '-N', '4'],
    timeout: 60)
else
  fc = disabler()
endif

flow2ps_exe = executable('flow2ps', 'SOFTWARE/psflow.c',
  dependencies: math_lib)

anandan_exe = executable('anandan', 'ANANDAN/anandan.c',
  dependencies: math_lib)
test('Anandan', anandan_exe,
  args: [meson.source_root() / 'TESTDATA/TREE_DATA/TRANS',
         meson.source_root() / 'out/',
         '/new2binarytreet.', '-N', '9', '10', '-W', '5', '-L', '3', '-I', '15'])

horn_exe = executable('horn', 'HORN/horn.c',
  dependencies: math_lib,
  include_directories: 'SOFTWARE')
test('HornSchunck', horn_exe,
  args: ['new2binarytreet.', '0.5', '1.5', '20', '100', meson.source_root() / 'TESTDATA/TREE_DATA/TRANS', meson.source_root() / 'out/', '-B', '150', '150', '-MH', '-T', '5.0']
)

lucas_exe = executable('lucas', 'LUCAS/lucas.c',
  dependencies: math_lib,
  include_directories: 'SOFTWARE')
test('LucasKanade', lucas_exe,
  args: ['new2binarytreet.', '1.5', '20', '1.0', meson.source_root() / 'TESTDATA/TREE_DATA/TRANS', meson.source_root() / 'out/', '-B', '150', '150']
)

nagel_exe = executable('nagel', 'NAGEL/nagel.c',
  dependencies: math_lib)
# test('Nagel', nagel_exe,
#  ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS
#  ${CMAKE_SOURCE_DIR}/out/
#  newbinarytreet. -SG 1.5 -TG 1.5 -M 20 -A 0.5 -I 100 -B 150 150 -F 5.0)

uras_exe = executable('uras', 'URAS/uras.c',
  dependencies: math_lib)
#test('Uras', uras_exe,
#  ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS
#  ${CMAKE_SOURCE_DIR}/out/
#  new2binarytreet. -SG 3.0 -TG 1.5 -M 9 -R 1 -F 0.5)
