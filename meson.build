project('optflow', 'c', 'fortran',
  default_options : ['default_library=static', 'buildtype=release', 'warning_level=0',
    'c_std=gnu99'])

cc = meson.get_compiler('c')
fc = meson.get_compiler('fortran')

if cc.get_id() == 'gcc'
  add_project_arguments('-Wno-unused-result', '-Wno-implicit-int', language: 'c')
endif

libm = cc.find_library('m')
blas = dependency('blas', required: false)
if not blas.found()
  blas = fc.find_library('blas', required: false)
endif

if blas.found()
  linpack = library('linpack', 'SOFTWARE/dsvdc.f',
    dependencies: blas)

  fleet_exe = executable('fleet', 'FLEET/fleet.c',
    dependencies: libm,
    link_with: linpack)
  test('Fleet', fleet_exe,
    args: ['new2rastertreet.', '2.', '9', '0.05', meson.source_root() / 'TESTDATA/TREE_DATA/TRANS', meson.source_root() / 'out', meson.source_root() / 'out'],
    timeout: 30,
    is_parallel: false)

  singh_exe = executable('singh', 'SINGH/singh.c',
    link_with: linpack,
    dependencies: libm,
    include_directories: 'SOFTWARE')

  test('Singh', singh_exe,
    args: ['new2binarytreet.', '20', meson.source_root() / 'TESTDATA/TREE_DATA/TRANS', meson.source_root() /'out', '-B', '150', '150', '-i', '25', '-n', '2', '-w', '2', '-N', '4'],
    timeout: 60)
endif

flow2ps_exe = executable('flow2ps', 'SOFTWARE/psflow.c',
  dependencies: libm)

anandan_exe = executable('anandan', 'ANANDAN/anandan.c',
  dependencies: libm)
# test('Anandan', anandan_exe
#  ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS
#  ${CMAKE_SOURCE_DIR}/out/
#  new2binarytreet. -N 9 10 -W 5 -L 3 -I 15)

horn_exe = executable('horn', 'HORN/horn.c',
  dependencies: libm,
  include_directories: 'SOFTWARE')
test('HornSchunck', horn_exe,
  args: ['new2binarytreet.', '0.5', '1.5', '20', '100', meson.source_root() / 'TESTDATA/TREE_DATA/TRANS', meson.source_root() / 'out/', '-B', '150', '150', '-MH', '-T', '5.0']
)

lucas_exe = executable('lucas', 'LUCAS/lucas.c',
  dependencies: libm,
  include_directories: 'SOFTWARE')
test('LucasKanade', lucas_exe,
  args: ['new2binarytreet.', '1.5', '20', '1.0', meson.source_root() / 'TESTDATA/TREE_DATA/TRANS', meson.source_root() / 'out/', '-B', '150', '150']
)

nagel_exe = executable('nagel', 'NAGEL/nagel.c',
  dependencies: libm)
# test('Nagel', nagel_exe,
#  ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS
#  ${CMAKE_SOURCE_DIR}/out/
#  newbinarytreet. -SG 1.5 -TG 1.5 -M 20 -A 0.5 -I 100 -B 150 150 -F 5.0)

uras_exe = executable('uras', 'URAS/uras.c',
  dependencies: libm)
#test('Uras', uras_exe,
#  ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS
#  ${CMAKE_SOURCE_DIR}/out/
#  new2binarytreet. -SG 3.0 -TG 1.5 -M 9 -R 1 -F 0.5)
