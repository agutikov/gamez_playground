project('gamez-playground', 'c', 'cpp',
  default_options: ['c_std=c11', 'cpp_std=c++20']
)

# From https://github.com/gilzoide/raylib-meson-template



cc = meson.get_compiler('c')
cxx = meson.get_compiler('cpp')

c_args = []
cpp_args = []
link_args = ['-lfmt']
name_suffix = []
dependencies = []



# Raylib configuration
cmake = import('cmake')
raylib_opts = cmake.subproject_options()
raylib_opts.set_install(false)
raylib_opts.add_cmake_defines({
  'BUILD_EXAMPLES': 'OFF',
  # Add other CMake options for Raylib, if needed
})
if cc.get_id() == 'emscripten'
  raylib_opts.add_cmake_defines({
    'PLATFORM': 'Web',
  })
else
  raylib_opts.add_cmake_defines({
    'PLATFORM': 'Desktop',
    'USE_EXTERNAL_GLFW': 'ON'
  })
endif
raylib_subproject = cmake.subproject('raylib', options: raylib_opts)
dependencies += raylib_subproject.dependency('raylib')


# GLFW configuration
#glfw_opts = cmake.subproject_options()
#glfw_opts.set_install(false)
#glfw_subproject = cmake.subproject('glfw', options: glfw_opts)
#dependencies += glfw_subproject.dependency('glfw3')


# General configuration
if host_machine.system() == 'windows'
  dependencies += [
    cc.find_library('winmm'),
  ]
elif host_machine.system() == 'darwin'
  # TODO
elif host_machine.system() == 'linux'
  dependencies += [
    cc.find_library('m'),
    cc.find_library('GL'),
    cc.find_library('dl'),
    cc.find_library('rt'),
    cc.find_library('X11')
  ]
elif host_machine.system() == 'emscripten'
  link_args += [
    '-s', 'ENVIRONMENT=web',
    '-s', 'USE_GLFW=3',
  ]
  name_suffix = 'html'
endif


inc_dir = include_directories(
  'src',
  'subprojects/raygui/src',
  'subprojects/PerlinNoise'
)

sources = files(
  'src/core_2d_camera_platformer.c',
)

executable('playground',
  sources,
  c_args: c_args,
  cpp_args: cpp_args,
  dependencies: dependencies,
  include_directories: inc_dir,
  link_args: link_args,
  name_suffix: name_suffix,
)

particles_sources = files(
  'src/2d_particles.cpp',
)

executable('2d_particles',
  particles_sources,
  c_args: c_args,
  cpp_args: cpp_args,
  dependencies: dependencies,
  include_directories: inc_dir,
  link_args: link_args,
  name_suffix: name_suffix,
)

space_gen_sources = files(
  'src/2d_space_gen.cpp',
)

executable('2d_space_gen',
  space_gen_sources,
  c_args: c_args,
  cpp_args: cpp_args,
  dependencies: dependencies,
  include_directories: inc_dir,
  link_args: link_args,
  name_suffix: name_suffix,
)

