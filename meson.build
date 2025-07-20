project('read-acq', 'c',
    version: run_command('python', '-m', 'setuptools_scm').stdout().strip(),
)

py = import('python').find_installation(pure: false)

# We need -lm for all C code (assuming it uses math functions, which is safe to
# assume for SciPy). For C++ it isn't needed, because libstdc++/libc++ is
# guaranteed to depend on it. For Fortran code, Meson already adds `-lm`.
cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

py.install_sources(
    [
        'src/read_acq/__init__.py',
        'src/read_acq/_coordinates.py',
        'src/read_acq/codec.py', 
        'src/read_acq/gsdata.py', 
        'src/read_acq/cli.py',
        'src/read_acq/read_acq.py',
    ],
    subdir: 'read_acq'
)
py.extension_module(
    'decode',
    'src/read_acq/decode.c',
    install: true,
    subdir: 'read_acq',
)