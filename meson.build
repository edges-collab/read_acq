project('read-acq', 'c',
    version: run_command('python', '-m', 'setuptools_scm').stdout().strip(),
)

py = import('python').find_installation(pure: false)

py.extension_module(
    'decode',
    'src/read_acq/decode.c',
    install: true,
)