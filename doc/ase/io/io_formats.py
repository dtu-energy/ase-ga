# creates: formatoptions.rst
from ase.io.formats import all_formats

with open('formatoptions.rst', 'w') as fd:
    print('============================================', file=fd)
    print('Format Specific Options', file=fd)
    print('============================================', file=fd)
    format_names = sorted(all_formats.keys())
    for name in format_names:
        fmt = all_formats[name]
        print(f".. _{name}:\n", file=fd)
        print(name, file=fd)
        print('----------------------------------------', file=fd)
        if fmt.can_read:
            print('\n .. autofunction:: {:}.read_{:}\n\n'.format(
                fmt.module_name, fmt._formatname), file=fd)
        if fmt.can_write:
            print('\n .. autofunction:: {:}.write_{:}\n\n'.format(
                fmt.module_name, fmt._formatname), file=fd)
        if (not fmt.can_read) and (not fmt.can_write):
            print('\n No automatic documentation of this module found.',
                  file=fd)
            print(
                '\n Visit '
                'https://gitlab.com/ase/ase/tree/master/ase/io/{:}.py'
                ' to see if you find the information needed in '
                'the source code.\n\n'.format(
                    fmt.module_name),
                file=fd)
