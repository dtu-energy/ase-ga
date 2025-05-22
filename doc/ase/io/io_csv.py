# creates: io.csv
from ase.io.formats import all_formats, get_ioformat

with open('io.csv', 'w') as fd:
    print('format, description, capabilities', file=fd)
    for format in sorted(all_formats):
        io = get_ioformat(format)
        c = ''
        if io.can_read:
            c = 'R'
        if io.can_write:
            c += 'W'
        if not io.single:
            c += '+'
        print(f':ref:`{format}`, {all_formats[format][0]}, {c}', file=fd)
