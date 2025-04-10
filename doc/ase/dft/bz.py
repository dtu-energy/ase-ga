# creates: bztable.rst
# creates: 00.CUB.svg 01.FCC.svg 02.BCC.svg 03.TET.svg 04.BCT1.svg
# creates: 05.BCT2.svg 06.ORC.svg 07.ORCF1.svg 08.ORCF2.svg 09.ORCF3.svg
# creates: 10.ORCI.svg 11.ORCC.svg 12.HEX.svg 13.RHL1.svg 14.RHL2.svg
# creates: 15.MCL.svg 16.MCLC1.svg 17.MCLC3.svg 18.MCLC5.svg 19.TRI1a.svg
# creates: 20.TRI1b.svg 21.TRI2a.svg 22.TRI2b.svg
# creates: 23.OBL.svg 24.RECT.svg 25.CRECT.svg 26.HEX2D.svg 27.SQR.svg

from ase.lattice import all_variants

header = """\

Brillouin zone data
-------------------

.. list-table::
    :widths: 10 15 45
"""


entry = """\
    * - {name} ({longname})
      - {bandpath}
      - .. image:: {fname}
            :width: 40 %
"""

with open('bztable.rst', 'w') as fd:
    print(header, file=fd)

    for i, lat in enumerate(all_variants()):
        id = f'{i:02d}.{lat.variant}'
        imagefname = f'{id}.svg'
        txt = entry.format(
            name=lat.variant,
            longname=lat.longname,
            bandpath=lat.bandpath().path,
            fname=imagefname,
        )
        print(txt, file=fd)
        ax = lat.plot_bz()
        fig = ax.get_figure()
        fig.savefig(imagefname, bbox_inches='tight')
        fig.clear()
