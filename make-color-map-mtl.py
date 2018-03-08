#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

def main():
    import sys
    import matplotlib.pylab as plt

    if len(sys.argv) not in (2, 3):
        print("usage: %s color-map-name [count]"%sys.argv[0])
        sys.exit(0)

    name = sys.argv[1]
    count = int(sys.argv[2] if len(sys.argv) == 3 else 256)
    cmap = plt.cm.get_cmap(name, count)
    with open('%s.mtl'%name, 'w') as f:
        print('# %s colormap with %d colors'%(name, count), file=f)
        print('', file=f)
        for i in range(cmap.N):
            rgb = cmap(i)[:3] # get RGB from RGBA
            print("""newmtl {name}{idx}
Ka {r:.3f} {g:.3f} {b:.3f}
Kd {r:.3f} {g:.3f} {b:.3f}
Ks 0.000 0.000 0.000
d 0.75
illum 2
""".format(name=name, idx=i, r=rgb[0], g=rgb[1], b=rgb[2]), file=f)

if __name__ == "__main__": main()
