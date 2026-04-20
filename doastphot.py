from astrolib import astronomy
from astrolib import photometry
from astrolib import visuals
import glob
import sys
import os

print("""
DATE/
|    \
BDF/  SCI_IMAGES

Örnek: python3 doastphot.py DATE/SCI_IMAGES/ R --skip-calib

""")

fitsdir = sys.argv[1]
filter = sys.argv[2]
# Ön indirgeme
if filter != "--plot-lc":
    target = os.path.basename(fitsdir).split('.')[0]

    if "--skip-calib" not in sys.argv:
        print("Ön indirgeme başladı!")
        ro = astronomy.RedOps()
        ro.ccdproc(fitsdir, filter=filter)
    
    # astrometry
    # Solve field with astrometry.net
    # Please provide image path, so ./ is important!
    if "--skip-astrometry" not in sys.argv:
        ac = astronomy.AstCalc()
        print("Astrometry başladı!")
        fitsfiles = glob.glob("./atmp/bf_*.fit?")
    
        for fitsfile in fitsfiles:
            ac.solve_field(fitsfile)
    # photometry
    if "--skip-photometry" not in sys.argv:
        ap = photometry.PhotOps()
        ap.asteroids_phot("./atmp/bf_*new.fits", multi_object=True, radius=11)
        # For plotting ligt curve of results.
# plot-lc
elif filter == "--plot-lc":
    # fitsdir burada dosya sonuç soyasını tutan değişken!
    plc = visuals.StarPlot()
    plc.lc_plot(fitsdir)
