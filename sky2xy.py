from astropy.coordinates import SkyCoord
from astropy import units as u
from astrolib import astronomy
import glob
import sys
import os
import numpy as np


# Görüntünün yolu
image_path = sys.argv[1]

# Kullanıcıdan silitin içine getirilecek hedef cismin koordinatları
ra = sys.argv[2]
dec = sys.argv[3]

# astrometry yapan astrolib modülü
ac = astronomy.AstCalc()
print("Astrometry başladı!")
ac.solve_field(image_path, ra_keyword="RA", dec_keyword="DEC")

# WCS çözümlü dosyanın yeni ismi.
image_with_wcs = "{0}_new.fits".format(os.path.splitext(image_path)[0])

# Kullanıcı tarafından talep edilen RA, DEC'in X, Y karşılığı
# astrolib modülünden
xy = ac.sky2xy(image_with_wcs, ra, dec)

# Seçilen cismin frame dışından olmadığının denetimi
fitsops = astronomy.FitsOps(image_with_wcs)
naxis1 = fitsops.get_header("naxis1")
naxis2 = fitsops.get_header("naxis2")

# string formatının işlem yapılabilir formata dönüştürülmesi
target_sky_coor = SkyCoord('{0} {1}'.format(ra, dec),
                           frame='icrs',
                           unit=(u.hourangle, u.deg))

# Hedef cisimin, merkez veya belirlenen bir noktaya uzaklığı
# Burada eğer X, Y kullanıcı tarafından özel olarak belirtilmezse
# hedef cismin görüntünün merkezine olan uzaklığı hesaplanıyor.
# ilk koşul X, Y komut satırında belirtilirse
if len(sys.argv) > 4:
    users_x, users_y = float(sys.argv[4]), float(sys.argv[5])

    if naxis1 < users_x or naxis2 < users_y or users_x < 0 or users_y < 0:
        print("Provided coordinates are out of frame!")
        raise SystemExit
    else:
        users_radec = ac.xy2sky2(image_with_wcs, users_x, users_y)
        users_sky_coor = SkyCoord(users_radec.ra.value,
                                  users_radec.dec.value,
                                  frame='icrs', unit='deg')
        # Hedef cismin kullanıcının talep ettiği X, Y noktasına uzaklığı
        dra, ddec = target_sky_coor.spherical_offsets_to(users_sky_coor)
else:
    # Eğer X, Y kullanıcı tarafından belirtilmezse
    # Hedef cisimin merkeze olan uzaklığı hesaplanıyor
    c_x, c_y = [float(naxis1) / 2, float(naxis2) / 2]

    image_center = ac.center_finder(image_with_wcs, wcs_ref=True)

    center_sky_coor = SkyCoord(image_center[0].value,
                               image_center[1].value, frame='icrs', unit='deg')
    # Hedef cismin merkeze olan uzaklığı
    dra, ddec = target_sky_coor.spherical_offsets_to(center_sky_coor)

if xy:
    print("X: {0}, Y: {1}".format(xy[0], xy[1]))
    # print("cX: {0}, cY: {1}".format(c_x, c_y))
    # print("dX: {0}, dY: {1}".format(float(xy[0] - c_x),
    # float(xy[1] - c_y)))
    print("{0}:{1}".format(dra.to(u.arcsec).value,
                           ddec.to(u.arcsec).value))