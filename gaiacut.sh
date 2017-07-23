#!/bin/bash
# Gaia Data Release 1 (for FITS table) catalogue column cut script for astrometry.net
# Developer: Yücel KILIÇ
# Usage: $ ./gaiacut.sh gaia_dr1_cat_path/ output_cat/ ra;ra_error;dec;dec_error;pmra;pmra_error;pmdec;pmdec_error;phot_g_mean_mag -c "[ra_error <= 1 && dec_error <= 1 && phot_g_mean_mag >=15 && phot_g_mean_mag <=18]"

time for file in $(ls $1/*.fits);
do
    fpath=$file
    basename=${fpath##*/}    
    echo "$file ==> $2/$basename"
    if [ ! -z "$4" ]; then
        time fitscopy $file"[1][col $3]" $2/$basename
    elif [ -z "$4" ] && [ "$4" == "-c" ]; then
        time fitscopy $file"[1][col $3]$5" $2/$basename
    fi
    du -h $2/$basename
done
