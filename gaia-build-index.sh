#!/bin/bash
# Gaia Data Release 1 (for FITS table) healpixes => index files script for astrometry.net
# Developer: Yücel KILIÇ
# Usage: $ ./gaia-build-index.sh healpix_folder/ index_output/ N_pix_count series_name
# N_pix => healpixes count
# Bash Progress Bar: https://gist.github.com/F1LT3R/fa7f102b08a514f2c535
# Thanks for the progress bar function ^.

progressBarWidth=20

# Function to draw progress bar
progressBar () {

    # Calculate number of fill/empty slots in the bar
    progress=$(echo "$progressBarWidth/$taskCount*$tasksDone" | bc -l)
    fill=$(printf "%.0f\n" $progress)
    if [ $fill -gt $progressBarWidth ]; then
        fill=$progressBarWidth
    fi
    empty=$(($fill-$progressBarWidth))

    # Percentage Calculation
    percent=$(echo "100/$taskCount*$tasksDone" | bc -l)
    percent=$(printf "%0.2f\n" $percent)
    if [ $(echo "$percent>100" | bc) -gt 0 ]; then
        percent="100.00"
    fi

    # Output to screen
    printf "\r["
    printf "%${fill}s" '' | tr ' ' '#'
    printf "%${empty}s" '' | tr ' ' '-'
    printf "] $percent%% - $text "
}

series=$4
sn=${series:0:2}

if [ $3 == 48 ]; then
    nside=2
else
    nside=1
fi

## Collect task count
taskCount=$(( $3 * 5 ))
tasksDone=0

for ((HP=0; HP<$3; HP++)); do
    for ((SCALE=0; SCALE<=4; SCALE++)); do
        HH=$(printf %02i $HP)
        SS=$(printf %02i $SCALE)
        build-astrometry-index -i $1/gaia-hp${HH}.fits -s $nside -H $HP -P $SCALE -E -S phot_g_mean_mag -o $2/index-$sn${SS}-${HH}.fits -I $sn${SS}-${HH}

        # Task done and progress the bar...
        (( tasksDone += 1 ))

        # Add some friendly output
        text=$(echo "index-$sn${SS}-${HH}.fits is indexed...")

        # Draw the progress bar
        progressBar $taskCount $taskDone $text
        
        sleep 0.01
    done
done
