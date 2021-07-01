#!/bin/csh 

set iter=4

set width=-5/5
set size=0.05

set Hrange=-85/110
set Vrange=-5/92.5
set AreaSize=15/-7.5
set ah=20
set fh=20
set av=20
set fv=20

makecpt -Chaxby -I -T0/4/0.2 > Rcpt

/home/yusui/001_Code/makeCutawayForGMT/makeCutawayForGMT param_V1.txt

# --- draw color section ---
gmtset ANOT_FONT Helvetica
gmtset HEADER_FONT Helvetica
gmtset LABEL_FONT Helvetica
gmtset ANOT_FONT_SIZE 20pt HEADER_FONT_SIZE 14pt LABEL_FONT_SIZE 20pt
gmtset FRAME_WIDTH 0.1c
gmtset TICK_LENGTH 0.05c ANOT_OFFSET 0.2c
gmtset COLOR_NAN 255/255/255
gmtset BASEMAP_TYPE PLAIN
gmtset DEGREE_FORMAT 100

set in_file=resistivity_GMT_iter${iter}.dat
set ps_file=resistivity_GMT.iter${iter}.V1.ps
echo $in_file
echo $ps_file

set shift=-0.5
set Boption=a${ah}f${fh}:"y(km)":/a${av}f${fv}:"z(km)":WSne

psxy $in_file -CRcpt -JX$AreaSize -R$Hrange/$Vrange -B$Boption -M -P -X4 -Y8 -K -V -L >! $ps_file
awk ' $1 > -40 && $1 < -35 {print $2, $3}' loc.txt | psxy -J -R -G0/0/0 -Ss6p -O -K -P >> $ps_file

psscale -CRcpt -D10/-1.5/7.0/0.5h -N  -Bf1a1:"Resistivity [log(@~W@~m)]": -U/0/-5/$ps_file -X2 -Y-2 -O -V >> $ps_file

#-- remove temp files --
rm -f Rcpt

