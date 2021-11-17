#!/bin/sh 

iter=4

width=-5/5
size=0.05

Hrange=-600/600
Vrange=-5/800
AreaSize=15/-10
ah=200
fh=200
av=100
fv=100

gmt makecpt -Chaxby -Z -I -T0/4/0.2 > Res.cpt

./makeCutawayForGMT param_V1.txt

# --- draw color section ---
#gmt set ANOT_FONT Helvetica
gmt set FONT_ANNOT_PRIMARY 20p,Helvetica,red
#gmt set HEADER_FONT Helvetica
gmt set FONT_HEADING 14p,Helvetica
#gmt set LABEL_FONT Helvetica
gmt set FONT_LABEL 20p,Helvetica
#gmt set ANOT_FONT_SIZE 20pt HEADER_FONT_SIZE 14pt LABEL_FONT_SIZE 20pt
gmt set MAP_FRAME_WIDTH 0.1c
gmt set MAP_TICK_LENGTH_PRIMARY 0.05c 
gmt set MAP_ANNOT_OFFSET_PRIMARY 0.2c
gmt set COLOR_NAN 255/255/255
gmt set MAP_FRAME_TYPE PLAIN
#gmt set DEGREE_FORMAT 100
gmt set FORMAT_GEO_MAP dddF
#set in_file= resistivity_GMT_iter${iter}.dat
#set resistivity_GMT.iter${iter}.H.ps
in_file=resistivity_GMT_iter13.dat
ps_file=resistivity_GMT.iter13.V1.ps
echo $in_file
echo $ps_file

shift=-0.5
Boption=a${ah}f${fh}:"y(km)":/a${av}f${fv}:"z(km)":WSne

grep Z $in_file | awk '{print $3}' > value.txt 
awk '{if ($1 == ">") print $0; else printf "%15.6e%15.6e\n", -$1, $2 }' $in_file > tmp.txt
gmt psxy tmp.txt -CRes.cpt -G+z -Zvalue.txt -JX$AreaSize -R$Hrange/$Vrange -B$Boption -P -X4 -Y8 -K -V -L > $ps_file
gmt psscale -CRes.cpt -Dx10/-1.5/7.0/0.5h -N  -Bf1a1:"Resistivity [log(@~W@~m)]": -U/0/-5/$ps_file -X2 -Y-2 -O -P -V >> $ps_file

# -- remove temp files --
rm -f Res.cpt
rm -f value.txt
rm -f tmp.txt
