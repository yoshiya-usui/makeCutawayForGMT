#!/bin/sh 

iter=$1
skip=0
rot=0

width=-5/5
size=0.05

Hrange=-10/10
Vrange=-0.5/5
AreaSize=15/-4.125
ah=5
fh=5
av=2
fv=2

gmt makecpt -Cwysiwyg -Z -I -T1/3/0.2 > Res.cpt
gmt makecpt -Chot -Z -I -T0/2/0.2 > Aniso.cpt
gmt makecpt -Cseis -Z -I -T0/60/0.2 > Ang.cpt

if [ $skip = 0 ]; then
  echo 2 > param.txt
  echo ${iter} >> param.txt
  echo 0 >> param.txt
  echo 0 0 0 >> param.txt
  echo ${rot} >> param.txt
  echo 1 >> param.txt
  echo 0 >> param.txt
  echo 0.1 >> param.txt
  ../makeCutawayForGMT param.txt -aniso
fi

# --- draw color section ---
gmt set FONT_ANNOT_PRIMARY 20p,Helvetica,black
gmt set FONT_HEADING 14p,Helvetica
gmt set FONT_LABEL 20p,Helvetica
gmt set MAP_FRAME_WIDTH 0.2c
gmt set MAP_TICK_LENGTH_PRIMARY 0.05c 
gmt set MAP_ANNOT_OFFSET_PRIMARY 0.2c
gmt set COLOR_NAN 255/255/255
gmt set MAP_FRAME_TYPE PLAIN
gmt set FORMAT_GEO_MAP dddF

black=0/0/0
gray=100/100/100
shift=-0.5
Boption=a${ah}f${fh}:"x(km)":/a${av}f${fv}:"Depth(km)":WSne
obsloc=../obs_loc.txt

in_file=rhoXX_GMT_iter${iter}.dat
ps_file=rhoXX_NS_iter${iter}.ps
gmt psxy ${in_file} -CRes.cpt -JX$AreaSize -R$Hrange/$Vrange -B$Boption -P -X4 -Y8 -K -L > $ps_file
awk '{print $1, -0.2}' ${obsloc} | gmt psxy -Si0.3 -G${gray} -JX -R -X0 -Y0 -O -K >> $ps_file
gmt psscale -CRes.cpt -Dx10/-1.5/7.0/0.5h -N  -Bf1a1:"Resistivity [log(@~W@~m)]": -X0 -Y-2 -O >> $ps_file
rm ${in_file}

in_file=rhoYY_GMT_iter${iter}.dat
ps_file=rhoYY_NS_iter${iter}.ps
gmt psxy ${in_file} -CRes.cpt -JX$AreaSize -R$Hrange/$Vrange -B$Boption -P -X4 -Y8 -K -L > $ps_file
awk '{print $1, -0.2}' ${obsloc} | gmt psxy -Si0.3 -G${gray} -JX -R -X0 -Y0 -O -K >> $ps_file
gmt psscale -CRes.cpt -Dx10/-1.5/7.0/0.5h -N  -Bf1a1:"Resistivity [log(@~W@~m)]": -X0 -Y-2 -O -P >> $ps_file
rm ${in_file}

in_file=anisotropy_GMT_iter${iter}.dat
ps_file=anisotropy_NS_iter${iter}.ps
gmt psxy ${in_file} -CAniso.cpt -JX$AreaSize -R$Hrange/$Vrange -B$Boption -P -X4 -Y8 -K -L > $ps_file
awk '{print $1, -0.2}' ${obsloc} | gmt psxy -Si0.3 -G${gray} -JX -R -X0 -Y0 -O -K >> $ps_file
gmt psscale -CAniso.cpt -Dx10/-1.5/7.0/0.5h -N  -Bf1a1:"Anisotropy": -X0 -Y-2 -O -P >> $ps_file
rm ${in_file}

in_file=strike_GMT_iter${iter}.dat
ps_file=strike_NS_iter${iter}.ps
gmt psxy ${in_file} -CAng.cpt -JX$AreaSize -R$Hrange/$Vrange -B$Boption -P -X4 -Y8 -K -L > $ps_file
awk '{print $1, -0.2}' ${obsloc} | gmt psxy -Si0.3 -G${gray} -JX -R -X0 -Y0 -O -K >> $ps_file
gmt psscale -CAng.cpt -Dx10/-1.5/7.0/0.5h -N  -Bf15a15:"Strike [deg.]": -X0 -Y-2 -O -P >> $ps_file
rm ${in_file}

in_file=dip_GMT_iter${iter}.dat
ps_file=dip_NS_iter${iter}.ps
gmt psxy ${in_file} -CAng.cpt -JX$AreaSize -R$Hrange/$Vrange -B$Boption -P -X4 -Y8 -K -L > $ps_file
awk '{print $1, -0.2}' ${obsloc} | gmt psxy -Si0.3 -G${gray} -JX -R -X0 -Y0 -O -K >> $ps_file
gmt psscale -CAng.cpt -Dx10/-1.5/7.0/0.5h -N  -Bf15a15:"Dip [deg.]": -X0 -Y-2 -O -P >> $ps_file
rm ${in_file}

rm *.cpt
