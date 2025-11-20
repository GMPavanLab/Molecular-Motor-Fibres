#!/bin/bash
SCALE=0.95
AXIS=xyz
RES=COR

while [[ $# -gt 0 ]]; do
    case $1 in 
        -s|--scale)
            shift
            SCALE=$1
            shift
            ;;
        -ax|--axis)
            shift
            AXIS=$1
            shift
            ;;
        -f|-file_name)
            shift
            FILE_PATH=$1
            shift
            ;;
        -on|--residue)
            shift
            RES=$1
            shift
            ;;
        -*|--*)
            echo "Unknown option $1"
            shift;shift;
            ;;
        esac
done

cp $FILE_PATH tmp.gro
box=$(tail -n 1 $FILE_PATH)
x=$(awk '{print $1}' <<< $box)
y=$(awk '{print $2}' <<< $box)
z=$(awk '{print $3}' <<< $box)

if [[ $AXIS == *'x'* ]];then 
     awk -v factor=$SCALE 'NR > 2 {xcoord = substr($0,21,8); xscaled = sprintf("%8.3f", xcoord * factor); $0=substr($0,1,20) xscaled substr($0,29);}{ print }' tmp.gro > tmp2.gro
    cp tmp2.gro tmp.gro
    x=$(echo ''"${x}"' * '"${SCALE}"'' | bc -l)
fi
if [[ $AXIS == *'y'* ]]; then 
     awk -v factor=$SCALE 'NR > 2 {ycoord = substr($0,29,8); yscaled = sprintf("%8.3f", ycoord * factor); $0=substr($0,1,28) yscaled substr($0,37);}{ print }' tmp.gro > tmp2.gro
    cp tmp2.gro tmp.gro
    y=$(echo ''"${y}"' * '"${SCALE}"'' | bc -l)
fi
if [[ $AXIS == *'z'* ]]; then
     awk -v factor=$SCALE 'NR > 2 {zcoord = substr($0,37,8); zscaled = sprintf("%8.3f", zcoord * factor); $0=substr($0,1,36) zscaled substr($0,45);}{ print }' tmp.gro > tmp2.gro
    cp tmp2.gro tmp.gro
    z=$(echo ''"${z}"' * '"${SCALE}"'' | bc -l)
fi
sed -i '$s/.*/'"$( printf "%10.5f%10.5f%10.5f\n" $x $y $z)"'/' tmp.gro

cp tmp.gro $FILE_PATH