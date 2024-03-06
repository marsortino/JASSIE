#!/bin/bash
FILENAME="log_compiler.txt"

echo 'Please make sure the terminal is in the same directory of compiler.sh'
echo 'Starting...'
date +%F_%T >> "$FILENAME"
echo 'Starting...' >> "$FILENAME"
gcc -g -Wall lib/c_routines/codes/raytrace.c -o raytrace -lm >> "$FILENAME" 2>&1

if [ -f raytrace ]; then
    mv raytrace lib/c_routines/exec/
    echo 'Compile successful.'
    echo 'Compile successful.' >> "$FILENAME"
else
    echo 'Compile failed. Check' "$FILENAME"
fi


echo '////////**********////////' >> "$FILENAME"