#!/bin/bash

# This file is part of the RHybrid simulation platform.
#
# Copyright 2018- Aalto University
# Copyright 2014- Finnish Meteorological Institute
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# check needed commands are available

if ! [ -x "$(command -v /bin/ls)" ]; then
 echo "ERROR: /bin/ls not found"
 exit 1
fi

if command -v rhb_log2pdf.sh >/dev/null 2>&1; then
 LOG_CMD="rhb_log2pdf.sh"
elif command -v ./rhb_log2pdf.sh >/dev/null 2>&1; then
 LOG_CMD="./rhb_log2pdf.sh"
else
 echo "ERROR: rhb_log2pdf.sh not in path and ./rhb_log2pdf.sh not found"
 exit 1
fi

# multiple dirs
if [ "$#" -le "0" ]; then
 echo "USAGE: rhb_log2pdf_all.sh lin/log [dir1/] [dir2/] [dir3/] [...]"
elif [ "$#" -eq "1" ]; then
 sdirs=""
 firstdir=./
else
 sdirs=$(echo $@ | cut -d " " -f 2-)
 firstdir=$2
fi

if ! [[ "$1" =~ ^("lin"|"log")$ ]]; then
 echo "USAGE: rhb_log2pdf_all.sh lin/log [dir1/] [dir2/] [dir3/] [...]"
 exit 1
fi

SCRIPT_DIR="$(dirname "$0")"

for i in $(/bin/ls ${firstdir}pop*.log)
do
 i=$(basename $i)
 echo "plotting $i"
 $LOG_CMD $i $1 $sdirs
done
echo "plotting field.log"
$LOG_CMD field.log $1 $sdirs