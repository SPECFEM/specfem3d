#!/bin/bash
#
# see:
# https://www.gnu.org/software/gettext/manual/html_node/config_002eguess.html
#
#
# updates config.guess from GNU config.git repository
echo
wget -O config.guess 'https://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.guess;hb=HEAD'
echo

# updates config.sub from GNU config.git repository
echo
wget -O config.sub 'https://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.sub;hb=HEAD'
echo

# updates install-sh from GNU automake.git repository
echo
wget -O install-sh 'https://git.savannah.gnu.org/gitweb/?p=automake.git;a=blob_plain;f=lib/install-sh;hb=HEAD'

