#!/bin/sh

# changes to wiki/ directory
if [ -d ../../wiki ]; then
  cd ../../wiki/
else
  echo "wiki directory does not exist, nothing to do..."
  exit 1
fi

# checks if translation file exists
if [ ! -f translate_user_manual_to_markdown.pl ]; then
  echo "translation script is not available, exiting..."
  exit 1
fi

# runs translation script
./translate_user_manual_to_markdown.pl

echo
echo

