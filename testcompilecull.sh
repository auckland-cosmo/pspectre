#!/bin/sh

compiler=$1
shift

if [ -z "$compiler" ]; then
	exit 1
fi

flags=
for flag in $@; do
	if [ "$(sh ./testcompile.sh $compiler $flags $flag)" = "yes" ]; then
		flags="$flags $flag"
	fi
done

echo $flags
exit 0

