#!/bin/sh

compiler=$1
shift

if [ -z "$compiler" ]; then
	exit 1
fi

tmpfl=$(mktemp /tmp/testcxx.XXXXXX)

if [ -z "$tmpfl" ]; then
	exit 1
fi

echo 'int main(int, char*[]) { return 0; }' > $tmpfl
echo >> $tmpfl

mv $tmpfl $tmpfl.cpp
if $compiler -Werror -o $tmpfl $@ $tmpfl.cpp > $tmpfl.out 2>&1; then
	if test $(grep -i 'error\|warning' $tmpfl.out | wc -l) -eq 0; then
		echo yes
	else
		echo no
	fi
else
	echo no
fi

rm -f $tmpfl $tmpfl.cpp $tmpfl.out
exit 0

