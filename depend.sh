#!/bin/bash
#From Peter Miller's "Recursive Make Considered Harmful".
#Adapted for to work for c or c++, by accepting gcc or g++ as its first argument.
cc="$1"
shift 1
dir="$1"
shift 1
case "$dir" in
	"" | ".")
		"$cc" -MM -MG "$@" | sed -e 's@^\(.*\)\.o:@.PRECIOUS \1.d \1.o:@'
		;;
	*)
		"$cc" -MM -MG "$@" | sed -e "s@^\(.*\)\.o:@.PRECIOUS $dir/\1.d $dir/\1.o:@"
		;;
esac
