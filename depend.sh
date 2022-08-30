#!/bin/bash
#From Peter Miller's "Recursive Make Considered Harmful".
#Adapted for to work for c or c++, by accepting gcc or g++ as its first argument.
CC="$1"
shift 1
DIR="$1"
shift 1
case "$DIR" in
	"" | ".")
		"$CC" -MM -MG "$@" | sed -e 's@^\(.*\)\.o:@.PRECIOUS \1.d \1.o:@'
		;;
	*)
		"$CC" -MM -MG "$@" | sed -e "s@^\(.*\)\.o:@.PRECIOUS $DIR/\1.d $DIR/\1.o:@"
		;;
esac
