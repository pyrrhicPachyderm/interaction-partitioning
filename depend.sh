#!/bin/bash
#From Peter Miller's "Recursive Make Considered Harmful".
#Adapted for to work for c or c++, by accepting gcc or g++ as its first argument.
#Also adapted to use $(here), taking that as its second argument.
cc="$1"
shift 1
here="$1"
shift 1
dir="$1"
shift 1

"$cc" -MM -MG "$@" |
	case "$dir" in
		"" | ".")
			sed -e 's@^\(.*\)\.o:@.PRECIOUS \1.d \1.o:@'
			;;
		*)
			sed -e "s@^\(.*\)\.o:@.PRECIOUS $dir/\1.d $dir/\1.o:@"
			;;
	esac |
	case "$here" in
		#Convert all relative paths to use $(here).
		"" | ".")
			sed -e 's@\([[:space:]]\)\([^/[:space:]][^[:space:]]\+\)\.\([^[:space:]]\+\)@\1$(here)/\2\.\3@g'
			;;
		*)
			sed -e "s@\([^/]\)$here@\1\$(here)@g"
			;;
	esac
