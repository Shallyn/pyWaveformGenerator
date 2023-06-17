#! /bin/sh

{
{ echo "running aclocal (please ignore \"underquoted\" warnings)..." ; aclocal; } &&
{ echo "running libtoolize..." ; libtoolize -i -c -f || glibtoolize $LIBTOOLIZE_FLAGS ;} &&
{ echo "running autoheader..."; autoheader ; } &&
{ echo "running automake..." ; automake -a ; } &&
{ echo "running autoconf..." ; autoconf ; } &&
echo "$0 complete." ;
} || { echo "$0 failed." ; false ; }

