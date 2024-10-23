#/bin/bash
do=echo
do=
for i in $*;do
  $do rm -rf $i *-$i.out
done
