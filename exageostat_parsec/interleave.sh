#!/bin/bash
# find the NUMA nodes that matches the cpu mask
nodes=$(numactl --hardware|grep available|awk '{print $2}')
n=0
while [ $n -lt $nodes ] ; do
  for cpu in $(numactl --hardware|grep "node $n cpus:"|sed "s/\(.*\): //");do
    numa[$cpu]=$n
  done
  nused[$n]=0
  ((n=n+1))
done
mask=$(taskset -cp $$|sed "s/\(.*\): //"|sed "s/^/[/"|sed "s/$/]/")
for i in $(scontrol show hostname $mask);do
  n=${numa[$i]}
  nused[$n]=1
done
#echo ${nused[*]}
interleaving=
sep=
n=0
while [ $n -lt $nodes ];do
  if [ ${nused[$n]} -eq 1 ] ; then
    interleaving="${interleaving}${sep}${n}"
    sep=,
  fi
  ((n=n+1))
done
#echo numactl --interleave=$interleaving "$@"
numactl --interleave=$interleaving "$@"
