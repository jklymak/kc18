
qstat -u jklymak
for dd in `ls -d ../results/^*.bak*/(mh-20) `; do
  echo ${dd}
  grep -E 'advcfl_wvel_max|time_seconds|dynstat_theta_mean|dynstat_uvel_max' ${dd}/input/STDOUT.0000 | tail -n 12
  echo 'done'
done
