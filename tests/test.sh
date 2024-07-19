export PATH=.:$PATH
resfile=$(mktemp -p/tmp tctestres.XXXXXXXXX)
bin/tcenum tests/$gr 0 0 $resfile >/dev/null 2>&1

echo 1..1
resd=$(diff ${resfile} tests/$gr.res)
if test "x$resd" = "x"; then 
   echo 'ok 1'
else
   echo 'not ok 1 - output differs'
   echo $resd
fi
