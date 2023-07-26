# place in main corsair folder for rhybrid compilation
folder=$(pwd)
logfile=$folder/compilation_$(date "+%Y_%m_%d_%H_%M_%S").txt
date -R >$logfile 2>&1
echo "host: " $(hostname) >>$logfile 2>&1
echo "user: " $(whoami) >>$logfile 2>&1
echo "log : " $logfile >>$logfile 2>&1
echo "pwd : " $(pwd) >>$logfile 2>&1

echo " " >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "module list" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
module list >>$logfile 2>&1

cd ../rhybrid/
echo "pwd : " $(pwd) >>$logfile 2>&1

echo " " >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "git show (rhybrid)" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
git show -q >>$logfile 2>&1

echo " " >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "git status (rhybrid)" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
git status >>$logfile 2>&1

echo " " >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "git diff (rhybrid)" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
git diff >>$logfile 2>&1

cd ../corsair/
echo "pwd : " $(pwd) >>$logfile 2>&1

echo " " >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "git show (corsair)" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
git show -q >>$logfile 2>&1

echo " " >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "git status (corsair)" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
git status >>$logfile 2>&1

echo " " >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "git diff (corsair)" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
git diff >>$logfile 2>&1

echo " " >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "make clean" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
make clean >>$logfile 2>&1

echo " " >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "make -j 5" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
echo "========================================================================" >>$logfile 2>&1
make -j 5 >>$logfile 2>&1

echo " " >>$logfile 2>&1
date -R >>$logfile 2>&1

