#!/bin/bash

echo "Condor job status"

echo -e -n "2016pre QCD jobs:   "; condor_q -af:j Cmd|grep 2016pre|grep -c QCD
echo -e -n "2016post QCD jobs:  "; condor_q -af:j Cmd|grep 2016post|grep -c QCD
echo -e -n "2017 QCD jobs:      "; condor_q -af:j Cmd|grep 2017|grep -c QCD
echo -e -n "2018 QCD jobs:      "; condor_q -af:j Cmd|grep 2018|grep -c QCD
echo "---------------------------------------------------------------------------"
echo -e -n "2016pre DEle jobs:  "; condor_q -af:j Cmd|grep 2016pre|grep -cv QCD
echo -e -n "2016post DEle jobs: "; condor_q -af:j Cmd|grep 2016post|grep -cv QCD
echo -e -n "2017 DEle jobs:     "; condor_q -af:j Cmd|grep 2017|grep -cv QCD
echo -e -n "2018 DEle jobs:     "; condor_q -af:j Cmd|grep 2018|grep -cv QCD

condor_q -totals
condor_release scooper
