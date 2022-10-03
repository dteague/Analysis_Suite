Full Run-Through
================

In this section, we give an example set of code used for the 3top analysis. These are the commands that are used to go from NanoAOD files to a result from combine. This set of instructions tries to be as thorough as possible and will be tested periodically to ensure it should work on all setups.

Currently, the code base is only useable on the wisconsin cluster, but efforts to expand to using on other machines, namely lxplus are on the way. The only thing needed at the Wisconsin cluster is a way to submit jobs and a place to store the output. Wisconsin has a script called "farmoutAnalysis" which submits the jobs to condor and writes the output to a hadoop area. These files are hadd'ed and stored on the hdfs area for running. To genericize the code, this step is what needs to be fixed

Each code section does a distinct part of the analysis needed in whole, but each is done such that one can run each section individually. This is also done for 2018, but easily expanded to 2016pre, 2016post, 2017.

********************
Initial Job Creation
********************

*************************
Make Nonprompt Fake Rates
*************************

.. code-block:: bash
   voms-proxy-init -rfc -valid 144:00 -voms cms -bits 2048
   create_lists.py -y 2018 -t mc -r output -g ttbar,DY,wjets,qcd -f nonprompt_measurement
   create_lists.py -y 2018 -t data -r output -g MM,E -f nonprompt_measurement
   create_lists.py -y 2018 -t mc -r output -g ttbar_lep,wjet_ht -f nonprompt_closure
   create_lists.py -y 2018 -t data -r output -g MM,EM,EE -f nonprompt_closure
   farmoutJob.py -f nonprompt_measurement -y 2018 -a FakeRate
   farmoutJob.py -f nonprompt_closure -y 2018 -a Nonprompt_Closure
   # After jobs are finished (monitor with condor_q)
   hadd.py -f nonprompt_measurement -y 2018 -a FakeRate -d run_through
   hadd.py -f nonprompt_closure -y 2018 -a Nonprompt_Closure -d run_through
   # creates and moves files to directory in /hdfs/store/user/<user>/workspace/charge_misId/2018/run_through
   calculate_fr.py -y 2018 -w run_through -r sideband,measurement,closure -d run_through
   # Run after make both fake rates for charge_misId and Nonprompt
   convert_fakerates.py -y 2018 --nonprompt run_through --charge run_through

****************************************
Make Charge Misidentification Fake Rates
****************************************

.. code-block:: bash
   voms-proxy-init -rfc -valid 144:00 -voms cms -bits 2048
   create_lists.py -y 2018 -t mc -r output -g dy,vv,ttbar_lep -f charge_misId
   create_lists.py -y 2018 -t data -r output -g MM,EM,EE -f charge_misId
   farmoutJob.py -f charge_misId -y 2018 -a Closure_MisId
   # After jobs are finished (monitor with condor_q)
   hadd.py -f charge_misId -y 2018 -a Closure_MisId -d run_through
   # creates and moves files to directory in /hdfs/store/user/<user>/workspace/charge_misId/2018/run_through
   calculate_misId.py -y 2018 -w run_through -r measurement,closure -d run_through
   # Run after make both fake rates for charge_misId and Nonprompt
   convert_fakerates.py -y 2018 --nonprompt run_through --charge run_through

***************************
Make Btagging Efficiencies
**************************
This snippet of code runs the basic btagging on signal like events and then calculates the MC btagging efficiencies for different flavor jets. This is used for getting the correct btagging scale factors. These scale factors are found as ``data/POG/USER/<year>/beff.json.gz``

.. code-block:: bash
   voms-proxy-init -rfc -valid 144:00 -voms cms -bits 2048
   create_lists.py -y 2018 -t mc -r output -g ttt,tttt,ttX,ttXY,xg,vv,vvv,ttbar,dy,wjet,other -f befficiency
   farmoutJob.py -f befficiency -y 2018 -a BEfficiency -t mc
   # After jobs are finished (monitor with condor_q)
   hadd.py -f befficiency -y 2018 -a BEfficiency -t mc -d run_through
   # creates and moves files to directory in /hdfs/store/user/<user>/workspace/befficiency/2018/run_through
   beff.py -y 2018 -w run_through
