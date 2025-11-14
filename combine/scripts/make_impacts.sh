#!/usr/bin/env sh
#

cardname=$1
filename="${cardname/txt/root}"
if [ -z "$2" ]; then
    outputname=''
else
    outputname="_$2"
fi
bounds=${3:-20}

# blind='-t -1 --expectSignal=1'
# blind='--expectSignal=0'
blind=''

if [ ${cardname} != ${filename} ]; then
    text2workspace.py ${cardname}
fi
combineTool.py -M Impacts -d ${filename} -m 100 --rMin -${bounds} --rMax ${bounds} --doInitialFit --robustFit 1 --parallel 40
combineTool.py -M Impacts -d ${filename} -m 100 --rMin -${bounds} --rMax ${bounds} --robustFit 1 --doFits --parallel 40
combineTool.py -M Impacts -d ${filename} -m 100 -o impacts.json
plotImpacts.py -i impacts.json -o impacts${outputname}
# combineTool.py -M Impacts -d ${filename} ${blind} -m 125 --doInitialFit --robustFit 1 --rMin -${bounds} --rMax ${bounds} --parallel 10 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH
# combineTool.py -M Impacts -d ${filename} ${blind} -m 125 --doFits --robustFit 1 --rMin -${bounds} --rMax ${bounds} --parallel 10 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH
# combineTool.py -M Impacts -d ${filename} ${blind} -m 125 -o impacts.json
# plotImpacts.py -i impacts.json -o impacts${outputname}
