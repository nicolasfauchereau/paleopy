#!/bin/bash

./proxy_oper.py --djsons '../jsons' --pjsons '../jsons/proxies' --opath '../tmp' --pfname Crow_Glacier.json --name 'Crow Glacier' --ptype 'Tree rings' --longitude 170.998122 --latitude -43.346222 --dataset 'vcsn' --variable 'TMean' --season 'DJF' --qualitative 0 --value -1 --period '1979-2012' --climatology '1981-2010' --calc_anoms 1 --detrend 1 --aspect 12 --elevation 3500 --dating 'Absolute' --calendar '1000 CE' --chronology 'Dendrochronology, Surface exposure age data (SED)' --measurement 'Ring width' --verbose 1
