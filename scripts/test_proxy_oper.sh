#!/bin/bash

djsons='../jsons'
pjsons='../jsons/proxies'
opath='../tmp'
pfname='Crow_Glacier.json'

./proxy_oper.py --djsons ${djsons} --pjsons ${pjsons} --opath ${opath} --pfname ${pfname} --name 'Crow Glacier' --ptype 'Tree rings' --longitude 170.998122 --latitude -43.346222 --dataset 'vcsn' --variable 'TMean' --season 'DJF' --qualitative 1 --value 'WB' --period '1979-2012' --climatology '1981-2010' --calc_anoms 0 --detrend 0 --aspect 12 --elevation 3500 --dating 'Absolute' --calendar '1000 CE' --chronology 'Dendrochronology, Surface exposure age data (SED)' --measurement 'Ring width' --verbose 1
