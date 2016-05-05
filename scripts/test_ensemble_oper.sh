#!/bin/bash

djsons='../jsons'
pjsons='../jsons/proxies'
opath='../tmp'
season='DJF'

./ensemble_oper.py --djsons ${djsons} --pjsons ${pjsons} --opath ${opath}  --season ${season}
