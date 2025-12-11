#!/bin/bash -x



# RECO steps for Summer20UL18
cmsDriver.py \
  --python_filename "RunIISummer20UL18RECO_cfg.py" \
  --eventcontent AODSIM \
  --customise Configuration/DataProcessing/Utils.addMonitoring \
  --datatier AODSIM \
  --filein "file:RunIISummer20UL18HLT.root" \
  --fileout "file:RunIISummer20UL18RECO.root" \
  --conditions 106X_upgrade2018_realistic_v11_L1v1 \
  --step RAW2DIGI,L1Reco,RECO,RECOSIM \
  --geometry DB:Extended \
  --era Run2_2018 \
  --runUnscheduled \
  --mc \
  --no_exec \
  --nThreads 1 \
  -n 100 || exit $?




# MiniAODv2 
cmsDriver.py \
  --python_filename "RunIISummer20UL18MINIAODSIM_cfg.py" \
  --eventcontent MINIAODSIM \
  --customise Configuration/DataProcessing/Utils.addMonitoring \
  --datatier MINIAODSIM \
  --filein "file:RunIISummer20UL18RECO.root" \
  --fileout file:mini.root \
  --conditions 106X_upgrade2018_realistic_v16_L1v1 \
  --step PAT \
  --procModifiers run2_miniAOD_UL \
  --geometry DB:Extended \
  --era Run2_2018 \
  --runUnscheduled \
  --no_exec \
  --nThreads 1 \
  --mc \
  -n 1


