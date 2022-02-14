#!/bin/bash

. /kb/deployment/user-env.sh

python ./scripts/prepare_deploy_cfg.py ./deploy.cfg ./work/config.properties

if [ -f ./work/token ] ; then
  export KB_AUTH_TOKEN=$(<./work/token)
fi

if [ $# -eq 0 ] ; then
  sh ./scripts/start_server.sh
elif [ "${1}" = "test" ] ; then
  TESTING=$(grep ^testing ./deploy.cfg | awk '{print $3}')
  if [ "${TESTING}" = "1" ] ; then
    echo "Initializing test data"
    TEST_DATA=$(grep ^test_data ./deploy.cfg | awk '{print $3}')
    mkdir -p /kb/module/work/test_data
    cd /kb/module/work/test_data
    if [ -d ${TEST_DATA} ] && [ -f ${TEST_DATA}/SpeciesIDs.txt ] ; then
      echo "Test data already present, skipping initialization"
    else
      PROTEIN_FAMILY_FILE="${TEST_DATA}.tar.gz"
      if ! [ -f ${PROTEIN_FAMILY_FILE} ] ; then
	echo "Fetching ${PROTEIN_FAMILY_FILE}"
	wget https://web.cels.anl.gov/~seaver/KBase_App_Files/${PROTEIN_FAMILY_FILE}
      fi
      echo "Gunzipping ${PROTEIN_FAMILY_FILE}"
      tar -zxf $PROTEIN_FAMILY_FILE
      if [ -f ${TEST_DATA}/SpeciesIDs.txt ] ; then
	echo "Test data initialized"
      else
	echo "Test data initialization failed"
	exit 1
      fi
    fi
  fi
  echo "Running Tests"
  cd /kb/module
  make test
elif [ "${1}" = "async" ] ; then
  sh ./scripts/run_async.sh
elif [ "${1}" = "init" ] ; then
  echo "Initialize module"
  cd /data
  PROTEIN_FAMILY_FILE="OrthoFinder_Phytozome_Reference.tar.gz"
  wget https://web.cels.anl.gov/~seaver/KBase_App_Files/${PROTEIN_FAMILY_FILE}
  echo "Gunzipping file"
  tar -zxf $PROTEIN_FAMILY_FILE
  if [ -f OrthoFinder_Phytozome_Reference/SpeciesIDs.txt ] ; then
      touch __READY__
      echo "Init succeeded"
  else
      echo "Init failed"
  fi
elif [ "${1}" = "bash" ] ; then
  bash
elif [ "${1}" = "report" ] ; then
  export KB_SDK_COMPILE_REPORT_FILE=./work/compile_report.json
  make compile
else
  echo Unknown
fi
