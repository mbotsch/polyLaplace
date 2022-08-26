#!/bin/bash

# Exit script on any error
set -e 

#=====================================
# Color Settings:
#=====================================
NC='\033[0m'
OUTPUT='\033[0;32m'
WARNING='\033[0;93m'

BUILDPATH="build-cppcheck"
mkdir -p "${BUILDPATH}"
cd "${BUILDPATH}"
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..
cd ..

echo -e "${OUTPUT}"
echo "=============================================================================="
echo "Running cppcheck"
echo -n "Version: "
cppcheck --version
echo "=============================================================================="
echo -e "${NC}"
echo "Please Wait ..."

# Run cppcheck and output into file
# we do not enable 'style' and 'unusedFunction' (the latter gives false positive
# for the public library interface)
cppcheck \
    --project=${BUILDPATH}/compile_commands.json \
    --force \
    --enable=warning,performance,portability,information,missingInclude \
    --suppress=missingIncludeSystem \
    --suppress=unmatchedSuppression \
    --inline-suppr \
    --quiet \
    2>&1 | tee cppcheck.log

COUNT=$(wc -l < cppcheck.log )

echo -e "${OUTPUT}"
echo "=============================================================================="
echo "CPPCHECK Summary"
echo "=============================================================================="
echo -e "${NC}"

if [ $COUNT -gt 3 ]; then
  echo -e ${WARNING}
  echo "Total CPPCHECK error Count is $COUNT, which is too High! CPPCHECK Run failed";
  echo -e "${NC}"
  exit 1;
else
  echo "Total CPPCHECK error Count is $COUNT ... OK"
fi

 
