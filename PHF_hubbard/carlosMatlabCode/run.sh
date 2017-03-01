#!/bin/bash

HUBOUT='test'
HUBLOAD=''

export HUBOUT
export HUBLOAD

matlab -nodisplay -r input_hfopt > test.log 2>&1 &
wait


exit 0

