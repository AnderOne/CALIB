#!/bin/bash

for run in $(find $(dirname "$0")/../bin/test/ -executable -type f); do $run --result_code=0 --report_level=short --log_level=message --color_output=1; done
