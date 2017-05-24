#!/bin/bash

. docker-tag.sh

docker run -v `pwd`:/src/chipvis -t jakubniwa/chipvis:$TAG \
  Rscript -e "setwd('/src/chipvis/');\
              suppressMessages(devtools::load_all());\
              devtools::test()"

