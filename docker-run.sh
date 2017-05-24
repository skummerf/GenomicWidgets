#!/bin/bash

. docker-tag.sh

docker run -d -p 8787:8787 -v `pwd`:/src/chipvis jakubniwa/chipvis:$TAG
