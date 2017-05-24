#!/bin/bash

. docker-tag.sh

docker build -t jakubniwa/chipvis:$TAG .
