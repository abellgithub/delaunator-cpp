#!/usr/bin/env bash

mkdir -p ./includes
rm -rf ./includes/*

git clone --branch v1.1.0 --depth=1 git@github.com:Tencent/rapidjson.git ./includes/rapidjson