#!/bin/bash

for i in `seq 1 13`; do
  cp ../$i.in $i.dat.in
  cp ../$i.in $i-klouz10.dat.in
  cp ../$i.in $i-klouz15.dat.in
  cp ../$i.in $i-klouz20.dat.in
  cp ../$i.in $i-klouz3.dat.in
  cp ../$i.in $i-klouz4.dat.in
  cp ../$i.in $i-klouz5.dat.in
done
