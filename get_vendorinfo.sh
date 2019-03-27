#!/bin/bash
##Get cids with similar 3d conformers with compounds in "similar3d_67891773group1.csv"-- first pipeline compounds with cids

while read  CID
do

curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$CID/sids/JSON" >>  similar3d_67891773group1_sid.Json

echo $CID
sleep 3

done <similar3d_67891773group1.csv 
