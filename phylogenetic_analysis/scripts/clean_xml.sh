#!/bin/bash

for var in "$@"
do 
	xml_name=$(echo $var | sed 's/\.xml//g')
	cat $var | sed 's/totalcount="4" value=".*"/totalcount="4" value=""/g' > ${xml_name}_clean.xml
done
