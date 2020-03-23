
# Openglobus DDM file adapter

## Introduction
    
This tool generates ddm tile files for https://www.openglobus.org/ .

## To build it:

```shell
git clone https://github.com/screwt/openglobusddm.git
mkdir build
cd build
cmake ../openglobusddm
make
```

## Usage


You will need data source files from:

* http://www.viewfinderpanoramas.org/Coverage%20map%20viewfinderpanoramas_org3.htm

Unzip all files in the same directory

```shell
./HeightAdapter -o /output/ -i /tmp/file_fome_viewfinderpanoramas/ -z 8 -s -5 51 -e 9 42
```
* -z zoom level (exemple: 8)
* -s lonlat starting point (exemple: -5 51)
* -e lonlat ending point (exemple: 9 42)
* -x quadsize (resolution) default: 33

# INFO

This tool is based on this project from the openglobus team:

* https://github.com/openglobus/tools/tree/master/HeightsAdapter

