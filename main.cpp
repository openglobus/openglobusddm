#include <iostream>

#ifdef _WIN32
#include <direct.h>
#define mkdir(dir, mode) _mkdir(dir)
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include "clipp.h"
#include "HgtFormat.h"
#include "HgtFilesGrid.h"

using namespace clipp;

const double PI = 3.141592653589793238462;
const double POLE = 20037508.34;

vec_t getHeightFromGrid(HgtFilesGrid& grid, const HgtFormat& format, int i, int j) {
    int demFileIndex_i = (int)floor((double)i / (double)(format.nrows - 1));
    int demFileIndex_j = (int)floor((double)j / (double)(format.ncols - 1));
    int i00 = i + demFileIndex_i - format.nrows * demFileIndex_i;
    int j00 = j + demFileIndex_j - format.ncols * demFileIndex_j;
    vec_t h = grid.GetHeight(demFileIndex_i, demFileIndex_j, i00, j00);
    return h;
}

double Merc2Lon(double x) {
    return 180.0 * x / POLE;
}

double Merc2Lat(double y) {
    return 180.0 / PI * (2 * atan(exp((y / POLE) * PI)) - PI / 2);
}

double Lon2Merc(double lon) {
    return lon * POLE / 180.0;
}

double Lat2Merc(double lat) {
    return log(tan((90.0 + lat) * PI / 360.0)) / PI * POLE;
}


double DegTail(double deg) {
	if (deg >= 0) {
		return deg - floor(deg);
	}
	return (-1)*floor(deg) + deg;
}

int main(int argc, char* argv[]){
	std::string outputdir;
	std::string inputdir;
	int zoom = 11;
	int quadSize = 33;
	float topleftlon = -180;
	float topleftlat = 85;
	float bottomrightlon = 180;
	float bottomrightlat = 0;
	std::vector<float> bottomrightlonlat{180, 0};
	vec_t incLat_merc, incLon_merc;

	auto cli = (
		required("-o", "--outputdir").doc("output directory root") \
		& value("outputdir", outputdir),
		required("-i", "--inputdir").doc("input directory containing hgt files from \
			http://www.viewfinderpanoramas.org/Coverage%20map%20viewfinderpanoramas_org3.htm") \
		& value("inputdir", inputdir),
		option("-z", "--zoomlevel").doc("target zoom level") \
		& value("zoom", zoom),
		option("-s", "--start").doc("Start lon lat coordonates") \
		& value("lon", topleftlon) & value("lat", topleftlat),
		option("-e", "--end").doc("End lon lat coordonates") \
		& value("lon", bottomrightlon) & value("lat", bottomrightlat),
		option("-x", "--quadsize").doc("Quadsize resolution default 33") \
		& value("size", quadSize)
		);

	if (!parse(argc, argv, cli)) {
		std::cout << make_man_page(cli, argv[0])\
			.prepend_section("DESCRIPTION", 
                "    This tools generates ddm file for openglobus.org \n"
				"    It needs data from  \n"
				"    - http://www.viewfinderpanoramas.org/Coverage%20map%20viewfinderpanoramas_org3.htm\n"
				"    it is based on the fallowing project: \n"
				"    - https://github.com/openglobus/tools/tree/master/HeightsAdapter"
			    "");
		return 0;
	}

	std::cout << "zoom: " << zoom << "\n";
	std::cout << "topleftlonlat: " << topleftlon << "|" << topleftlat << "\n";
	std::cout << "bottomrightlonlat: " << bottomrightlon << "|" << bottomrightlat << "\n";
	std::cout << "quadsize: " << quadSize << "\n";

	int qn_start = Lon2Merc(topleftlon);
	int qn_end = Lon2Merc(bottomrightlon);
	int qm_start = Lat2Merc(topleftlat);
	int qm_end = Lat2Merc(bottomrightlat);

	int quadsCount = pow(2, zoom); //(dstFieldSize - 1) / 32;
	const int dstFieldSize = quadsCount * (quadSize - 1) + 1;//524289;// = 2^19 + 1 > 1201 * 360 - (360 - 1)
	incLat_merc = incLon_merc = 2.0 * POLE / (vec_t)(dstFieldSize - 1);

	int qm_start_in_grid = (POLE - qm_start) / (incLat_merc * (quadSize - 1));
	int qm_end_in_grid = (POLE - qm_end) / (incLat_merc * (quadSize - 1));
	int qn_start_in_grid = (qn_start - POLE) / (incLon_merc * (quadSize - 1));
	int qn_end_in_grid = (qn_end - POLE) / (incLon_merc * (quadSize - 1));

	std::cout << "qm_start_in_grid: " << qm_start_in_grid << "<->" << qm_end_in_grid << "\n";
	std::cout << "qn_start_in_grid: " << qn_start_in_grid << "<->" << qn_end_in_grid << "\n";

	HgtFilesGrid* demGrid = new HgtFilesGrid();
	
	demGrid->Init(4, inputdir.c_str());
	
	HgtFormat srcHgtFormat(1201, 1201, 1.0 / 1200.0);// 3 / 3600 == 1 / ( 1201 - 1 ) deg.
	HgtFormat srcField(srcHgtFormat.nrows * 180 - (180 - 1), srcHgtFormat.ncols * 360 - (360 - 1));
	HgtFormat dstField(dstFieldSize, dstFieldSize);

	int quadSize2 = quadSize * quadSize;
	vec_t* quadHeightData = new vec_t[quadSize2];

	vec3_t tr[3];
	vec3_t line[2];
	vecSet(line[0], 0.0, 1000000.0, 0.0);
	vecSet(line[1], 0.0, -1000000.0, 0.0);
	//LogAll("Adaptation prepared to proceed...\n");

	double lon_d = 0, lat_d = 0;
	double coordi = 0, coordj = 0;

	for (int qm = 0; qm < quadsCount; qm++) {

		if (qm >= qm_start_in_grid && qm <= qm_end_in_grid) {
			std::cout << qm << "/" << quadsCount << "\n";
			for (int qn = 0; qn < quadsCount; qn++)
			{
				bool isZeroHeight = true;
				//std::cout << "    " << quadsCount << "/" << qn << "\n";
				for (int i = 0; i < quadSize; i++) {
					//std::cout << i << "/" << quadSize << "\n";
					for (int j = 0; j < quadSize; j++) {
						coordi = POLE - ((quadSize - 1) * qm + i) * incLat_merc;
						coordj = (-1) * POLE + ((quadSize - 1) * qn + j) * incLon_merc;

						lat_d = Merc2Lat(coordi);
						lon_d = Merc2Lon(coordj);
						//int test = 100000 % (i+1);
						//std::cout << "    " << test << "\n";
						if (qn == 0 && i == 0 && j == 0 ||
							qn == quadsCount - 1 && i == quadSize - 1 && j == quadSize - 1
							) {
							std::cout << "    lon: " << lon_d << " lat: " << lat_d << "\n";
						}

						int demFileIndex_i = (int)ceil(90.0 - lat_d);
						int demFileIndex_j = (int)floor(180.0 + lon_d);

						vec_t onedlat = DegTail(lat_d);
						vec_t onedlon = DegTail(lon_d);

						int indLat = (int)floor(onedlat / srcHgtFormat.cellsize);
						int i00 = 1200 - 1 - indLat;
						int j00 = (int)floor(onedlon / srcHgtFormat.cellsize);

						vec_t h00 = demGrid->GetHeight(demFileIndex_i, demFileIndex_j, i00, j00);
						vec_t h01 = demGrid->GetHeight(demFileIndex_i, demFileIndex_j, i00, j00 + 1);
						vec_t h10 = demGrid->GetHeight(demFileIndex_i, demFileIndex_j, i00 + 1, j00);

						vec_t cornerLat = 90 - demFileIndex_i;
						vec_t cornerLon = -180 + demFileIndex_j;

						vecSet(tr[0],
							cornerLon + j00 * srcHgtFormat.cellsize,
							h00,
							cornerLat + (indLat + 1) * srcHgtFormat.cellsize);

						vecSet(tr[2],
							cornerLon + (j00 + 1) * srcHgtFormat.cellsize,
							j00 < srcHgtFormat.ncols - 1 ? h01 : h00,
							cornerLat + (indLat + 1) * srcHgtFormat.cellsize);

						vecSet(tr[1],
							cornerLon + j00 * srcHgtFormat.cellsize,
							i00 < srcHgtFormat.nrows - 1 ? h10 : h00,
							cornerLat + indLat * srcHgtFormat.cellsize);

						line[0][X] = line[1][X] = lon_d;
						line[0][Z] = line[1][Z] = lat_d;

						vec_t h11 = 0;

						vec_t edge = ((line[0][X] - tr[2][X]) * (tr[1][Z] - tr[2][Z]) - (line[0][Z] - tr[2][Z]) * (tr[1][X] - tr[2][X]));

						if (edge < 0.0)
						{
							h11 = demGrid->GetHeight(demFileIndex_i, demFileIndex_j, i00 + 1, j00 + 1);

							vecSet(tr[0],
								cornerLon + (j00 + 1) * srcHgtFormat.cellsize,
								i00 < srcHgtFormat.nrows - 1 && j00 < srcHgtFormat.ncols - 1 ? h11 :
								i00 < srcHgtFormat.nrows - 1 ? h10 :
								j00 < srcHgtFormat.nrows - 1 ? h01 :
								h00,
								cornerLat + indLat * srcHgtFormat.cellsize);
						}

						quadHeightData[i * quadSize + j] = 0;

						if (h00 != 0 || h01 != 0 || h10 != 0 || h11 != 0) {
							vec_t h = LineIntersectPlane(tr, line);
							quadHeightData[i * quadSize + j] = h;
							if (h > 0)
								isZeroHeight = false;
						}
					}
				}


				if (!isZeroHeight) {
					//std::cout << "    isZeroHeight : " << isZeroHeight << "\n";
					//std::cout << qm << "\n";
					char ccn[10];
					FILE* fp;

					std::string zoomDir(outputdir);
					sprintf(ccn, "%d", zoom);
					zoomDir.append(ccn);

					mkdir(zoomDir.c_str(), S_IRWXU);

					sprintf(ccn, "%d", qm);
					mkdir(zoomDir.append("/").append(ccn).c_str(), S_IRWXU);

					std::string fileName(zoomDir);
					fileName.append("/");
					sprintf(ccn, "%d", qn);
					fileName.append(ccn).append(".ddm");

					std::cout << "droping: " << fileName << "\n";
					fp = fopen(fileName.c_str(), "wb");

					if (fp == NULL) {
						//LogAll(std::string("Error: ").append(fileName).append("\n").c_str());
						return 1;
					}
					else {
						//LogAll(fileName.append("\n").c_str());
					}

					float* quadHeightData_fl = new float[quadSize2];
					for (int i = 0; i < quadSize2; i++) {
						quadHeightData_fl[i] = quadHeightData[i];
					}

					fwrite(quadHeightData_fl, sizeof(float), quadSize2, fp);
					delete[] quadHeightData_fl;
					fclose(fp);
				}

			}
		}
	}

	delete[] quadHeightData;
    return 0;    
}

