//
//  main.cpp
//  SphericalHarmonics
//  See https://en.wikipedia.org/wiki/Spherical_harmonics
//  Created by larryhou on 2023/10/30.
//

#include <fbxsdk.h>
#include <fbxsdk/core/fbxdatatypes.h>

#include <stdio.h>
#include <stdint.h>
#include <iostream>
#include <vector>
#include <map>

#define epsilon 1e-14

typedef double (*function)(double,double);
using band = std::vector<function>;

struct triangle { int v[3]; };

struct vectex {
    double x,y,z,r;
    
    bool operator == (const vectex &v) const {
        return
        abs(x - v.x) < epsilon &&
        abs(y - v.y) < epsilon &&
        abs(z - v.z) < epsilon && r == v.r;
    }
};

//MARK: 00
double L0M0(double t, double p) { return 1.0/sqrt(M_PI)/2; }

//MARK: 01
double L1P0(double t, double p) { return sqrt(3/M_PI)/2*cos(t); }
double L1P1(double t, double p) { return sqrt(3/M_PI)/2*sin(t)*cos(p); }
double L1M1(double t, double p) { return sqrt(3/M_PI)/2*sin(t)*sin(p); }

//MARK: 02
double L2P0(double t, double p) { return sqrt(5/M_PI)/8*(1 + 3*cos(2*t)); }
double L2P1(double t, double p) { return sqrt(15/M_PI/2)/2*cos(t)*sin(t)*cos(p); }
double L2M1(double t, double p) { return sqrt(15/M_PI/2)/2*cos(t)*sin(t)*sin(p); }
double L2P2(double t, double p) { return sqrt(15/M_PI/2)/4*pow(sin(t),2)*cos(2*p); }
double L2M2(double t, double p) { return sqrt(15/M_PI/2)/4*pow(sin(t),2)*sin(2*p); }

//MARK: 03
double L3P0(double t, double p) {return sqrt(7/M_PI)/16*(3*cos(t)+5*cos(3*t));}
double L3P1(double t, double p) {return sqrt(21/M_PI)/8*(5*pow(cos(t),2)-1)*sin(t)*cos(p);}
double L3M1(double t, double p) {return sqrt(21/M_PI)/8*(5*pow(cos(t),2)-1)*sin(t)*sin(p);}
double L3P2(double t, double p) {return sqrt(105/M_PI/2)/4*cos(t)*pow(sin(t),2)*cos(2*p);}
double L3M2(double t, double p) {return sqrt(105/M_PI/2)/4*cos(t)*pow(sin(t),2)*sin(2*p);}
double L3P3(double t, double p) {return sqrt(35/M_PI)/8*pow(sin(t),3)*cos(3*p);}
double L3M3(double t, double p) {return sqrt(35/M_PI)/8*pow(sin(t),3)*sin(3*p);}

//MARK: 04
double L4P0(double t, double p) {return 3*sqrt(1/M_PI)/128*(9+20*cos(2*t)+35*cos(4*t));}
double L4P1(double t, double p) {return 3*sqrt(5/M_PI)/8*cos(t)*(7*pow(cos(t),2)-3)*sin(t)*cos(p);}
double L4M1(double t, double p) {return 3*sqrt(5/M_PI)/8*cos(t)*(7*pow(cos(t),2)-3)*sin(t)*sin(p);}
double L4P2(double t, double p) {return 3*sqrt(5/M_PI/2)/16*(7*cos(2*t)+5)*pow(sin(t),2)*cos(2*p);}
double L4M2(double t, double p) {return 3*sqrt(5/M_PI/2)/16*(7*cos(2*t)+5)*pow(sin(t),2)*sin(2*p);}
double L4P3(double t, double p) {return 3*sqrt(35/M_PI)/8*cos(t)*pow(sin(t),3)*cos(3*p);}
double L4M3(double t, double p) {return 3*sqrt(35/M_PI)/8*cos(t)*pow(sin(t),3)*sin(3*p);}
double L4P4(double t, double p) {return 3*sqrt(35/M_PI/2)/16*pow(sin(t),4)*cos(4*p);}
double L4M4(double t, double p) {return 3*sqrt(35/M_PI/2)/16*pow(sin(t),4)*sin(4*p);}

//MARK: 05
double L5P0(double t, double p) {return sqrt(11/M_PI)/256*(30*cos(t)+35*cos(3*t)+63*cos(5*t));}
double L5P1(double t, double p) {return sqrt(165/M_PI/2)/128*(15+28*cos(2*t)+21*cos(4*t))*sin(t)*cos(p);}
double L5M1(double t, double p) {return sqrt(165/M_PI/2)/128*(15+28*cos(2*t)+21*cos(4*t))*sin(t)*sin(p);}
double L5P2(double t, double p) {return sqrt(1155/M_PI/2)/16*cos(t)*(1+3*cos(2*t))*pow(sin(t),2)*cos(2*p);}
double L5M2(double t, double p) {return sqrt(1155/M_PI/2)/16*cos(t)*(1+3*cos(2*t))*pow(sin(t),2)*sin(2*p);}
double L5P3(double t, double p) {return sqrt(385/M_PI)/64*(7+9*cos(2*t))*pow(sin(t),3)*cos(3*p);}
double L5M3(double t, double p) {return sqrt(385/M_PI)/64*(7+9*cos(2*t))*pow(sin(t),3)*sin(3*p);}
double L5P4(double t, double p) {return sqrt(385/M_PI/2)*3/16*cos(t)*pow(sin(t),4)*cos(4*p);}
double L5M4(double t, double p) {return sqrt(385/M_PI/2)*3/16*cos(t)*pow(sin(t),4)*sin(4*p);}
double L5P5(double t, double p) {return sqrt(77/M_PI)*3/32*pow(sin(t),5)*cos(5*p);}
double L5M5(double t, double p) {return sqrt(77/M_PI)*3/32*pow(sin(t),5)*sin(5*p);}

//MARK: 06
double L6P0(double t, double p) {return sqrt(13/M_PI)/32*(-5+105*pow(cos(t),2)-315*pow(cos(t),4)+231*pow(cos(t),6));}
double L6P1(double t, double p) {return sqrt(273/M_PI/2)/128*cos(t)*(19+12*cos(2*t)+33*cos(4*t))*sin(t)*cos(p);}
double L6M1(double t, double p) {return sqrt(273/M_PI/2)/128*cos(t)*(19+12*cos(2*t)+33*cos(4*t))*sin(t)*sin(p);}
double L6P2(double t, double p) {return sqrt(1365/M_PI)/512*(35+60*cos(2*t)+33*cos(4*t))*pow(sin(t),2)*cos(2*p);}
double L6M2(double t, double p) {return sqrt(1365/M_PI)/512*(35+60*cos(2*t)+33*cos(4*t))*pow(sin(t),2)*sin(2*p);}
double L6P3(double t, double p) {return sqrt(1365/M_PI)/64*cos(t)*(5+11*cos(2*t))*pow(sin(t),3)*cos(3*p);}
double L6M3(double t, double p) {return sqrt(1365/M_PI)/64*cos(t)*(5+11*cos(2*t))*pow(sin(t),3)*sin(3*p);}
double L6P4(double t, double p) {return sqrt(91/M_PI/2)*3/64*(9+11*cos(2*t))*pow(sin(t),4)*cos(4*p);}
double L6M4(double t, double p) {return sqrt(91/M_PI/2)*3/64*(9+11*cos(2*t))*pow(sin(t),4)*sin(4*p);}
double L6P5(double t, double p) {return sqrt(1001/M_PI)*3/32*cos(t)*pow(sin(t),5)*cos(5*p);}
double L6M5(double t, double p) {return sqrt(1001/M_PI)*3/32*cos(t)*pow(sin(t),5)*sin(5*p);}
double L6P6(double t, double p) {return sqrt(3003/M_PI)/64*pow(sin(t),6)*cos(6*p);}
double L6M6(double t, double p) {return sqrt(3003/M_PI)/64*pow(sin(t),6)*sin(6*p);}

//MARK: 07
double L7P0(double t, double p) {return sqrt(15/M_PI)/32*cos(t)*(-35+315*pow(cos(t),2)-693*pow(cos(t),4)+429*pow(cos(t),6));}
double L7P1(double t, double p) {return sqrt(105/M_PI/2)/2048*(350+675*cos(2*t)+594*cos(4*t)+429*cos(6*t))*sin(t)*cos(p);}
double L7M1(double t, double p) {return sqrt(105/M_PI/2)/2048*(350+675*cos(2*t)+594*cos(4*t)+429*cos(6*t))*sin(t)*sin(p);}
double L7P2(double t, double p) {return sqrt(35/M_PI)*3/512*cos(t)*(109+132*cos(2*t)+143*cos(4*t))*pow(sin(t),2)*cos(2*p);}
double L7M2(double t, double p) {return sqrt(35/M_PI)*3/512*cos(t)*(109+132*cos(2*t)+143*cos(4*t))*pow(sin(t),2)*sin(2*p);}
double L7P3(double t, double p) {return sqrt(35/M_PI/2)*3/512*(189+308*cos(2*t)+143*cos(4*t))*pow(sin(t),3)*cos(3*p);}
double L7M3(double t, double p) {return sqrt(35/M_PI/2)*3/512*(189+308*cos(2*t)+143*cos(4*t))*pow(sin(t),3)*sin(3*p);}
double L7P4(double t, double p) {return sqrt(385/M_PI/2)*3/64*cos(t)*(7+13*cos(2*t))*pow(sin(t),4)*cos(4*p);}
double L7M4(double t, double p) {return sqrt(385/M_PI/2)*3/64*cos(t)*(7+13*cos(2*t))*pow(sin(t),4)*sin(4*p);}
double L7P5(double t, double p) {return sqrt(385/M_PI/2)*3/128*(11+13*cos(2*t))*pow(sin(t),5)*cos(5*p);}
double L7M5(double t, double p) {return sqrt(385/M_PI/2)*3/128*(11+13*cos(2*t))*pow(sin(t),5)*sin(5*p);}
double L7P6(double t, double p) {return sqrt(5005/M_PI)*3/64*cos(t)*pow(sin(t),6)*cos(6*p);}
double L7M6(double t, double p) {return sqrt(5005/M_PI)*3/64*cos(t)*pow(sin(t),6)*sin(6*p);}
double L7P7(double t, double p) {return sqrt(715/M_PI/2)*3/64*pow(sin(t),7)*cos(7*p);}
double L7M7(double t, double p) {return sqrt(715/M_PI/2)*3/64*pow(sin(t),7)*sin(7*p);}

//MARK: 08
double L8P0(double t, double p) {return sqrt(17/M_PI)/256*(35-1260*pow(cos(t),2)+6930*pow(cos(t),4)-12012*pow(cos(t),6)+6435*pow(cos(t),8));}
double L8P1(double t, double p) {return sqrt(17/M_PI/2)*3/2048*cos(t)*(178+869*cos(2*t)+286*cos(4*t)+715*cos(6*t))*sin(t)*cos(p);}
double L8M1(double t, double p) {return sqrt(17/M_PI/2)*3/2048*cos(t)*(178+869*cos(2*t)+286*cos(4*t)+715*cos(6*t))*sin(t)*sin(p);}
double L8P2(double t, double p) {return sqrt(595/M_PI)*3/4096*(210+385*cos(2*t)+286*cos(4*t)+143*cos(6*t))*pow(sin(t),2)*cos(2*p);}
double L8M2(double t, double p) {return sqrt(595/M_PI)*3/4096*(210+385*cos(2*t)+286*cos(4*t)+143*cos(6*t))*pow(sin(t),2)*sin(2*p);}
double L8P3(double t, double p) {return sqrt(19635/M_PI/2)/512*cos(t)*(37+52*cos(2*t)+39*cos(4*t))*pow(sin(t),3)*cos(3*p);}
double L8M3(double t, double p) {return sqrt(19635/M_PI/2)/512*cos(t)*(37+52*cos(2*t)+39*cos(4*t))*pow(sin(t),3)*sin(3*p);}
double L8P4(double t, double p) {return sqrt(1309/M_PI/2)*3/1024*(99+156*cos(2*t)+65*cos(4*t))*pow(sin(t),4)*cos(4*p);}
double L8M4(double t, double p) {return sqrt(1309/M_PI/2)*3/1024*(99+156*cos(2*t)+65*cos(4*t))*pow(sin(t),4)*sin(4*p);}
double L8P5(double t, double p) {return sqrt(17017/M_PI/2)*3/128*cos(t)*(3+5*cos(2*t))*pow(sin(t),5)*cos(5*p);}
double L8M5(double t, double p) {return sqrt(17017/M_PI/2)*3/128*cos(t)*(3+5*cos(2*t))*pow(sin(t),5)*sin(5*p);}
double L8P6(double t, double p) {return sqrt(7293/M_PI)/256*(13+15*cos(2*t))*pow(sin(t),6)*cos(6*p);}
double L8M6(double t, double p) {return sqrt(7293/M_PI)/256*(13+15*cos(2*t))*pow(sin(t),6)*sin(6*p);}
double L8P7(double t, double p) {return sqrt(12155/M_PI/2)*3/64*cos(t)*pow(sin(t),7)*cos(7*p);}
double L8M7(double t, double p) {return sqrt(12155/M_PI/2)*3/64*cos(t)*pow(sin(t),7)*sin(7*p);}
double L8P8(double t, double p) {return sqrt(12155/M_PI/2)*3/256*pow(sin(t),8)*cos(8*p);}
double L8M8(double t, double p) {return sqrt(12155/M_PI/2)*3/256*pow(sin(t),8)*cos(8*p);}


FbxMesh* generate(fbxsdk::FbxScene *scene, function fn, const char* name, int l);

int main(int argc, const char * argv[]) {
    
    std::vector<band> library {
        {[](double a,double b)->double {return 0; }}
    };
    
    auto manager = FbxManager::Create();
    manager->SetIOSettings(FbxIOSettings::Create(manager, IOSROOT));
    auto scene = FbxScene::Create(manager, "scene");
    const FbxSystemUnit::ConversionOptions options = {
        true, /* mConvertRrsNodes */
        true, /* mConvertLimits */
        true, /* mConvertClusters */
        true, /* mConvertLightIntensity */
        true, /* mConvertPhotometricLProperties */
        true  /* mConvertCameraClipPlanes */
    };
    FbxSystemUnit::m.ConvertScene(scene, options);
    
    std::vector<band> bands = {
        {&L0M0},
        {&L1M1, &L1P0, &L1P1},
        {&L2M2, &L2M1, &L2P0, &L2P1, &L2P2},
        {&L3M3, &L3M2, &L3M1, &L3P0, &L3P1, &L3P2, &L3P3},
        {&L4M4, &L4M3, &L4M2, &L4M1, &L4P0, &L4P1, &L4P2, &L4P3, &L4P4},
        {&L5M5, &L5M4, &L5M3, &L5M2, &L5M1, &L5P0, &L5P1, &L5P2, &L5P3, &L5P4, &L5P5},
        {&L6M6, &L6M5, &L6M4, &L6M3, &L6M2, &L6M1, &L6P0, &L6P1, &L6P2, &L6P3, &L6P4, &L6P5, &L6P6},
        {&L7M7, &L7M6, &L7M5, &L7M4, &L7M3, &L7M2, &L7M1, &L7P0, &L7P1, &L7P2, &L7P3, &L7P4, &L7P5, &L7P6, &L7P7},
        {&L8M8, &L8M7, &L8M6, &L8M5, &L8M4, &L8M3, &L8M2, &L8M1, &L8P0, &L8P1, &L8P2, &L8P3, &L8P4, &L8P5, &L8P6, &L8P7, &L8P8},
    };
    
    char name[5];
    for (auto l = 0; l < bands.size(); l++) {
        name[0] = 'L';
        name[1] = '0' + l;
        for (auto m = -l; m <= l; m++) {
            name[2] = m >= 0 ? 'P' : 'M';
            name[3] = '0' + abs(m);
            
            auto geom = generate(scene, bands[l][m+l], name, l);
            geom->GetNode()->LclTranslation.Set(FbxVectorTemplate3<double>(m, 0.0, l - static_cast<int>(bands.size())/2));
            std::cout << name << std::endl;
        }
    }
    
    auto exporter = FbxExporter::Create(manager, "");
    if (exporter->Initialize("SphericalHarmonics.fbx"))
        if (exporter->Export(scene)) { return 0; }
    return 1;
}

FbxMesh* generate(fbxsdk::FbxScene *scene, function fn, const char* name, int l) {
    auto node = FbxNode::Create(scene, name);
    scene->GetRootNode()->AddChild(node);
    
    auto geom = FbxMesh::Create(scene, node->GetName());
    node->SetNodeAttribute(geom);
    
    static const int MESH_SEGMENT_COUNT = 20;
    auto n = MESH_SEGMENT_COUNT + l * 20;
    auto h = n >> 1;
    
    std::vector<vectex> vertice;
    vertice.reserve((h+1)*(n+1));
    std::vector<triangle> faces;
    
    for (auto j = 0; j <= h; j++) {
        double theta(j * M_PI / h);
        for (auto i = 0; i <= n; i++) {
            double phi(i * 2 * M_PI / n);
            auto sign = 1;
            auto radius = fn(theta, phi);
            if (radius < 0){ sign = -1; }
            radius *= sign;
            
            auto x = radius * sin(theta) * cos(phi);
            auto z = radius * sin(theta) * sin(phi);
            auto y = radius * cos(theta);
            
            if (j > 0) {
                int an = static_cast<int>(vertice.size());
                int a_ = an - 1;
                int bn = an - (n + 1);
                int bf = bn + 1;
                
                if (i == 0) {
                    faces.push_back({an,bn,bf});
                } else if (i == n) {
                    faces.push_back({an,a_,bn});
                } else {
                    faces.push_back({an,a_,bn});
                    faces.push_back({an,bn,bf});
                }
            }
            
            vertice.push_back({x,y,z,sign*radius});
        }
    }
    
    auto m = 0;
    std::map<int,int> remap;
    
    for (auto i = 0; i < vertice.size(); i++) {
        auto c = vertice[i];
        auto f = false;
        for (auto j = i - 1; j >= 0 && i - j <= 2 * n; j--) {
            if (vertice[j] == c) {
                auto iter = remap.find(j);
                if (iter == remap.end()) { break; } else {
                    remap[i] = iter->second;
                    f = true;
                }
            }
        }
        
        if (!f) remap[i] = m++;
    }
    
    for (auto i = 0, j = 0; i < vertice.size(); i++) {
        auto r = remap[i];
        if (r >= j) { vertice[r] = vertice[i]; j++; }
    }
    
    geom->InitControlPoints(m);
    for (auto i = 0; i < m; i++) {
        auto &v = vertice[i];
        geom->SetControlPointAt({v.x, v.y, v.z}, i);
    }
    
    for (auto i = 0; i < faces.size(); i++) {
        auto f = (int *)&faces[i];
        geom->BeginPolygon();
        geom->AddPolygon(remap[f[0]]);
        geom->AddPolygon(remap[f[1]]);
        geom->AddPolygon(remap[f[2]]);
        geom->EndPolygon();
    }
    
    auto color = geom->CreateElementVertexColor();
    color->SetMappingMode(FbxLayerElement::eByControlPoint);
    color->SetReferenceMode(FbxLayerElement::eDirect);
    for (auto i = 0; i < m; i++) {
        FbxColor c(0,0,0,1);
        auto r = vertice[i].r;
        r >= 0 ? (c.mRed = r) : (c.mGreen = -r);
        color->GetDirectArray().Add(c);
    }
    
    geom->GenerateNormals();
    return geom;
}


