//
//  main.cpp
//  SphericalHarmonics
//  See https://en.wikipedia.org/wiki/Spherical_harmonics
//  Created by larryhou on 2023/10/30.
//

#include <fbxsdk.h>
#include <fbxsdk/core/fbxdatatypes.h>

#include <stdint.h>
#include <iostream>
#include <vector>
#include <map>

#define epsilon 1e-14

using kernel = std::function<double(double,double)>;
using band = std::vector<kernel>;

struct triangle { int data[3]; };

struct vectex {
    double x,y,z;
    int sign;
    double raidus;
    
    bool operator == (const vectex &v) const {
        return
        abs(x - v.x) < epsilon &&
        abs(y - v.y) < epsilon &&
        abs(z - v.z) < epsilon && sign == v.sign;
    }
};

double f31(double theta, double phi) {
    auto costheta = cos(theta);
    return sqrt(21/M_PI/2)/4*(5*costheta*costheta-1)*sin(theta)*(cos(phi));
}

FbxMesh* generate(fbxsdk::FbxScene *scene, kernel fn, const char* name);

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
    generate(scene, &f31, "L3M1");
    
    auto exporter = FbxExporter::Create(manager, "");
    if (exporter->Initialize("SphericalHarmonics.fbx"))
        if (exporter->Export(scene)) { return 0; }
    return 1;
}

FbxMesh* generate(fbxsdk::FbxScene *scene, kernel fn, const char* name) {
    auto node = FbxNode::Create(scene, name);
    scene->GetRootNode()->AddChild(node);
    
    auto geom = FbxMesh::Create(scene, node->GetName());
    node->SetNodeAttribute(geom);
    
    auto h = 40;
    auto n  = h * 2;
    
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
            
            vertice.push_back({x,y,z,sign, radius});
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
        FbxColor c;
        if (vertice[i].sign == 1) {
            c.mRed   = vertice[i].raidus;
        } else {
            c.mGreen = vertice[i].raidus;
        }
        color->GetDirectArray().Add(c);
    }
    
    return geom;
}


