// comparison/OperaMap.h
//
// Loader + trilinear interpolator for an OPERA Cartesian field map
// (TNtuple "fieldmap", branches x,y,z in cm and bx,by,bz in Tesla).  The grid
// extent and step are read from the ntuple, so any uniform Cartesian OPERA file
// works.  Positions are in cm.
#ifndef OPERAMAP_H
#define OPERAMAP_H

#include <TFile.h>
#include <TNtuple.h>
#include <vector>
#include <set>
#include <cmath>
#include <algorithm>
#include <cstdio>

struct OperaMap {
    int    nx=0, ny=0, nz=0;
    double x0=0, y0=0, z0=0, dx=0, dy=0, dz=0; // cm
    std::vector<float> bx, by, bz;
    int idx(int ix,int iy,int iz) const { return (ix*ny+iy)*nz+iz; }

    bool load(const char* file) {
        TFile* f = TFile::Open(file);
        if (!f || f->IsZombie()) { fprintf(stderr,"OperaMap: cannot open %s\n",file); return false; }
        TNtuple* nt = (TNtuple*)f->Get("fieldmap");
        if (!nt) { fprintf(stderr,"OperaMap: no TNtuple 'fieldmap' in %s\n",file); f->Close(); return false; }
        float X,Y,Z,BX,BY,BZ;
        nt->SetBranchAddress("x",&X); nt->SetBranchAddress("y",&Y);
        nt->SetBranchAddress("z",&Z); nt->SetBranchAddress("bx",&BX);
        nt->SetBranchAddress("by",&BY); nt->SetBranchAddress("bz",&BZ);
        std::set<float> sx,sy,sz;
        Long64_t N = nt->GetEntries();
        for (Long64_t i=0;i<N;++i){ nt->GetEntry(i); sx.insert(X); sy.insert(Y); sz.insert(Z);}
        nx=sx.size(); ny=sy.size(); nz=sz.size();
        x0=*sx.begin(); y0=*sy.begin(); z0=*sz.begin();
        dx=(*sx.rbegin()-x0)/(nx-1); dy=(*sy.rbegin()-y0)/(ny-1); dz=(*sz.rbegin()-z0)/(nz-1);
        bx.assign((size_t)nx*ny*nz,0); by=bx; bz=bx;
        for (Long64_t i=0;i<N;++i){
            nt->GetEntry(i);
            int ix=lround((X-x0)/dx), iy=lround((Y-y0)/dy), iz=lround((Z-z0)/dz);
            bx[idx(ix,iy,iz)]=BX; by[idx(ix,iy,iz)]=BY; bz[idx(ix,iy,iz)]=BZ;
        }
        printf("OPERA: %lld pts  grid %dx%dx%d  x[%.1f,%.1f] y[%.1f,%.1f] z[%.1f,%.1f] cm (step %.2f)\n",
               (long long)N,nx,ny,nz,x0,*sx.rbegin(),y0,*sy.rbegin(),z0,*sz.rbegin(),dx);
        f->Close();
        return true;
    }

    // Trilinear interpolation; xc,yc,zc in cm.  Returns false (and zeros) outside.
    bool get(double xc,double yc,double zc,double&BX,double&BY,double&BZ) const {
        double fx=(xc-x0)/dx, fy=(yc-y0)/dy, fz=(zc-z0)/dz;
        if(fx<0||fx>nx-1||fy<0||fy>ny-1||fz<0||fz>nz-1){BX=BY=BZ=0;return false;}
        int ix=std::min((int)fx,nx-2), iy=std::min((int)fy,ny-2), iz=std::min((int)fz,nz-2);
        double tx=fx-ix,ty=fy-iy,tz=fz-iz;
        auto I=[&](const std::vector<float>&G){
            return G[idx(ix,iy,iz)]*(1-tx)*(1-ty)*(1-tz)+G[idx(ix+1,iy,iz)]*tx*(1-ty)*(1-tz)
                  +G[idx(ix,iy+1,iz)]*(1-tx)*ty*(1-tz)+G[idx(ix+1,iy+1,iz)]*tx*ty*(1-tz)
                  +G[idx(ix,iy,iz+1)]*(1-tx)*(1-ty)*tz+G[idx(ix+1,iy,iz+1)]*tx*(1-ty)*tz
                  +G[idx(ix,iy+1,iz+1)]*(1-tx)*ty*tz+G[idx(ix+1,iy+1,iz+1)]*tx*ty*tz;
        };
        BX=I(bx);BY=I(by);BZ=I(bz);return true;
    }
};

#endif // OPERAMAP_H
