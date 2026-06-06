// comparison/compareNewVsOpera.C
//
// Quick overlay of the CERN-measured solenoid map (2022-12-02 pointCloud CSVs,
// read through sPHENIXFieldMap) against the official OPERA map
// (sphenixoperamaps/sphenix3dmapxyz.root; TNtuple "fieldmap"; cm, Tesla).
//
//   root -l -b -q 'comparison/compareNewVsOpera.C+("data","data/sphenix3dmapxyz.root","plots")'
//
// The OPERA grid extent/step are read from the ntuple, so this works for any
// uniform Cartesian OPERA file.  For the rigorous Fourier/tilt comparison see
// compareFieldMaps.C.

#include "../sPHENIXFieldMap.cxx"

#include <TFile.h>
#include <TNtuple.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TString.h>
#include <vector>
#include <set>
#include <cmath>
#include <cstdio>

// ── OPERA Cartesian map (units cm, Tesla); grid read from the ntuple ──────────
struct OperaMap {
    int    nx=0, ny=0, nz=0;
    double x0=0, y0=0, z0=0, dx=0, dy=0, dz=0; // cm
    std::vector<float> bx, by, bz;
    int idx(int ix,int iy,int iz) const { return (ix*ny+iy)*nz+iz; }

    void load(const char* file) {
        TFile* f = TFile::Open(file);
        TNtuple* nt = (TNtuple*)f->Get("fieldmap");
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
        printf("OPERA: %lld pts  grid %dx%dx%d  x[%.1f,%.1f] y[%.1f,%.1f] z[%.1f,%.1f] cm (step %.1f)\n",
               (long long)N,nx,ny,nz,x0,*sx.rbegin(),y0,*sy.rbegin(),z0,*sz.rbegin(),dx);
        f->Close();
    }
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

void compareNewVsOpera(const char* datadir   = "data",
                       const char* operaFile = "data/sphenix3dmapxyz.root",
                       const char* outdir    = "plots")
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    sPHENIXFieldMap fmap(Form("%s/pointCloudFineFullField.csv",  datadir),
                         Form("%s/pointCloudRoughFullField.csv", datadir));
    OperaMap op; op.load(operaFile);

    const double zlo=op.z0*10, zhi=(op.z0+(op.nz-1)*op.dz)*10;       // mm
    const double rmax=std::min(800.0, (op.x0+(op.nx-1)*op.dx)*10);   // mm

    TGraph *gMz=new TGraph, *gOz=new TGraph, *gDz=new TGraph;
    for(double z=zlo; z<=zhi; z+=20){
        double Br,Bphi,Bz; fmap.GetField(0,0,z,Br,Bphi,Bz);
        double ox,oy,oz; op.get(0,0,z/10.0,ox,oy,oz);
        gMz->AddPoint(z,Bz); gOz->AddPoint(z,oz); gDz->AddPoint(z,(Bz-oz)*1000);
    }
    TGraph *gMr=new TGraph, *gOr=new TGraph, *gMrr=new TGraph, *gOrr=new TGraph;
    for(double r=0;r<=rmax;r+=10){
        double Br,Bphi,Bz; fmap.GetField(r,0,0,Br,Bphi,Bz);
        double ox,oy,oz; op.get(r/10.0,0,0,ox,oy,oz);
        gMr->AddPoint(r,Bz); gOr->AddPoint(r,oz);
        gMrr->AddPoint(r,Br); gOrr->AddPoint(r,ox);
    }
    TH2D* hD=new TH2D("hD","B_{z} measured - OPERA [mT];z [mm];r [mm]",
                      (int)((zhi-zlo)/20),zlo,zhi,(int)(rmax/25)+1,0,rmax);
    for(int iz=1;iz<=hD->GetNbinsX();++iz)for(int ir=1;ir<=hD->GetNbinsY();++ir){
        double z=hD->GetXaxis()->GetBinCenter(iz), r=hD->GetYaxis()->GetBinCenter(ir);
        double Br,Bphi,Bz; fmap.GetField(r,0,z,Br,Bphi,Bz);
        double ox,oy,oz; if(!op.get(r/10.0,0,z/10.0,ox,oy,oz)) continue;
        hD->SetBinContent(iz,ir,(Bz-oz)*1000);
    }

    auto style=[](TGraph*g,int c){g->SetLineColor(c);g->SetLineWidth(3);g->SetMarkerColor(c);};
    style(gMz,kRed+1);style(gOz,kBlue+1);style(gMr,kRed+1);style(gOr,kBlue+1);
    style(gMrr,kRed+1);style(gOrr,kBlue+1);style(gDz,kBlack);

    TCanvas* c=new TCanvas("c","New measured map vs OPERA",1500,1000);
    c->Divide(2,2);
    c->cd(1);
    TMultiGraph*m1=new TMultiGraph();
    m1->SetTitle("B_{z} on axis (r=0);z [mm];B_{z} [T]");
    m1->Add(gMz,"L");m1->Add(gOz,"L");m1->Draw("A");
    TLegend*l1=new TLegend(0.35,0.18,0.65,0.35);
    l1->AddEntry(gMz,"CERN measured (2022-12-02)","l");
    l1->AddEntry(gOz,"OPERA (official)","l");l1->Draw();
    c->cd(2);
    gDz->SetTitle("B_{z} difference on axis (measured - OPERA);z [mm];#DeltaB_{z} [mT]");
    gDz->Draw("AL");
    TLine*z0=new TLine(zlo,0,zhi,0);z0->SetLineStyle(2);z0->Draw();
    c->cd(3);
    TMultiGraph*m3=new TMultiGraph();
    m3->SetTitle("B_{z} vs r at z=0;r [mm];B_{z} [T]");
    m3->Add(gMr,"L");m3->Add(gOr,"L");m3->Draw("A");
    TLegend*l3=new TLegend(0.18,0.18,0.55,0.35);
    l3->AddEntry(gMr,"CERN measured","l");l3->AddEntry(gOr,"OPERA","l");l3->Draw();
    c->cd(4); hD->Draw("COLZ");
    c->SaveAs(Form("%s/compare_newVsOpera.pdf",outdir));
    c->SaveAs(Form("%s/compare_newVsOpera.png",outdir));

    TCanvas* c2=new TCanvas("c2","Br vs r",700,550);
    TMultiGraph*m2=new TMultiGraph();
    m2->SetTitle("B_{r} vs r at z=0;r [mm];B_{r} [T]");
    m2->Add(gMrr,"L");m2->Add(gOrr,"L");m2->Draw("A");
    TLegend*l2=new TLegend(0.18,0.7,0.55,0.88);
    l2->AddEntry(gMrr,"CERN measured","l");l2->AddEntry(gOrr,"OPERA","l");l2->Draw();
    c2->SaveAs(Form("%s/compare_Br_z0.png",outdir));

    printf("Saved compare_newVsOpera.{pdf,png} and compare_Br_z0.png to %s/\n", outdir);
}
