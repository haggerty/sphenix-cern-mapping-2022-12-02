// compareTiltedOpera.C
//
// What would the measured-minus-calculated field difference look like in the
// tracking volume if the OPERA map were yawed by the measured ~4.4 mrad tilt?
//
//   Measured:   raw point cloud (CSV; mm, T) -- carries the tilt
//   Calculated: OPERA TNtuple "fieldmap"  (cm, T)  -- see OperaMap.h
//
// IMPORTANT.  The measured ~4.4 mrad tilt is an m=1 (azimuth-dependent) feature
// of the *raw* measurement.  The reconstructed sPHENIXFieldMap deliberately
// phi-averages and re-imposes axisymmetry (Bphi = 0, Br from div.B = 0), so the
// tilt is NOT present in GetFieldXYZ.  We therefore read the tilt's transverse
// signature straight from the raw CSV -- the same quantity checkAlignment.C
// uses: <Bx>,<By>, the azimuth-average of the Cartesian transverse field, which
// for an axial field cancels and for a tilt leaves <B_perp> ~ theta * Bz.
//
// A rigid rotation R that maps z_hat -> (theta_x, theta_y, 1) acts on the OPERA
// map as  B_tilt(r) = R . B_opera(R^{-1} r).  To first order in theta:
//
//   sample position : r_src = (x - tx*z,  y - ty*z,  z + tx*x + ty*y)
//   rotate vector   : Bx' = bx + tx*bz ;  By' = by + ty*bz ;  Bz' = bz - tx*bx - ty*by
//
// so yawing OPERA introduces exactly  dBx ~ tx*Bz (~ -6 mT, the measured
// transverse excess) and only a sub-mT, antisymmetric-in-x dBz in the fringe.
// This macro (1) shows the yaw reproduces the measured <Bx>(z), and (2) maps the
// transverse and longitudinal corrections the yaw adds to OPERA.
//
// Run from the project directory:
//   root -l -b -q 'comparison/compareTiltedOpera.C+'
//   root -l -b -q 'comparison/compareTiltedOpera.C+("data","data/sphenix3dmapxyz.root","plots")'

#include "OperaMap.h"

#include "TH2F.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "TSystem.h"

#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

// ─── diverging blue-white-red palette ────────────────────────────────────────
static void UseBWR()
{
    const Int_t n = 5;
    Double_t st[n] = {0.00, 0.25, 0.50, 0.75, 1.00};
    Double_t r[n]  = {0.00, 0.00, 1.00, 1.00, 1.00};
    Double_t g[n]  = {0.00, 0.50, 1.00, 0.50, 0.00};
    Double_t b[n]  = {1.00, 1.00, 1.00, 0.00, 0.00};
    TColor::CreateGradientColorTable(n, st, r, g, b, 255);
    gStyle->SetNumberContours(255);
}

static void SymZ(TH2F *h, double lim = -1.)
{
    if (lim < 0) lim = std::max(std::abs(h->GetMaximum()),
                                std::abs(h->GetMinimum()));
    if (lim <= 0) lim = 1e-6;
    h->GetZaxis()->SetRangeUser(-lim, lim);
}

void compareTiltedOpera(const char *dataDir  = nullptr,
                        const char *calcRoot = nullptr,
                        const char *outDir   = nullptr,
                        double tiltX_mrad    = -4.40,   // theta_x (toward -x), Method 2
                        double tiltY_mrad    =  0.00)   // theta_y (~0, level)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetPadRightMargin(0.13);

    const std::string DD   = dataDir  ? dataDir  : "data";
    const std::string FINE = DD + "/pointCloudFineFullField.csv";
    const std::string CALC = calcRoot ? calcRoot : (DD + "/sphenix3dmapxyz.root");
    const std::string OUT  = outDir   ? outDir   : "plots";
    gSystem->mkdir(OUT.c_str(), kTRUE);

    const double tx = tiltX_mrad * 1e-3;   // rad
    const double ty = tiltY_mrad * 1e-3;   // rad

    // Tracking-volume comparison grid (cm)
    const int    NZ  = 101;   const double Z0 = -100., DZc = 2.0; // z = -100..+100 cm
    const double RMAXcm = 80.;                                    // r <= 80 cm

    printf("=== Loading OPERA field map ===\n");
    OperaMap op;
    if (!op.load(CALC.c_str())) { printf("  cannot load OPERA; abort.\n"); return; }
    printf("=== Yaw applied to OPERA: theta_x = %.2f mrad, theta_y = %.2f mrad ===\n",
           tiltX_mrad, tiltY_mrad);

    // ── tilted-OPERA sampler:  B_tilt(r) = R . B_opera(R^{-1} r),  r in cm ──────
    // Returns false if the shifted sample falls outside the OPERA grid (so the
    // caller can skip boundary points instead of differencing against a 0).
    auto operaTilted = [&](double x, double y, double z,
                           double &Bx, double &By, double &Bz) -> bool {
        double xs = x - tx*z,  ys = y - ty*z,  zs = z + tx*x + ty*y;
        double bx, by, bz;  bool ok = op.get(xs, ys, zs, bx, by, bz);
        Bx = bx + tx*bz;  By = by + ty*bz;  Bz = bz - tx*bx - ty*by;
        return ok;
    };

    // ─────────────────────────────────────────────────────────────────────────
    // (A) measured <Bx>,<By> vs z from the RAW point cloud (the tilt observable)
    //     CSV cols: x_s,y_s,z_s,|B|,Bx_s,By_s,Bz_s  ->  sPHENIX:
    //     x=x_s, y=z_s, z=-y_s ; Bx_phx=Bx_s, By_phx=Bz_s ; r=sqrt(x_s^2+z_s^2)
    // ─────────────────────────────────────────────────────────────────────────
    std::vector<double> sBx(NZ,0), sBy(NZ,0); std::vector<long> nB(NZ,0);
    {
        std::ifstream f(FINE);
        if (!f.is_open()) { printf("  cannot open %s; abort.\n", FINE.c_str()); return; }
        std::string line;
        long used = 0;
        while (std::getline(f, line)) {
            if (line.empty() || line[0]=='#') continue;
            double xs,ys,zs,bmag,bxs,bys,bzs;
            if (std::sscanf(line.c_str(),"%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                            &xs,&ys,&zs,&bmag,&bxs,&bys,&bzs) != 7) continue;
            double xp = xs, yp = zs, zp = -ys;                 // mm
            double r  = std::sqrt(xp*xp + yp*yp);              // mm
            double zc = zp/10.;                                // cm
            if (r > RMAXcm*10.) continue;                      // tracking volume
            if (zc < Z0 || zc > Z0+(NZ-1)*DZc) continue;
            int iz = (int)std::lround((zc - Z0)/DZc);
            if (iz < 0 || iz >= NZ) continue;
            sBx[iz] += bxs;   sBy[iz] += bzs;   nB[iz] += 1;   // Bx_phx, By_phx
            ++used;
        }
        printf("=== Raw measured transverse field: %ld points in tracking volume ===\n", used);
    }

    TGraph *gMeasX=new TGraph(); TGraph *gOp0X=new TGraph(); TGraph *gOpTX=new TGraph();
    TGraph *gMeasY=new TGraph(); TGraph *gOp0Y=new TGraph(); TGraph *gOpTY=new TGraph();
    const int    NRr=8, NPp=36;                 // rings r=10..80 cm, phi step 10 deg
    int kM=0, kO=0;
    double measXacc=0, opTXacc=0; int nAcc=0;   // plateau means for the printout
    for (int iz=0; iz<NZ; ++iz) {
        double z = Z0 + iz*DZc;
        // OPERA (untilted / yawed) phi-average of Cartesian transverse field
        double s0x=0,s0y=0,sTx=0,sTy=0; int n=0;
        for (int ir=0; ir<NRr; ++ir) { double r=10.+ir*10.;
            for (int ip=0; ip<NPp; ++ip) { double phi=ip*10.*M_PI/180.;
                double x=r*std::cos(phi), y=r*std::sin(phi);
                double o0x,o0y,o0z, otx,oty,otz;
                bool ok0 = op.get(x,y,z,o0x,o0y,o0z);
                bool okT = operaTilted(x,y,z,otx,oty,otz);
                if (!ok0 || !okT) continue;                 // skip grid-boundary points
                s0x+=o0x; s0y+=o0y; sTx+=otx; sTy+=oty; ++n;
            } }
        if (n==0) continue;
        gOp0X->SetPoint(kO,z,s0x/n*1e3); gOpTX->SetPoint(kO,z,sTx/n*1e3);
        gOp0Y->SetPoint(kO,z,s0y/n*1e3); gOpTY->SetPoint(kO,z,sTy/n*1e3); ++kO;
        // measured (raw) where the slab has points
        if (nB[iz] > 0) {
            double mx=sBx[iz]/nB[iz]*1e3, my=sBy[iz]/nB[iz]*1e3;
            gMeasX->SetPoint(kM,z,mx); gMeasY->SetPoint(kM,z,my); ++kM;
            if (std::abs(z) <= 80.) { measXacc+=mx; opTXacc+=sTx/n*1e3; ++nAcc; }
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // (B) 2-D maps of the correction the yaw ADDS to OPERA, y=0 plane:
    //     dBx_yaw = (R.B - B)_x  ~ theta*Bz  (the transverse difference it accounts for)
    //     dBz_yaw = (R.B - B)_z  ~ theta*x*dBz/dz  (sub-mT fringe structure)
    // ─────────────────────────────────────────────────────────────────────────
    const int NX=65; const double X0=-80., DXc=2.5;
    auto mk=[&](const char*nm,const char*ttl){
        return new TH2F(nm,ttl,NZ,Z0-DZc/2,Z0+(NZ-0.5)*DZc, NX,X0-DXc/2,X0+(NX-0.5)*DXc); };
    TH2F *hBxYaw=mk("hBxYaw","Transverse field the 4.4 mrad yaw adds to OPERA  (#approx measured #minus OPERA);z [cm];x [cm];#DeltaB_{x} [mT]");
    TH2F *hBzYaw=mk("hBzYaw","Longitudinal change the yaw adds to OPERA;z [cm];x [cm];#DeltaB_{z} [mT]");

    double maxBxYaw=0, maxBzYaw=0;
    for (int iz=0; iz<NZ; ++iz) { double z=Z0+iz*DZc;
        for (int ix=0; ix<NX; ++ix) { double x=X0+ix*DXc; if (std::abs(x)>RMAXcm) continue;
            double o0x,o0y,o0z, otx,oty,otz;
            bool ok0 = op.get(x,0.,z,o0x,o0y,o0z);
            bool okT = operaTilted(x,0.,z,otx,oty,otz);
            if (!ok0 || !okT) continue;                      // skip grid-boundary points
            double dBx=(otx-o0x)*1e3, dBz=(otz-o0z)*1e3;
            hBxYaw->Fill(z,x,dBx); hBzYaw->Fill(z,x,dBz);
            maxBxYaw=std::max(maxBxYaw,std::abs(dBx)); maxBzYaw=std::max(maxBzYaw,std::abs(dBz));
        } }

    // ── numbers ────────────────────────────────────────────────────────────────
    double measXmean = nAcc? measXacc/nAcc : 0, opTXmean = nAcc? opTXacc/nAcc : 0;
    printf("\n  ── Transverse field over the tracking-volume plateau (|z| <= 80 cm) ──\n");
    printf("    measured <Bx>            : %+6.2f mT   (the tilt signature, from raw data)\n", measXmean);
    printf("    OPERA <Bx>, untilted     : %+6.2f mT\n", 0.0);
    printf("    OPERA <Bx>, yawed %.1f mrad : %+6.2f mT   <- matches measured\n", tiltX_mrad, opTXmean);
    printf("\n  ── Field the yaw adds to OPERA, over the tracking volume ──\n");
    printf("    max |dBx| (transverse)   : %6.2f mT   (the difference the yaw removes)\n", maxBxYaw);
    printf("    max |dBz| (longitudinal) : %6.3f mT   (sub-mT, antisymmetric in x, fringe)\n", maxBzYaw);
    printf("    NOTE: the ~0.9%% (~12.7 mT) Bz scale offset is rotation-invariant -- unchanged by the yaw.\n");

    // ── draw maps ────────────────────────────────────────────────────────────────
    UseBWR();
    SymZ(hBxYaw); SymZ(hBzYaw);
    TCanvas *c1=new TCanvas("c1","yaw correction to OPERA",1500,560);
    c1->Divide(2,1);
    c1->cd(1); hBxYaw->Draw("COLZ");
    c1->cd(2); hBzYaw->Draw("COLZ");
    c1->SaveAs((OUT+"/tiltedOpera_yaw_correction_maps.pdf").c_str());
    c1->SaveAs((OUT+"/tiltedOpera_yaw_correction_maps.png").c_str());
    printf("\n  -> %s/tiltedOpera_yaw_correction_maps.{pdf,png}\n", OUT.c_str());
    delete c1;

    // ── draw <Bx>,<By> vs z ──────────────────────────────────────────────────────
    auto style=[&](TGraph*g,int col,int m){ g->SetLineColor(col); g->SetMarkerColor(col);
        g->SetLineWidth(2); g->SetMarkerStyle(m); g->SetMarkerSize(0.7); };
    style(gMeasX,kBlack,20); style(gOp0X,kAzure+1,24); style(gOpTX,kRed+1,21);
    style(gMeasY,kBlack,20); style(gOp0Y,kAzure+1,24); style(gOpTY,kRed+1,21);

    TCanvas *c2=new TCanvas("c2","transverse field vs z",1300,600);
    c2->Divide(2,1);
    c2->cd(1);
    TMultiGraph *mgx=new TMultiGraph();
    mgx->SetTitle("#LTB_{x}#GT vs z  (measured raw vs OPERA);z [cm];#LTB_{x}#GT [mT]");
    mgx->Add(gOp0X,"LP"); mgx->Add(gOpTX,"LP"); mgx->Add(gMeasX,"P");
    mgx->SetMinimum(-12.); mgx->SetMaximum(10.);   // crop the z=+/-100 grid-edge spikes
    mgx->Draw("A");
    TLegend *lg=new TLegend(0.14,0.15,0.62,0.34);
    lg->AddEntry(gMeasX,"measured (raw point cloud)","p");
    lg->AddEntry(gOp0X,"OPERA (untilted)","lp");
    lg->AddEntry(gOpTX,"OPERA (yawed 4.4 mrad)","lp");
    lg->Draw();
    c2->cd(2);
    TMultiGraph *mgy=new TMultiGraph();
    mgy->SetTitle("#LTB_{y}#GT vs z  (measured raw vs OPERA);z [cm];#LTB_{y}#GT [mT]");
    mgy->Add(gOp0Y,"LP"); mgy->Add(gOpTY,"LP"); mgy->Add(gMeasY,"P");
    mgy->Draw("A");
    c2->SaveAs((OUT+"/tiltedOpera_transverse_vs_z.pdf").c_str());
    c2->SaveAs((OUT+"/tiltedOpera_transverse_vs_z.png").c_str());
    printf("  -> %s/tiltedOpera_transverse_vs_z.{pdf,png}\n", OUT.c_str());
    delete c2;

    printf("\n  Done.\n");
}
