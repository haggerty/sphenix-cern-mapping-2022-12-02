// analysis/estimateTilt.C
//
// Two INDEPENDENT estimates of the solenoid tilt, printed side by side, using
// only measured points (no r=0 extrapolation through sPHENIXFieldMap).
//
//   Method 1  "near-axis m=0 transverse field"
//     theta(r) = atan2(<B_transverse>, <Bz>) from the phi-averaged (m=0)
//     transverse field at each measured radius, in the flat central-z band.
//     A rigid tilt predicts theta(r) = const; r-dependence flags a systematic.
//     Reported at the innermost measured ring -> cleanest, least model-dependent.
//
//   Method 2  "magnetic-axis line + tilt fit" (fringe-driven, separates tilt
//     from translation)
//     At each z the transverse field is  Bx = theta_x*B - g*(x - x0 - theta_x*z),
//     g = -1/2 dB/dz.  A per-z linear fit Bx = px(z) + g(z)*x gives px(z), g(z);
//     a global least-squares  px(z) = theta_x*(B(z) - g(z) z) - x0*g(z)  then
//     separates the tilt slope theta_x from the axis offset x0 (same for y).
//     Uses the bore coverage where g != 0, i.e. real interpolation, not r=0.
//
// The two methods use physically distinct signatures (near-axis uniform field
// vs. fringe focusing gradient) and should agree if the tilt is a real rigid
// rotation.
//
// CSV (surveyor frame): x_s,y_s,z_s,|B|,Bx_s,By_s,Bz_s  [mm, T]
//   sPHENIX transform: xp=xs, yp=zs, zp=-ys; Bx=Bxs, By=Bzs, Bz=-Bys
//
// Run from the repository root:
//   root -l -b -q 'analysis/estimateTilt.C+("data","plots")'
//
// TODO(consolidation): fold the m=0 (checkAlignment.C) and m=1-Bz
//   (comparison/checkTilt.C) estimators in here so all tilt estimates live in
//   one macro with a single summary table.

#include <TH1D.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLatex.h>
#include <TString.h>

#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

// ── grid matching sPHENIXFieldMap ────────────────────────────────────────────
static const int    kNR = 37,   kNZ = 241;
static const double kRMin = 0.,  kRMax = 900., kdR = 25.;   // mm
static const double kZMin = -2700., kZMax = 2100., kdZ = 20.; // mm

// Solve a 2x2 normal-equation system [[a,b],[b,c]]{x;y} = {d;e}. Returns false
// if singular (no spread in the regressor).
static bool Solve2(double a, double b, double c, double d, double e,
                   double& x, double& y)
{
    double det = a * c - b * b;
    if (std::abs(det) < 1e-30) return false;
    x = ( c * d - b * e) / det;
    y = (-b * d + a * e) / det;
    return true;
}

void estimateTilt(const char* datadir = "data", const char* outdir = "plots")
{
    const std::string fineCSV  = std::string(datadir) + "/pointCloudFineFullField.csv";
    const std::string roughCSV = std::string(datadir) + "/pointCloudRoughFullField.csv";

    // ── per-(ir,iz) m=0 accumulators (Method 1) ──────────────────────────────
    std::vector<double> mBx(kNR * kNZ, 0), mBy(kNR * kNZ, 0), mBz(kNR * kNZ, 0);
    std::vector<int>    mN (kNR * kNZ, 0);
    auto IDX = [&](int ir, int iz){ return ir * kNZ + iz; };

    // ── per-iz linear-fit accumulators (Method 2): Bx = px + gx*x, By = py+gy*y
    std::vector<double> Sn(kNZ,0);                       // count
    std::vector<double> Sx(kNZ,0), Sxx(kNZ,0), SBx(kNZ,0), SxBx(kNZ,0);
    std::vector<double> Sy(kNZ,0), Syy(kNZ,0), SBy(kNZ,0), SyBy(kNZ,0);
    std::vector<double> SBz(kNZ,0);

    auto readFile = [&](const std::string& fname){
        FILE* f = std::fopen(fname.c_str(), "r");
        if (!f) { std::fprintf(stderr, "Cannot open %s\n", fname.c_str()); return; }
        std::printf("Reading %s ...\n", fname.c_str());
        long long nread = 0, nused = 0;
        char line[512];
        while (std::fgets(line, sizeof(line), f)) {
            if (line[0] == '#' || line[0] == '\n') continue;
            double xs, ys, zs, Bmag, Bxs, Bys, Bzs;
            if (std::sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                            &xs,&ys,&zs,&Bmag,&Bxs,&Bys,&Bzs) != 7) continue;
            ++nread;
            const double xp = xs, yp = zs, zp = -ys;
            const double r  = std::sqrt(xp*xp + yp*yp);
            const double Bx =  Bxs;     // sPHENIX transverse x
            const double By =  Bzs;     // sPHENIX transverse y
            const double Bz = -Bys;     // sPHENIX longitudinal
            if (r < kRMin || r > kRMax + kdR)   continue;
            if (zp < kZMin || zp > kZMax + kdZ) continue;
            int ir = (int)std::lround((r  - kRMin) / kdR);
            int iz = (int)std::lround((zp - kZMin) / kdZ);
            if (ir < 0) ir = 0; if (ir >= kNR) ir = kNR - 1;
            if (iz < 0) iz = 0; if (iz >= kNZ) iz = kNZ - 1;
            int k = IDX(ir, iz);
            mBx[k] += Bx; mBy[k] += By; mBz[k] += Bz; mN[k] += 1;

            Sn[iz]  += 1;
            Sx[iz]  += xp; Sxx[iz] += xp*xp; SBx[iz] += Bx; SxBx[iz] += xp*Bx;
            Sy[iz]  += yp; Syy[iz] += yp*yp; SBy[iz] += By; SyBy[iz] += yp*By;
            SBz[iz] += Bz;
            ++nused;
        }
        std::fclose(f);
        std::printf("  %lld rows read, %lld used\n", nread, nused);
    };
    readFile(fineCSV);
    readFile(roughCSV);

    // ═══════════════════════════════════════════════════════════════════════
    // Method 1: near-axis m=0 transverse field, theta(r) in the central band
    // ═══════════════════════════════════════════════════════════════════════
    const double zBand = 150.0;                  // |z| < 150 mm = flat central field
    const int izLo = (int)std::lround((-zBand - kZMin)/kdZ);
    const int izHi = (int)std::lround(( zBand - kZMin)/kdZ);

    std::printf("\n=== Method 1: near-axis m=0 transverse field (|z| < %.0f mm) ===\n",
                zBand);
    std::printf("  r[mm]   N     <Bx>[mT]  <By>[mT]   <Bz>[T]   theta_x  theta_y  |theta|[mrad]\n");

    const double rInnerMax = 300.0;              // mm: "near-axis" combined region
    TGraph* gThetaR = new TGraph;
    int    irInner = -1;
    double thInner = 0, txInner = 0, tyInner = 0;
    double cBx=0, cBy=0, cBz=0; long long cN=0;  // combined inner-region sums
    std::vector<double> txRing, tyRing;          // per-ring angles (r<=rInnerMax)
    for (int ir = 0; ir < kNR; ++ir) {
        double sBx=0, sBy=0, sBz=0; long long n=0;
        for (int iz = izLo; iz <= izHi; ++iz) {
            int k = IDX(ir, iz);
            if (mN[k] == 0) continue;
            // weight each (r,z) node by its own measurement count
            sBx += mBx[k]; sBy += mBy[k]; sBz += mBz[k]; n += mN[k];
        }
        if (n < 4) continue;                     // require real coverage
        double bx = sBx/n, by = sBy/n, bz = sBz/n;
        double tx = 1e3*std::atan2(bx, bz);
        double ty = 1e3*std::atan2(by, bz);
        double th = 1e3*std::sqrt(bx*bx + by*by)/std::abs(bz);
        double r  = kRMin + ir*kdR;
        std::printf("  %5.0f  %4lld  %+8.3f  %+8.3f  %8.4f  %+7.3f  %+7.3f   %7.3f\n",
                    r, n, bx*1e3, by*1e3, bz, tx, ty, th);
        gThetaR->AddPoint(r, th);
        if (irInner < 0 && r > 0.) {             // smallest measured ring with r>0
            irInner = ir; thInner = th; txInner = tx; tyInner = ty;
        }
        if (r > 0. && r <= rInnerMax) {
            cBx += sBx; cBy += sBy; cBz += sBz; cN += n;
            txRing.push_back(tx); tyRing.push_back(ty);
        }
    }
    double phiInner = std::atan2(tyInner, txInner) * 180.0/M_PI;
    std::printf("  --> innermost measured ring r=%.0f mm:  |theta| = %.3f mrad  "
                "(theta_x=%+.3f, theta_y=%+.3f)  direction phi=%.1f deg\n",
                kRMin + irInner*kdR, thInner, txInner, tyInner, phiInner);
    // Error propagation from (theta_x +/- s_x, theta_y +/- s_y) to (|theta|, phi).
    // All angles in mrad; phi and its error in degrees.
    auto propagate = [](double tx, double ty, double stx, double sty,
                        double& th, double& sth, double& phi, double& sphi){
        double r2 = tx*tx + ty*ty;
        th  = std::sqrt(r2);
        sth = (th  > 0) ? std::sqrt(tx*tx*stx*stx + ty*ty*sty*sty)/th : 0;
        phi = std::atan2(ty, tx)*180.0/M_PI;
        sphi= (r2  > 0) ? std::sqrt(tx*tx*sty*sty + ty*ty*stx*stx)/r2 * 180.0/M_PI : 0;
    };

    // Method-1 estimate = mean of the per-ring angles (r<=rInnerMax); the
    // error is the standard error of that mean across the independent rings,
    // which captures the ring-to-ring scatter (the dominant uncertainty here).
    double cTx=0, cTy=0, cTh=0, cPhi=0, sTx=0, sTy=0, cSth=0, cSphi=0;
    int nRing = (int)txRing.size();
    if (nRing > 1) {
        for (int i=0;i<nRing;++i){ cTx += txRing[i]; cTy += tyRing[i]; }
        cTx /= nRing; cTy /= nRing;
        double vx=0, vy=0;
        for (int i=0;i<nRing;++i){ vx += (txRing[i]-cTx)*(txRing[i]-cTx);
                                   vy += (tyRing[i]-cTy)*(tyRing[i]-cTy); }
        vx /= (nRing-1); vy /= (nRing-1);
        sTx = std::sqrt(vx/nRing); sTy = std::sqrt(vy/nRing);   // SEM across rings
        propagate(cTx, cTy, sTx, sTy, cTh, cSth, cPhi, cSphi);
        std::printf("  --> mean over %d rings (r<=%.0f mm):  |theta| = %.3f +/- %.3f mrad  "
                    "(theta_x=%+.3f+/-%.3f, theta_y=%+.3f+/-%.3f)  phi=%.1f +/- %.1f deg\n",
                    nRing, rInnerMax, cTh, cSth, cTx, sTx, cTy, sTy, cPhi, cSphi);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Method 2: magnetic-axis line + tilt fit (separates tilt from translation)
    // ═══════════════════════════════════════════════════════════════════════
    // Per-z: Bx = px + gx*x, By = py + gy*y, B(z) = <Bz>.
    // Global: px(z) = theta_x*(B - gx*z) - x0*gx   (linear in theta_x, x0).
    std::vector<double> Zc, Bc, gxc, gyc, pxc, pyc, Wc;
    const double zFitMax = 2000.0;               // mm; paraxial region
    for (int iz = 0; iz < kNZ; ++iz) {
        if (Sn[iz] < 6) continue;
        double z = kZMin + iz*kdZ;
        if (std::abs(z) > zFitMax) continue;
        double n = Sn[iz];
        // need transverse spread to determine the slope
        double varx = Sxx[iz] - Sx[iz]*Sx[iz]/n;
        double vary = Syy[iz] - Sy[iz]*Sy[iz]/n;
        if (varx < 1e3 || vary < 1e3) continue;  // < ~30 mm spread: skip
        double px, gx, py, gy;
        if (!Solve2(n, Sx[iz], Sxx[iz], SBx[iz], SxBx[iz], px, gx)) continue;
        if (!Solve2(n, Sy[iz], Syy[iz], SBy[iz], SyBy[iz], py, gy)) continue;
        Zc.push_back(z);
        Bc.push_back(SBz[iz]/n);
        gxc.push_back(gx); gyc.push_back(gy);
        pxc.push_back(px); pyc.push_back(py);
        Wc.push_back(n);
    }

    // Global 2-parameter fits.  Design columns: c1 = B - g*z, c2 = -g.
    // Errors from the least-squares covariance Cov = s^2 (X^T W X)^-1, with the
    // residual variance s^2 = sum_i w_i r_i^2 / (N-2) setting the scale (the
    // per-slab weights w_i = N_i are treated as relative weights).
    auto fitTilt = [&](const std::vector<double>& g, const std::vector<double>& p,
                       double& theta, double& off,
                       double& sTheta, double& sOff)->bool{
        double a=0,b=0,c=0,d=0,e=0;            // [[a,b],[b,c]]{theta;off}={d;e}
        for (size_t i=0;i<Zc.size();++i){
            double w = Wc[i];
            double u =  Bc[i] - g[i]*Zc[i];    // regressor for theta
            double v = -g[i];                  // regressor for off (=x0/y0)
            a += w*u*u; b += w*u*v; c += w*v*v;
            d += w*u*p[i]; e += w*v*p[i];
        }
        if (!Solve2(a,b,c,d,e,theta,off)) return false;
        double chi2=0; int N=0;
        for (size_t i=0;i<Zc.size();++i){
            double u =  Bc[i] - g[i]*Zc[i], v = -g[i];
            double res = p[i] - (theta*u + off*v);
            chi2 += Wc[i]*res*res; ++N;
        }
        double det = a*c - b*b;
        double s2  = (N>2) ? chi2/(N-2) : 0.0;     // weighted residual variance
        sTheta = (det!=0) ? std::sqrt(s2 * c/det) : 0.0;  // [Cov]_theta,theta
        sOff   = (det!=0) ? std::sqrt(s2 * a/det) : 0.0;  // [Cov]_off,off
        return true;
    };
    double thx=0,thy=0,x0=0,y0=0, sthx=0,sthy=0,sx0=0,sy0=0;
    bool okx = fitTilt(gxc, pxc, thx, x0, sthx, sx0);
    bool oky = fitTilt(gyc, pyc, thy, y0, sthy, sy0);

    std::printf("\n=== Method 2: magnetic-axis line + tilt fit "
                "(|z| < %.0f mm, %zu z-slabs) ===\n", zFitMax, Zc.size());
    double th2=0, sth2=0, phi2=0, sphi2=0;
    if (okx && oky) {
        propagate(thx*1e3, thy*1e3, sthx*1e3, sthy*1e3, th2, sth2, phi2, sphi2);
        std::printf("  theta_x = %+8.3f +/- %.3f mrad   theta_y = %+8.3f +/- %.3f mrad\n",
                    thx*1e3, sthx*1e3, thy*1e3, sthy*1e3);
        std::printf("  |theta| = %8.3f +/- %.3f mrad    tilt direction phi = %.1f +/- %.1f deg\n",
                    th2, sth2, phi2, sphi2);
        std::printf("  axis offset at z=0:  x0 = %+.2f +/- %.2f mm   y0 = %+.2f +/- %.2f mm\n",
                    x0, sx0, y0, sy0);
    } else {
        std::printf("  fit failed (insufficient transverse spread)\n");
    }

    std::printf("\n--- Side by side ---\n");
    std::printf("  Method 1 (near-axis, mean of %d rings): |theta| = %.2f +/- %.2f mrad, "
                "phi = %.1f +/- %.1f deg\n", nRing, cTh, cSth, cPhi, cSphi);
    if (okx && oky)
        std::printf("  Method 2 (axis-line fit)             : |theta| = %.2f +/- %.2f mrad, "
                    "phi = %.1f +/- %.1f deg\n", th2, sth2, phi2, sphi2);

    // ── plots ────────────────────────────────────────────────────────────────
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("cTilt","Tilt estimates",1300,560);
    c->Divide(2,1);

    c->cd(1);
    gThetaR->SetTitle("Method 1: |#theta|(r) per ring in central band (|z|<150 mm);"
                      "r [mm];|#theta| [mrad]");
    gThetaR->SetMarkerStyle(20); gThetaR->SetMarkerColor(kBlue+1);
    gThetaR->SetLineColor(kBlue+1); gThetaR->SetLineWidth(2);
    gThetaR->Draw("ALP");
    {
        // headline = combined inner-region average (single rings are noise-limited)
        TLine* lin = new TLine(0, cTh, kRMax, cTh);
        lin->SetLineStyle(2); lin->SetLineColor(kRed+1); lin->SetLineWidth(2);
        lin->Draw();
        TLatex t; t.SetTextSize(0.034); t.SetTextColor(kRed+1);
        t.DrawLatex(0.04*kRMax, cTh + 0.04*gThetaR->GetYaxis()->GetXmax(),
                    Form("mean of %d rings: |#theta| = %.2f #pm %.2f mrad (#phi=%.0f#circ)",
                         nRing, cTh, cSth, cPhi));
        TLatex t2; t2.SetTextSize(0.030); t2.SetTextColor(kGray+2);
        t2.DrawLatex(0.04*kRMax, 0.90*gThetaR->GetYaxis()->GetXmax(),
                     "(single rings are noise-limited; ~5 mT signal)");
    }

    c->cd(2);
    TGraph* gPx = new TGraph; TGraph* gPy = new TGraph;
    for (size_t i=0;i<Zc.size();++i){ gPx->AddPoint(Zc[i],pxc[i]*1e3);
                                      gPy->AddPoint(Zc[i],pyc[i]*1e3); }
    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle("Method 2: transverse-field intercept p(z) (data) + model fit;"
                 "z [mm];p [mT]");
    gPx->SetMarkerStyle(20); gPx->SetMarkerColor(kBlue+1); gPx->SetMarkerSize(0.6);
    gPy->SetMarkerStyle(21); gPy->SetMarkerColor(kRed+1);  gPy->SetMarkerSize(0.6);
    mg->Add(gPx,"P"); mg->Add(gPy,"P");
    if (okx && oky) {
        TGraph* fPx = new TGraph; TGraph* fPy = new TGraph;
        for (size_t i=0;i<Zc.size();++i){
            fPx->AddPoint(Zc[i], (thx*(Bc[i]-gxc[i]*Zc[i]) - x0*gxc[i])*1e3);
            fPy->AddPoint(Zc[i], (thy*(Bc[i]-gyc[i]*Zc[i]) - y0*gyc[i])*1e3);
        }
        fPx->SetLineColor(kBlue+1); fPx->SetLineWidth(2);
        fPy->SetLineColor(kRed+1);  fPy->SetLineWidth(2);
        mg->Add(fPx,"L"); mg->Add(fPy,"L");
    }
    mg->Draw("A");
    TLegend* leg = new TLegend(0.15,0.75,0.45,0.88);
    leg->AddEntry(gPx,"p_{x}(z) data","p");
    leg->AddEntry(gPy,"p_{y}(z) data","p");
    if (okx && oky) {
        leg->AddEntry((TObject*)nullptr,"lines = model fit","");
    }
    leg->Draw();
    if (okx && oky) {
        TLatex t; t.SetNDC(); t.SetTextSize(0.034); t.SetTextColor(kBlack);
        t.DrawLatex(0.45, 0.86,
            Form("fit: |#theta| = %.2f #pm %.2f mrad  (#phi = %.0f#circ)", th2, sth2, phi2));
        t.DrawLatex(0.45, 0.81,
            Form("#theta_{x}=%+.2f#pm%.2f, #theta_{y}=%+.2f#pm%.2f mrad",
                 thx*1e3, sthx*1e3, thy*1e3, sthy*1e3));
        t.DrawLatex(0.45, 0.76,
            Form("axis offset (x_{0},y_{0}) = (%+.1f#pm%.1f, %+.1f#pm%.1f) mm",
                 x0, sx0, y0, sy0));
    }

    c->SaveAs(Form("%s/tilt_estimates.pdf", outdir));
    c->SaveAs(Form("%s/tilt_estimates.png", outdir));
    std::printf("\nSaved tilt_estimates.{pdf,png} to %s/\n", outdir);
}
