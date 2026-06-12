// analysis/checkHalfField.C
// Estimate the solenoid tilt from the HALF-FIELD (half-excitation) mapping
// campaign and compare it to the FULL-FIELD rough map.
//
// Motivation.  Fine vs. rough at full field is not an independent cross-check of
// the ~4.4 mrad tilt: both scans share the same probe alignment, mapper
// registration and any field-dependent deflection (see README "Systematics").
// The half-field rough map is a genuinely new lever -- a SECOND EXCITATION taken
// on the *same 10 cm grid* as the full-field rough map (verified: identical
// 41760-point sampling, plateau B ratio = 0.50).  Comparing the apparent tilt
// angle theta = atan2(<B_transverse>, <Bz>) between the two excitations tests how
// the signal scales with field, which discriminates:
//
//   hypothesis                              <B_perp>   theta(half) vs theta(full)
//   --------------------------------------  ---------  --------------------------
//   rigid geometric rotation (magnet yaw,   ~ Bz       unchanged
//     frame yaw, probe-triad mount, regn.)
//   probe transverse<->axial cross-talk     ~ Bz       unchanged
//   additive external / ambient field       ~ const    ~doubles
//   magnetic-force mapping-arm deflection    field-dep  decreases
//
//   What this CAN do: confirm the tilt is a linear-in-Bz rotation/cross-talk and
//     detect/exclude an additive external field or a field-dependent deflection.
//   What this CANNOT do: break the magnet-yaw vs measurement-frame-yaw
//     degeneracy -- both are field-independent rotations (theta unchanged).  Only
//     the survey resolves that.
//
// Method is identical to analysis/checkAlignment.C (global phi-averaged
// transverse field), applied to each map through the same code path so any
// method bias cancels in the comparison.  The only adaptation is that the
// strong-field z-slice cut used for the dispersion error is scaled to each map's
// own plateau (0.70 x max|Bz|), so the half map's ~0.7 T plateau is not rejected
// by the full-field's fixed 1.0 T cut.  The full-field analysis is untouched.
//
// CSV (surveyor frame): x_s,y_s,z_s,|B|,Bx_s,By_s,Bz_s  [mm, T]
//   sPHENIX transform: xp=xs, yp=zs, zp=-ys; Bx=Bxs, By=Bzs, Bz=-Bys
//
// Run from the repository root:
//   root -l -b -q 'analysis/checkHalfField.C+("data","plots")'

#include <TH1D.h>
#include <TGraph.h>
#include <TMarker.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TString.h>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>

// z-binning matches checkAlignment.C / estimateTilt.C
static const double kZMin = -2700.0, kZMax = 2100.0, kdZ = 20.0;
static const int    kNZ   = int((kZMax - kZMin) / kdZ) + 1;   // 241

struct TiltResult {
    bool   ok       = false;
    long long gN    = 0;
    double gBx=0, gBy=0, gBz=0;       // global phi-averaged field [T]
    double bperp    = 0;              // sqrt(gBx^2+gBy^2) [T]
    double thx=0, thy=0, th=0, phi=0; // global tilt [mrad] / [deg]
    double sth=0, sphi=0;             // z-slice dispersion error [mrad]/[deg]
    int    nS       = 0;              // # strong-field slices
    double platMax  = 0;             // max |<Bz>| over z-slices [T]
};

// Run the global phi-averaged transverse-field estimator on one CSV.  Fills the
// optional per-z histograms (theta_x [mrad] and |B_perp| [mT]) if supplied.
static TiltResult analyzeMap(const std::string& csv, double strongFrac,
                             TH1D* hThx = nullptr, TH1D* hBperp = nullptr)
{
    TiltResult R;
    std::vector<double> sumBx(kNZ,0), sumBy(kNZ,0), sumBz(kNZ,0);
    std::vector<long long> npts(kNZ,0);

    std::ifstream f(csv);
    if (!f.is_open()) { printf("Cannot open %s\n", csv.c_str()); return R; }
    printf("Reading %s ...\n", csv.c_str());

    long long nread = 0;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        double xs, ys, zs, Bmag, Bxs, Bys, Bzs;
        if (std::sscanf(line.c_str(), "%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                        &xs,&ys,&zs,&Bmag,&Bxs,&Bys,&Bzs) != 7) continue;
        ++nread;
        const double zp   = -ys;          // sPHENIX z
        const double Bx_p =  Bxs;          // Bx_phx = Bx_s
        const double By_p =  Bzs;          // By_phx = Bz_s
        const double Bz_p = -Bys;          // Bz_phx = -By_s
        int iz = int((zp - kZMin) / kdZ);
        if (iz < 0 || iz >= kNZ) continue;
        sumBx[iz] += Bx_p; sumBy[iz] += By_p; sumBz[iz] += Bz_p; npts[iz] += 1;
    }
    printf("  read %lld rows\n", nread);

    // global means over all measured points
    double gBx=0,gBy=0,gBz=0; long long gN=0;
    for (int iz=0; iz<kNZ; ++iz) {
        if (npts[iz]==0) continue;
        gBx += sumBx[iz]; gBy += sumBy[iz]; gBz += sumBz[iz]; gN += npts[iz];
        double bz = sumBz[iz]/npts[iz];
        if (std::abs(bz) > R.platMax) R.platMax = std::abs(bz);
        if (hThx || hBperp) {
            double bx = sumBx[iz]/npts[iz], by = sumBy[iz]/npts[iz];
            if (hThx && bz != 0)
                hThx->SetBinContent(iz+1, 1e3*std::atan2(bx, bz));
            if (hBperp)
                hBperp->SetBinContent(iz+1, 1e3*std::sqrt(bx*bx + by*by));
        }
    }
    if (gN == 0) return R;
    gBx/=gN; gBy/=gN; gBz/=gN;
    R.ok=true; R.gN=gN; R.gBx=gBx; R.gBy=gBy; R.gBz=gBz;
    R.bperp = std::sqrt(gBx*gBx + gBy*gBy);
    R.thx = 1e3*std::atan2(gBx,gBz);
    R.thy = 1e3*std::atan2(gBy,gBz);
    R.th  = 1e3*std::sqrt(gBx*gBx+gBy*gBy)/std::abs(gBz);
    R.phi = std::atan2(gBy,gBx)*180.0/M_PI;

    // dispersion error from independent strong-field z-slices, threshold scaled
    // to this map's own plateau so half field is not rejected
    const double cut = strongFrac * R.platMax;
    std::vector<double> txS, tyS; double mtx=0, mty=0; int nS=0;
    for (int iz=0; iz<kNZ; ++iz) {
        if (npts[iz]==0) continue;
        double bx=sumBx[iz]/npts[iz], by=sumBy[iz]/npts[iz], bz=sumBz[iz]/npts[iz];
        if (std::abs(bz) < cut) continue;
        double tx=1e3*std::atan2(bx,bz), ty=1e3*std::atan2(by,bz);
        txS.push_back(tx); tyS.push_back(ty); mtx+=tx; mty+=ty; ++nS;
    }
    R.nS = nS;
    if (nS>1) {
        mtx/=nS; mty/=nS;
        double vx=0, vy=0;
        for (int i=0;i<nS;++i){ vx+=(txS[i]-mtx)*(txS[i]-mtx);
                                vy+=(tyS[i]-mty)*(tyS[i]-mty); }
        vx/=(nS-1); vy/=(nS-1);
        double semx=std::sqrt(vx/nS), semy=std::sqrt(vy/nS);
        double th=std::sqrt(mtx*mtx+mty*mty);
        R.sth  = (th>0)?std::sqrt(mtx*mtx*semx*semx+mty*mty*semy*semy)/th:0;
        double r2=mtx*mtx+mty*mty;
        R.sphi = (r2>0)?std::sqrt(mtx*mtx*semy*semy+mty*mty*semx*semx)/r2*180.0/M_PI:0;
    }
    return R;
}

static void report(const char* tag, const TiltResult& R, double cutFrac)
{
    if (!R.ok) { printf("  %s: no data\n", tag); return; }
    printf("\n== %s ==\n", tag);
    printf("  N points          : %lld\n", R.gN);
    printf("  plateau max|<Bz>| : %.4f T   (strong-slice cut = %.2f x plateau = %.3f T, %d slices)\n",
           R.platMax, cutFrac, cutFrac*R.platMax, R.nS);
    printf("  <Bx> <By> <Bz>    : %+.4f  %+.4f mT   %.4f T\n",
           R.gBx*1e3, R.gBy*1e3, R.gBz);
    printf("  <B_perp>          : %.4f mT\n", R.bperp*1e3);
    printf("  theta_x theta_y   : %+.3f  %+.3f mrad\n", R.thx, R.thy);
    printf("  |theta|           : %.3f +/- %.3f mrad   (phi = %.1f +/- %.1f deg)\n",
           R.th, R.sth, R.phi, R.sphi);
}

// ── Magnetic centre: transverse axis (x0,y0) + tilt, and axial centre z0 ──────
// Same model as analysis/estimateTilt.C Method 2 (which is left untouched): per
// z-slice fit Bx = px + gx*x, By = py + gy*y; then a global fit
//   px(z) = thx*(B - gx*z) - x0*gx
// separates the tilt thx from the transverse axis offset x0 (likewise y).  A
// field-INDEPENDENT additive transverse field projects onto the tilt regressor
// (mean <B> != 0) but NOT onto the focusing regressor g = -1/2 dBz/dz (g is odd
// about the centre, <g> ~ 0).  So it biases the fitted TILT by ~1/<B> -- doubling
// at half field -- while leaving the focusing-derived axis (x0,y0) essentially
// unchanged.  Hence a stable centre together with a growing apparent tilt is the
// signature of an additive offset; a real geometric centre is field-independent.
static bool Solve2(double a,double b,double c,double d,double e,double& x,double& y){
    double det=a*c-b*b; if(std::abs(det)<1e-300) return false;
    x=( c*d - b*e)/det; y=(-b*d + a*e)/det; return true;
}

struct CenterResult {
    bool   ok=false;
    double thx=0,thy=0,x0=0,y0=0,sx0=0,sy0=0,sthx=0,sthy=0;
    double z0=0;        // axial centre from <Bz>(z) parabola vertex [mm]
    double bzPk=0;      // peak all-r <Bz> [T]
};

static CenterResult analyzeCenter(const std::string& csv)
{
    CenterResult R;
    const double rMax=900., dR=25.;
    std::vector<double> Sn(kNZ,0),Sx(kNZ,0),Sxx(kNZ,0),SBx(kNZ,0),SxBx(kNZ,0);
    std::vector<double> Sy(kNZ,0),Syy(kNZ,0),SBy(kNZ,0),SyBy(kNZ,0),SBz(kNZ,0);

    FILE* f=std::fopen(csv.c_str(),"r");
    if(!f){ printf("Cannot open %s\n",csv.c_str()); return R; }
    char line[512];
    while(std::fgets(line,sizeof(line),f)){
        if(line[0]=='#'||line[0]=='\n') continue;
        double xs,ys,zs,Bmag,Bxs,Bys,Bzs;
        if(std::sscanf(line,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                       &xs,&ys,&zs,&Bmag,&Bxs,&Bys,&Bzs)!=7) continue;
        const double xp=xs, yp=zs, zp=-ys;
        const double r=std::sqrt(xp*xp+yp*yp);
        const double Bx=Bxs, By=Bzs, Bz=-Bys;
        if(r>rMax+dR) continue;
        int iz=int((zp-kZMin)/kdZ);
        if(iz<0||iz>=kNZ) continue;
        Sn[iz]+=1; Sx[iz]+=xp; Sxx[iz]+=xp*xp; SBx[iz]+=Bx; SxBx[iz]+=xp*Bx;
        Sy[iz]+=yp; Syy[iz]+=yp*yp; SBy[iz]+=By; SyBy[iz]+=yp*By; SBz[iz]+=Bz;
    }
    std::fclose(f);

    std::vector<double> Zc,Bc,gxc,gyc,pxc,pyc,Wc;
    const double zFitMax=2000.0;
    for(int iz=0;iz<kNZ;++iz){
        if(Sn[iz]<6) continue;
        double z=kZMin+iz*kdZ;
        if(std::abs(z)>zFitMax) continue;
        double n=Sn[iz];
        double varx=Sxx[iz]-Sx[iz]*Sx[iz]/n;
        double vary=Syy[iz]-Sy[iz]*Sy[iz]/n;
        if(varx<1e3||vary<1e3) continue;
        double px,gx,py,gy;
        if(!Solve2(n,Sx[iz],Sxx[iz],SBx[iz],SxBx[iz],px,gx)) continue;
        if(!Solve2(n,Sy[iz],Syy[iz],SBy[iz],SyBy[iz],py,gy)) continue;
        Zc.push_back(z); Bc.push_back(SBz[iz]/n);
        gxc.push_back(gx); gyc.push_back(gy);
        pxc.push_back(px); pyc.push_back(py); Wc.push_back(n);
    }
    auto fitTilt=[&](const std::vector<double>& g,const std::vector<double>& p,
                     double& theta,double& off,double& sTheta,double& sOff)->bool{
        double a=0,b=0,c=0,d=0,e=0;
        for(size_t i=0;i<Zc.size();++i){
            double w=Wc[i], u=Bc[i]-g[i]*Zc[i], v=-g[i];
            a+=w*u*u; b+=w*u*v; c+=w*v*v; d+=w*u*p[i]; e+=w*v*p[i];
        }
        if(!Solve2(a,b,c,d,e,theta,off)) return false;
        double chi2=0; int N=0;
        for(size_t i=0;i<Zc.size();++i){
            double u=Bc[i]-g[i]*Zc[i], v=-g[i];
            double res=p[i]-(theta*u+off*v); chi2+=Wc[i]*res*res; ++N;
        }
        double det=a*c-b*b, s2=(N>2)?chi2/(N-2):0.0;
        sTheta=(det!=0)?std::sqrt(s2*c/det):0.0;
        sOff  =(det!=0)?std::sqrt(s2*a/det):0.0;
        return true;
    };
    bool okx=fitTilt(gxc,pxc,R.thx,R.x0,R.sthx,R.sx0);
    bool oky=fitTilt(gyc,pyc,R.thy,R.y0,R.sthy,R.sy0);
    R.ok = okx && oky;

    // axial centre: parabola vertex of the all-r <Bz>(z) near its peak
    int izPk=-1; double bzPk=-1e9;
    for(int iz=0;iz<kNZ;++iz){ if(Sn[iz]<1) continue;
        double bz=SBz[iz]/Sn[iz]; if(bz>bzPk){bzPk=bz;izPk=iz;} }
    R.bzPk=bzPk;
    if(izPk>=0){
        double zPk=kZMin+izPk*kdZ, W=600.0;
        double S0=0,S1=0,S2=0,S3=0,S4=0,T0=0,T1=0,T2=0;
        for(int iz=0;iz<kNZ;++iz){ if(Sn[iz]<1) continue;
            double z=kZMin+iz*kdZ; if(std::abs(z-zPk)>W) continue;
            double y=SBz[iz]/Sn[iz], z2=z*z;
            S0+=1;S1+=z;S2+=z2;S3+=z2*z;S4+=z2*z2; T0+=y;T1+=y*z;T2+=y*z2;
        }
        double m[3][4]={{S0,S1,S2,T0},{S1,S2,S3,T1},{S2,S3,S4,T2}};
        for(int col=0;col<3;++col){
            int piv=col; for(int rr=col+1;rr<3;++rr) if(std::abs(m[rr][col])>std::abs(m[piv][col])) piv=rr;
            for(int j=0;j<4;++j) std::swap(m[col][j],m[piv][j]);
            double dd=m[col][col]; for(int j=col;j<4;++j) m[col][j]/=dd;
            for(int rr=0;rr<3;++rr) if(rr!=col){ double fac=m[rr][col];
                for(int j=col;j<4;++j) m[rr][j]-=fac*m[col][j]; }
        }
        double a1=m[1][3],a2=m[2][3];
        R.z0=(a2!=0.0)?-a1/(2.0*a2):0.0;
    }
    return R;
}

// ── Pooled <B_perp> vs Bz slice fit (the "through-zero linearity" test) ───────
// Per z-slice, the phi-averaged transverse field obeys <Bx> = a*Bz + c, where
// a = real field-PROPORTIONAL tilt and c = a current-INDEPENDENT transverse term
// (probe pedestal, remanence, or ambient field).  A genuine rigid tilt has c=0
// (the line passes through the origin).  Restricting to the flat-field region
// kills the focusing term g*x0 (g = -1/2 dBz/dz ~ 0 there), so the intercept is
// clean.  Pooling full + half field gives two Bz clusters that pin the slope and
// the intercept simultaneously -- which a single excitation cannot do.
struct Slice { double z, bx, by, bz; double n; };

static std::vector<Slice> collectSlices(const std::string& csv)
{
    std::vector<double> sBx(kNZ,0), sBy(kNZ,0), sBz(kNZ,0), nn(kNZ,0);
    std::ifstream f(csv);
    if (!f.is_open()) { printf("Cannot open %s\n", csv.c_str()); return {}; }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0]=='#') continue;
        double xs,ys,zs,Bmag,Bxs,Bys,Bzs;
        if (std::sscanf(line.c_str(), "%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                        &xs,&ys,&zs,&Bmag,&Bxs,&Bys,&Bzs) != 7) continue;
        double zp = -ys; int iz = int((zp - kZMin)/kdZ);
        if (iz<0 || iz>=kNZ) continue;
        sBx[iz]+=Bxs; sBy[iz]+=Bzs; sBz[iz]+=-Bys; nn[iz]+=1;
    }
    std::vector<Slice> out;
    for (int iz=0; iz<kNZ; ++iz) if (nn[iz]>0)
        out.push_back({kZMin+iz*kdZ, sBx[iz]/nn[iz], sBy[iz]/nn[iz], sBz[iz]/nn[iz], nn[iz]});
    return out;
}

// weighted least squares y = a*x + c; errors from the residual scatter
static void linfit(const std::vector<double>& x, const std::vector<double>& y,
                   const std::vector<double>& w,
                   double& a, double& c, double& sa, double& sc, double& rms)
{
    double Sw=0,Sx=0,Sy=0,Sxx=0,Sxy=0;
    for (size_t i=0;i<x.size();++i){ double wi=w[i];
        Sw+=wi; Sx+=wi*x[i]; Sy+=wi*y[i]; Sxx+=wi*x[i]*x[i]; Sxy+=wi*x[i]*y[i]; }
    double det = Sw*Sxx - Sx*Sx;
    a = (Sw*Sxy - Sx*Sy)/det;
    c = (Sxx*Sy - Sx*Sxy)/det;
    double chi2=0, Wsum=0; int N=(int)x.size();
    for (size_t i=0;i<x.size();++i){ double r=y[i]-(a*x[i]+c); chi2+=w[i]*r*r; Wsum+=w[i]; }
    double s2 = (N>2)? chi2/(N-2) : 0.0;
    rms = std::sqrt(chi2/Wsum);
    sa = std::sqrt(s2 * Sw  / det);
    sc = std::sqrt(s2 * Sxx / det);
}

void checkHalfField(const char* datadir = "data", const char* outdir = "plots")
{
    const double strongFrac = 0.70;   // same fraction of plateau for both maps
    const std::string fullCSV = std::string(datadir) + "/pointCloudRoughFullField.csv";
    const std::string halfCSV = std::string(datadir) + "/pointCloudRoughHalfField.csv";

    TH1D* hThxF = new TH1D("hThxF",";z [mm];#theta_{x} [mrad]", kNZ, kZMin, kZMax);
    TH1D* hThxH = new TH1D("hThxH",";z [mm];#theta_{x} [mrad]", kNZ, kZMin, kZMax);
    TH1D* hBpF  = new TH1D("hBpF", ";z [mm];#LT|B_{#perp}|#GT [mT]", kNZ, kZMin, kZMax);
    TH1D* hBpH  = new TH1D("hBpH", ";z [mm];#LT|B_{#perp}|#GT [mT]", kNZ, kZMin, kZMax);

    printf("\n========== tilt from full- vs half-field rough maps ==========\n");
    TiltResult Rf = analyzeMap(fullCSV, strongFrac, hThxF, hBpF);
    TiltResult Rh = analyzeMap(halfCSV, strongFrac, hThxH, hBpH);

    report("FULL field (rough)", Rf, strongFrac);
    report("HALF field (rough)", Rh, strongFrac);

    if (Rf.ok && Rh.ok) {
        const double rBz    = Rh.platMax / Rf.platMax;            // excitation ratio
        const double rBperp = Rh.bperp   / Rf.bperp;             // <B_perp> ratio
        const double rTheta = Rh.th      / Rf.th;                // tilt-angle ratio
        printf("\n== comparison (half / full) ==\n");
        printf("  excitation (plateau Bz) ratio : %.3f\n", rBz);
        printf("  <B_perp> ratio                : %.3f   (rotation/cross-talk -> ~%.3f)\n",
               rBperp, rBz);
        printf("  |theta| ratio                 : %.3f   (rotation/cross-talk -> ~1.00)\n",
               rTheta);
        printf("\n  Interpretation guide:\n");
        printf("    theta ratio ~1, B_perp ratio ~Bz ratio  -> linear-in-Bz rotation or cross-talk\n");
        printf("    theta ratio >1 (B_perp ~constant)        -> additive external/ambient field\n");
        printf("    theta ratio <1 (B_perp falls faster)     -> field-dependent mechanical deflection\n");
        printf("  NB: a 'linear-in-Bz rotation' result still cannot separate a magnet\n");
        printf("      yaw from a measurement-frame yaw -- only the survey does that.\n");
    }

    // ── magnetic centre: transverse axis (x0,y0) and axial centre z0 ──────────
    printf("\n========== magnetic centre: full vs half field ==========\n");
    CenterResult Cf = analyzeCenter(fullCSV);
    CenterResult Ch = analyzeCenter(halfCSV);
    auto repC=[&](const char* tag,const CenterResult& C){
        if(!C.ok){ printf("  %s: axis fit failed\n",tag); return; }
        printf("  %s  (peak <Bz> = %.4f T):\n", tag, C.bzPk);
        printf("    transverse axis at z=0 : x0 = %+.2f +/- %.2f mm   y0 = %+.2f +/- %.2f mm\n",
               C.x0,C.sx0,C.y0,C.sy0);
        printf("    tilt (M2 axis-line fit): theta_x = %+.3f +/- %.3f mrad   theta_y = %+.3f +/- %.3f mrad\n",
               C.thx*1e3,C.sthx*1e3,C.thy*1e3,C.sthy*1e3);
        printf("    axial centre z0        : %.1f mm\n",C.z0);
    };
    repC("FULL field (rough)",Cf);
    repC("HALF field (rough)",Ch);
    if(Cf.ok && Ch.ok){
        printf("\n  centre shift (half - full):  dx0 = %+.2f mm   dy0 = %+.2f mm   dz0 = %+.1f mm\n",
               Ch.x0-Cf.x0, Ch.y0-Cf.y0, Ch.z0-Cf.z0);
        printf("    A geometric centre is field-independent (dx0,dy0,dz0 ~ 0 within errors).\n");
        printf("    An additive transverse offset biases the fitted TILT (not the centre): it\n");
        printf("    grows ~1/Bz, so theta increasing while x0,y0,z0 hold steady is its signature.\n");
        printf("    M2 theta: %+.3f -> %+.3f mrad (theta_x), centre stable -> additive term ~%.2f mrad at full field.\n",
               Cf.thx*1e3, Ch.thx*1e3, (Ch.thx-Cf.thx)*1e3);
    }

    // ── pooled <B_perp> vs Bz slice fit (through-zero linearity test) ─────────
    printf("\n========== pooled <B_perp> vs Bz flat-field slice fit ==========\n");
    const double zFlat = 700.0;     // |z| < zFlat mm: flat field, focusing g ~ 0
    std::vector<Slice> Sf = collectSlices(fullCSV);
    std::vector<Slice> Sh = collectSlices(halfCSV);
    std::vector<double> X, Yx, Yy, W;
    TGraph* gFx=new TGraph; TGraph* gHx=new TGraph;
    TGraph* gFy=new TGraph; TGraph* gHy=new TGraph;
    int nF=0,nH=0;
    auto addBand=[&](std::vector<Slice>& S, int& cnt, TGraph* gx, TGraph* gy){
        for (auto& s : S) {
            if (std::abs(s.z) > zFlat) continue;
            X.push_back(s.bz); Yx.push_back(s.bx); Yy.push_back(s.by); W.push_back(s.n);
            gx->AddPoint(s.bz, s.bx*1e3); gy->AddPoint(s.bz, s.by*1e3); ++cnt;
        }
    };
    addBand(Sf, nF, gFx, gFy);
    addBand(Sh, nH, gHx, gHy);

    double ax,cx,sax,scx,rmsx, ay,cy,say,scy,rmsy;
    linfit(X, Yx, W, ax, cx, sax, scx, rmsx);
    linfit(X, Yy, W, ay, cy, say, scy, rmsy);
    double bzlo=*std::min_element(X.begin(),X.end());
    double bzhi=*std::max_element(X.begin(),X.end());
    printf("  flat-field slices (|z|<%.0f mm): %d full + %d half = %zu points, Bz in [%.3f, %.3f] T\n",
           zFlat, nF, nH, X.size(), bzlo, bzhi);
    printf("  <Bx> = a*Bz + c :  a = %+.3f +/- %.3f mrad    c = %+.4f +/- %.4f mT   (resid rms %.4f mT)\n",
           ax*1e3, sax*1e3, cx*1e3, scx*1e3, rmsx*1e3);
    printf("  <By> = a*Bz + c :  a = %+.3f +/- %.3f mrad    c = %+.4f +/- %.4f mT   (resid rms %.4f mT)\n",
           ay*1e3, say*1e3, cy*1e3, scy*1e3, rmsy*1e3);
    double thP  = 1e3*std::sqrt(ax*ax+ay*ay);
    double phiP = std::atan2(ay,ax)*180.0/M_PI;
    double off  = 1e3*std::sqrt(cx*cx+cy*cy);
    double phiO = std::atan2(cy,cx)*180.0/M_PI;
    printf("\n  => field-PROPORTIONAL tilt (slope)      : |theta| = %.3f mrad   phi = %.1f deg\n",
           thP, phiP);
    printf("     current-INDEPENDENT offset (intercept) : |B|     = %.4f mT     phi = %.1f deg\n",
           off, phiO);
    printf("     intercept significance: c_x/sigma = %.1f, c_y/sigma = %.1f\n",
           (scx!=0?cx/scx:0), (scy!=0?cy/scy:0));
    printf("     (a rigid tilt predicts c = 0; nonzero c = current-independent transverse field)\n");

    // money plot: <Bx>, <By> vs Bz with both excitations and the fit line
    gStyle->SetOptStat(0);
    TCanvas* cF = new TCanvas("cFit","B_perp vs Bz",1200,500);
    cF->Divide(2,1);
    auto drawFit=[&](int pad,TGraph* gF,TGraph* gH,double a,double c,
                     const char* ylab){
        cF->cd(pad);
        gF->SetMarkerStyle(20); gF->SetMarkerColor(kBlue+1); gF->SetMarkerSize(0.5);
        gH->SetMarkerStyle(21); gH->SetMarkerColor(kRed+1);  gH->SetMarkerSize(0.5);
        gF->SetTitle(Form(";#LT B_{z} #GT [T];%s", ylab));
        gF->GetXaxis()->SetLimits(0.0, 1.5);
        double y0=c*1e3, y1=(a*1.5+c)*1e3;
        double ymin=std::min({y0,y1,0.0})-0.5, ymax=std::max({y0,y1,0.0})+0.5;
        gF->GetYaxis()->SetRangeUser(ymin,ymax);
        gF->Draw("AP"); gH->Draw("P SAME");
        TLine* lf=new TLine(0.0,y0,1.5,y1); lf->SetLineColor(kBlack); lf->SetLineWidth(2); lf->Draw();
        TLine* l0=new TLine(0.0,0.0,1.5,0.0); l0->SetLineStyle(3); l0->Draw();
        TMarker* mc=new TMarker(0.0,y0,29); mc->SetMarkerColor(kBlack); mc->SetMarkerSize(1.6); mc->Draw();
        // data fall in the lower half (<B_perp> grows negative with Bz), so the
        // legend goes top-right where the plot is empty
        TLegend* lg=new TLegend(0.58,0.70,0.88,0.88);
        lg->AddEntry(gF,"full field","p"); lg->AddEntry(gH,"half field","p");
        lg->AddEntry(mc,"intercept c (Bz=0)","p"); lg->Draw();
    };
    drawFit(1,gFx,gHx,ax,cx,"#LT B_{x} #GT [mT]");
    drawFit(2,gFy,gHy,ay,cy,"#LT B_{y} #GT [mT]");
    cF->SaveAs(Form("%s/halffield_BperpVsBz.pdf", outdir));
    cF->SaveAs(Form("%s/halffield_BperpVsBz.png", outdir));
    printf("Saved halffield_BperpVsBz.{pdf,png}\n");

    // overlay plots
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("cHalf","tilt: full vs half field",1200,500);
    c->Divide(2,1);
    auto style=[](TH1D* h,int col){ h->SetLineColor(col); h->SetLineWidth(2); };

    c->cd(1);
    style(hThxF,kBlue+1); style(hThxH,kRed+1);
    hThxF->SetTitle("#theta_{x}(z): full (blue) vs half (red) field;z [mm];#theta_{x} [mrad]");
    hThxF->GetXaxis()->SetRangeUser(-1500,1500);
    hThxF->SetMinimum(-10); hThxF->SetMaximum(10);
    hThxF->Draw("HIST"); hThxH->Draw("HIST SAME");
    TLine* l0=new TLine(-1500,0,1500,0); l0->SetLineStyle(2); l0->Draw();
    { TLegend* lg=new TLegend(0.60,0.75,0.88,0.88);
      lg->AddEntry(hThxF,"full field","l"); lg->AddEntry(hThxH,"half field","l");
      lg->Draw(); }

    c->cd(2);
    style(hBpF,kBlue+1); style(hBpH,kRed+1);
    hBpF->SetTitle("#LT|B_{#perp}|#GT(z): full (blue) vs half (red);z [mm];#LT|B_{#perp}|#GT [mT]");
    hBpF->GetXaxis()->SetRangeUser(-1500,1500);
    hBpF->Draw("HIST"); hBpH->Draw("HIST SAME");
    { TLegend* lg=new TLegend(0.60,0.75,0.88,0.88);
      lg->AddEntry(hBpF,"full field","l"); lg->AddEntry(hBpH,"half field","l");
      lg->Draw(); }

    c->SaveAs(Form("%s/halffield_tilt.pdf", outdir));
    c->SaveAs(Form("%s/halffield_tilt.png", outdir));
    printf("\nSaved halffield_tilt.{pdf,png}\n");
}
