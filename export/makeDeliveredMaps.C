// makeDeliveredMaps.C
//
// Build the two field maps delivered to tracking, both on the production
// 111^3 / +-110 cm Cartesian grid used by PHField3DCartesian, branch layout
// x:y:z:bx:by:bz:hz (cm, T), identical to the OPERA tracking map.
//
// Both maps carry the two physical results of the 2022-12-02 mapping — a
// uniform +0.9 % amplitude scale and a 4.4 mrad horizontal yaw (toward -x,
// "east"; theta_y ~ 0, i.e. the magnet is level so only the horizontal yaw is
// applied).  They differ only in the BASE field:
//
//   (2) opera_matched_to_mapping : the OPERA *calculation*, rescaled and yawed
//       to match the measurement.  Defined everywhere in the cube (no cutoff).
//       NOT a measured map — it is OPERA corrected in the spirit of the mapping.
//
//   (3) measured_smoothed        : the *measured* field — the phi-averaged,
//       div.B-enforced (r,z) reconstruction of the point cloud (which removes
//       the ~10 mT per-point transverse noise but keeps the real amplitude and
//       z-profile), with the measured 4.4 mrad yaw reinserted as a rigid
//       rotation.  Honestly cut off (zero) beyond the measured coverage r=90 cm.
//
// A rigid yaw R that maps z_hat -> (theta_x, theta_y, 1) acts on a field map as
// B'(r) = R . B(R^{-1} r); to first order:
//   sample at  r_src = (x - tx*z, y - ty*z, z + tx*x + ty*y)
//   rotate     Bx'=bx+tx*bz, By'=by+ty*bz, Bz'=bz-tx*bx-ty*by
//
// The yaw, amplitude scale and vertical tilt are ARGUMENTS, so the same code
// produces the menu of maps that bracket the (degenerate) tilt systematic:
//
//   root -l -b -q 'export/makeDeliveredMaps.C+(-4.40)'  # raw apparent yaw (now disfavoured)
//   root -l -b -q 'export/makeDeliveredMaps.C+(-2.20)'  # field-proportional yaw (half-field analysis)
//   root -l -b -q 'export/makeDeliveredMaps.C+(0.0)'    # no tilt (amplitude scale only)
//
// Signature: makeDeliveredMaps(yawMrad=-4.40, scale=1.0091, thetaYmrad=0.0, tag="").
// Output filenames carry a yaw tag (default "yaw<|yawMrad|>mrad") so the variants
// coexist; the yaw=0 map replaces the old "run OPERA with magfield_rescale"
// recipe with an explicit scale-only file.  Defaults reproduce the original
// 4.4 mrad / +0.9 % corrections.
//
// Run from the project directory.

#include "../sPHENIXFieldMap.cxx"
#include "../comparison/OperaMap.h"

#include "TFile.h"
#include "TNtuple.h"
#include "TSystem.h"
#include "TString.h"

#include <cmath>
#include <cstdio>
#include <string>

// ── Default physical corrections from the 2022-12-02 mapping (all overridable) ─
//   scale   1.0091 : +0.9 % amplitude (best-fit, == Bz(0,0,0) ratio)
//   yawMrad -4.40  : horizontal yaw toward -x ("east").  The half-field analysis
//                    splits this into ~2.2 mrad field-proportional + a ~1.4 mT
//                    current-independent offset, so -2.20 and 0.0 are the other
//                    defensible choices (see README half-field / global sections).
//   thetaY   0.0   : magnet level (theta_y ~ 0): no vertical tilt applied

// ── Output grid (production PHField3DCartesian grid) ─────────────────────────
static const int   N     = 111;            // nodes/axis
static const float GMIN  = -110.f;         // cm
static const float GMAX  =  110.f;         // cm
static const float DG    =    2.f;         // cm/step
static const float CM2MM =   10.f;

static inline void rotPos(double x,double y,double z,double tx,double ty,
                          double&xs,double&ys,double&zs){
    xs = x - tx*z;  ys = y - ty*z;  zs = z + tx*x + ty*y;
}
static inline void rotVec(double bx,double by,double bz,double tx,double ty,
                          double&Bx,double&By,double&Bz){
    Bx = bx + tx*bz;  By = by + ty*bz;  Bz = bz - tx*bx - ty*by;
}
static inline double clampg(double v){ return v<GMIN?GMIN:(v>GMAX?GMAX:v); }

// Write one 111^3 TNtuple given a per-point sampler f(x,y,z)->(Bx,By,Bz) [cm,T].
template<class F>
static void writeMap(const char* outFile, const char* title, F sample)
{
    gSystem->mkdir("output", kTRUE);
    TFile *fOut = TFile::Open(outFile, "RECREATE");
    if (!fOut || fOut->IsZombie()) { fprintf(stderr,"ERROR: cannot create %s\n",outFile); return; }
    TNtuple *nt = new TNtuple("fieldmap", title, "x:y:z:bx:by:bz:hz");
    long nZero=0; float v[7];
    double bz0=0;
    for (int ix=0; ix<N; ++ix){ float x=GMIN+ix*DG;
        for (int iy=0; iy<N; ++iy){ float y=GMIN+iy*DG;
            for (int iz=0; iz<N; ++iz){ float z=GMIN+iz*DG;
                double Bx,By,Bz; sample(x,y,z,Bx,By,Bz);
                v[0]=x; v[1]=y; v[2]=z; v[3]=(float)Bx; v[4]=(float)By; v[5]=(float)Bz; v[6]=(float)Bz;
                nt->Fill(v);
                if (Bx==0&&By==0&&Bz==0) ++nZero;
                if (x==0&&y==0&&z==0) bz0=Bz;
            }}}
    fOut->cd(); nt->Write(); fOut->Close();
    printf("  -> %s\n     %d entries, Bz(0,0,0)=%.5f T, %ld zero-field points\n",
           outFile, N*N*N, bz0, nZero);
}

// The three delivered cases that bracket the (degenerate) tilt systematic.
// All share the +0.9 % amplitude; they differ only in the applied yaw:
//
//   case          yawMrad   meaning
//   ------------  -------   ---------------------------------------------------
//   no-tilt          0.0    amplitude scale only; consistent with the survey
//                           (~0.17 mrad) and with the tilt being a measurement
//                           systematic.  Lower bound of the menu.
//   field-prop      -2.20   the field-PROPORTIONAL yaw from the half-field
//                           analysis (slope of <B_perp> vs Bz), with the
//                           ~1.4 mT current-independent offset removed.  Best
//                           single estimate if a real rotation is assumed.
//   raw apparent    -4.40   the full-field apparent yaw (<B_perp>/Bz).  Now
//                           disfavoured: ~1 mrad of it is the additive offset,
//                           and the remainder still exceeds the survey.
//
// Direction: yaw is toward -x ("east"), theta_y ~ 0 (magnet level).  See the
// README "half-field" and "global" sections for the decomposition and errors.

void makeDeliveredMaps(double yawMrad = -4.40, double scale = 1.0091,
                       double thetaYmrad = 0.0, const char* tag = "")
{
    const double THETA_X = yawMrad    * 1e-3;
    const double THETA_Y = thetaYmrad * 1e-3;
    const double SCALE   = scale;
    const std::string sfx = (tag && tag[0]) ? std::string(tag)
                                            : Form("yaw%.1fmrad", std::fabs(yawMrad));
    const std::string fOpera =
        Form("output/sphenix_solenoid_opera_matched_to_mapping_2022-12-02_%s.root", sfx.c_str());
    const std::string fMeas =
        Form("output/sphenix_solenoid_measured_smoothed_2022-12-02_%s.root", sfx.c_str());

    const char *FINE  = "data/pointCloudFineFullField.csv";
    const char *ROUGH = "data/pointCloudRoughFullField.csv";
    const char *OPERA = "data/sphenix3dtrackingmapxyz.root";   // 111^3, +-110 cm

    printf("=== Loading measured field map ===\n");
    sPHENIXFieldMap meas(FINE, ROUGH);

    printf("=== Loading OPERA tracking map ===\n");
    OperaMap op;
    const bool haveOpera = op.load(OPERA);
    if (!haveOpera)
        printf("  WARNING: cannot load %s; skipping the opera_matched map.\n", OPERA);

    printf("=== Corrections: scale=%.4f, yaw theta_x=%.2f mrad, theta_y=%.2f mrad  [tag %s] ===\n",
           SCALE, THETA_X*1e3, THETA_Y*1e3, sfx.c_str());

    // (2) OPERA matched to the mapping: SCALE * R . B_opera(R^{-1} r).
    //     Sample position clamped to the grid edge so the full cube is defined.
    if (haveOpera) {
        printf("\n=== (2) opera_matched_to_mapping  [%s] ===\n", sfx.c_str());
        writeMap(fOpera.c_str(),
                 "sPHENIX solenoid: OPERA scaled+yawed to the 2022-12-02 mapping (Cartesian, cm/T)",
                 [&](double x,double y,double z,double&Bx,double&By,double&Bz){
                     double xs,ys,zs; rotPos(x,y,z,THETA_X,THETA_Y,xs,ys,zs);
                     double bx,by,bz; op.get(clampg(xs),clampg(ys),clampg(zs),bx,by,bz);
                     rotVec(bx,by,bz,THETA_X,THETA_Y,Bx,By,Bz);
                     Bx*=SCALE; By*=SCALE; Bz*=SCALE;
                 });
    }

    // (3) measured, smoothed: R . B_meas(R^{-1} r).  GetFieldXYZ returns 0 for
    //     r>900 mm, so the map is honestly cut off beyond the measured coverage.
    printf("\n=== (3) measured_smoothed  [%s] ===\n", sfx.c_str());
    writeMap(fMeas.c_str(),
             "sPHENIX solenoid: measured field, smoothed + measured yaw (Cartesian, cm/T)",
             [&](double x,double y,double z,double&Bx,double&By,double&Bz){
                 double xs,ys,zs; rotPos(x,y,z,THETA_X,THETA_Y,xs,ys,zs);
                 double bx,by,bz; meas.GetFieldXYZ(xs*CM2MM,ys*CM2MM,zs*CM2MM,bx,by,bz);
                 rotVec(bx,by,bz,THETA_X,THETA_Y,Bx,By,Bz);
             });

    // ── Verification: <Bx> on the |z|<=40 cm, r<=40 cm core should show the yaw ─
    //    (phi-average of Bx ~ theta_x*Bz for both maps).
    printf("\n=== Verification: <Bx> over the r<=40 cm, |z|<=40 cm core (expect ~ %.1f mT) ===\n",
           THETA_X*1.4*1e3);
    auto coreBx=[&](const std::string& p){
        TFile*f=TFile::Open(p.c_str()); if(!f||f->IsZombie()) return;
        TNtuple*nt=(TNtuple*)f->Get("fieldmap");
        float x,y,z,bx; nt->SetBranchAddress("x",&x); nt->SetBranchAddress("y",&y);
        nt->SetBranchAddress("z",&z); nt->SetBranchAddress("bx",&bx);
        double s=0; long n=0;
        for (Long64_t i=0;i<nt->GetEntries();++i){ nt->GetEntry(i);
            if (x*x+y*y<=40*40 && std::abs(z)<=40 && (x!=0||y!=0)){ s+=bx; ++n; } }
        printf("  %-72s <Bx> = %+5.2f mT\n", gSystem->BaseName(p.c_str()), n? s/n*1e3:0.);
        f->Close();
    };
    if (haveOpera) coreBx(fOpera);
    coreBx(fMeas);
    printf("\nDone.\n");
}
