// analysis/findCenter.C
// Locate the solenoid magnetic centre on axis.
//
// The on-axis Bz of a solenoid is flat to < 0.1 mT over several cm around the
// peak, so a plain argmax on the field grid cannot localise the centre (a 1e-6 T
// difference between two nodes flips the answer by a grid step).  This macro
// therefore reports three estimates:
//   1. naive argmax of on-axis Bz   (grid-limited; shown for reference only)
//   2. symmetry centre of Bz(z)     -- robust
//   3. parabola-vertex fit of Bz(z) -- robust cross-check
// plus the off-axis Br sign change (also robust).
//
// Run from the repository root:
//   root -l -b -q 'analysis/findCenter.C+("data")'

#include "../sPHENIXFieldMap.cxx"
#include <TString.h>
#include <cstdio>
#include <cmath>
#include <algorithm>

void findCenter(const char* datadir = "data")
{
    sPHENIXFieldMap fmap(
        Form("%s/pointCloudFineFullField.csv",  datadir),
        Form("%s/pointCloudRoughFullField.csv", datadir));

    const int    NZ   = fmap.GetNZ();
    const double zMin = fmap.GetZMin();
    const double zMax = fmap.GetZMax();
    const double dz   = fmap.GetDZ();

    // (1) Naive Bz maximum on axis (ir=0) -- grid-limited, unreliable for the
    //     flat top; kept only as a reference point for the table below.
    double BzMax = -1e9; int izBzMax = -1;
    for (int iz = 0; iz < NZ; ++iz) {
        double v = fmap.GetBzGrid(0, iz);
        if (v > BzMax) { BzMax = v; izBzMax = iz; }
    }
    double zBzMax = zMin + izBzMax * dz;
    printf("\nNaive Bz argmax on axis (grid-limited): Bz=%.5f T at iz=%d  z=%.1f mm\n",
           BzMax, izBzMax, zBzMax);

    // ── Robust on-axis centre (uses interpolation, so c is continuous) ─────────
    auto BzAxis = [&](double z) {
        double Br, Bphi, Bz; fmap.GetField(0.0, 0.0, z, Br, Bphi, Bz); return Bz;
    };

    // (2) Symmetry centre: minimise S(c) = Sum_d [Bz(c+d) - Bz(c-d)]^2 over a
    //     window where the central field has clear curvature.
    const double dMax = 600.0, dStep = 20.0;   // mm
    double cSym = zBzMax, sBest = 1e30;
    for (double c = zBzMax - 200.0; c <= zBzMax + 200.0; c += 1.0) {
        double s = 0.0; int npair = 0;
        for (double d = dStep; d <= dMax; d += dStep) {
            if (c - d < zMin || c + d > zMax) continue;
            double diff = BzAxis(c + d) - BzAxis(c - d);
            s += diff * diff; ++npair;
        }
        if (npair > 0 && s < sBest) { sBest = s; cSym = c; }
    }
    printf("On-axis centre (Bz symmetry, |d|<=%.0f mm):     z = %.1f mm\n",
           dMax, cSym);

    // (3) Parabola vertex: least-squares Bz = a0 + a1 z + a2 z^2 over the same
    //     window; centre = -a1 / (2 a2).
    {
        const double W = 600.0, step = 20.0;   // fit half-window, sampling step
        double S0=0,S1=0,S2=0,S3=0,S4=0, T0=0,T1=0,T2=0;
        for (double z = cSym - W; z <= cSym + W; z += step) {
            if (z < zMin || z > zMax) continue;
            double y = BzAxis(z), z2 = z*z;
            S0+=1; S1+=z; S2+=z2; S3+=z2*z; S4+=z2*z2;
            T0+=y; T1+=y*z; T2+=y*z2;
        }
        // Solve [S0 S1 S2; S1 S2 S3; S2 S3 S4] a = [T0; T1; T2] (Gauss-Jordan).
        double m[3][4] = {{S0,S1,S2,T0},{S1,S2,S3,T1},{S2,S3,S4,T2}};
        for (int col = 0; col < 3; ++col) {
            int piv = col;
            for (int r = col+1; r < 3; ++r)
                if (std::abs(m[r][col]) > std::abs(m[piv][col])) piv = r;
            for (int j = 0; j < 4; ++j) std::swap(m[col][j], m[piv][j]);
            double d = m[col][col];
            for (int j = col; j < 4; ++j) m[col][j] /= d;
            for (int r = 0; r < 3; ++r) if (r != col) {
                double fac = m[r][col];
                for (int j = col; j < 4; ++j) m[r][j] -= fac * m[col][j];
            }
        }
        double a1 = m[1][3], a2 = m[2][3];
        double cPar = (a2 != 0.0) ? -a1 / (2.0 * a2) : 0.0;
        printf("On-axis centre (parabola vertex, |z-c|<=%.0f mm): z = %.1f mm\n",
               W, cPar);
    }

    // Br sign change at r=100 mm (ir=4)
    int irTest = 4;
    printf("\nBr zero crossings at r=100 mm:\n");
    for (int iz = 1; iz < NZ; ++iz) {
        double b0 = fmap.GetBrGrid(irTest, iz-1);
        double b1 = fmap.GetBrGrid(irTest, iz  );
        if (b0 * b1 <= 0.0 && (b0 != 0.0 || b1 != 0.0)) {
            double z0 = zMin + (iz-1) * dz;
            double z1 = zMin +  iz    * dz;
            double zBrZero = z0 - b0 * (z1 - z0) / (b1 - b0);
            printf("  iz=%d..%d  z=%.1f..%.1f mm  Br=%.5f..%.5f T  -> zero at z=%.2f mm\n",
                   iz-1, iz, z0, z1, b0, b1, zBrZero);
        }
    }

    printf("\n  iz   z[mm]   Bz_axis[T]   Br(r=100mm)[T]\n");
    int i0 = std::max(0, izBzMax - 8);
    int i1 = std::min(NZ-1, izBzMax + 8);
    for (int iz = i0; iz <= i1; ++iz) {
        double z = zMin + iz * dz;
        printf("  %3d  %7.1f  %9.5f    %+10.6f\n",
               iz, z, fmap.GetBzGrid(0,iz), fmap.GetBrGrid(irTest,iz));
    }
}
