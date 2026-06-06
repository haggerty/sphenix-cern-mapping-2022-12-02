// analysis/findCenter.C
// Locate the solenoid magnetic centre (Bz maximum on axis / Br sign change).
//
// Run from the repository root:
//   root -l -b -q 'analysis/findCenter.C+("data")'

#include "../sPHENIXFieldMap.cxx"
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
    const double dz   = fmap.GetDZ();

    // Bz maximum on axis (ir=0)
    double BzMax = -1e9; int izBzMax = -1;
    for (int iz = 0; iz < NZ; ++iz) {
        double v = fmap.GetBzGrid(0, iz);
        if (v > BzMax) { BzMax = v; izBzMax = iz; }
    }
    double zBzMax = zMin + izBzMax * dz;
    printf("\nBz maximum on axis: Bz=%.5f T at iz=%d  z=%.1f mm\n",
           BzMax, izBzMax, zBzMax);

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
