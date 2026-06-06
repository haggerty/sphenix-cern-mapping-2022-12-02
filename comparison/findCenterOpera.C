// comparison/findCenterOpera.C
// Locate the OPERA solenoid magnetic centre on axis, using the same robust
// estimators as analysis/findCenter.C (symmetry centre + parabola vertex +
// off-axis Br zero crossing).  The on-axis Bz is flat to < 0.1 mT over several
// cm, so a plain argmax on the 2 cm grid is unreliable (it reports +40 mm even
// though the map is built with the 28.5 mm design coil offset).
//
// Run from the repository root:
//   root -l -b -q 'comparison/findCenterOpera.C+("data/sphenix3dmapxyz.root")'

#include "OperaMap.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

void findCenterOpera(const char* operaFile = "data/sphenix3dmapxyz.root")
{
    OperaMap op;
    if (!op.load(operaFile)) return;

    // on-axis Bz(z); off-axis Br = bx along +x at radius r.  z, r in mm.
    auto BzAxis = [&](double z){ double bx,by,bz; op.get(0.0,0.0,z/10.0,bx,by,bz); return bz; };
    auto BrAt   = [&](double r,double z){ double bx,by,bz; op.get(r/10.0,0.0,z/10.0,bx,by,bz); return bx; };

    const double zMin = op.z0*10.0;                       // mm
    const double zMax = (op.z0 + (op.nz-1)*op.dz)*10.0;   // mm
    const double dz   = op.dz*10.0;                       // mm (grid step)

    // (1) Naive argmax on grid nodes (grid-limited; reference only)
    double BzMax=-1e9, zArg=zMin;
    for (double z=zMin; z<=zMax+1e-6; z+=dz) { double v=BzAxis(z); if (v>BzMax){BzMax=v; zArg=z;} }
    printf("\nNaive Bz argmax on axis (grid-limited): Bz=%.5f T at z=%.1f mm\n", BzMax, zArg);

    // (2) Symmetry centre: minimise S(c) = Sum_d [Bz(c+d) - Bz(c-d)]^2
    const double dMax=600.0, dStep=20.0;   // mm
    double cSym=zArg, sBest=1e30;
    for (double c=zArg-200.0; c<=zArg+200.0; c+=1.0) {
        double s=0.0; int n=0;
        for (double d=dStep; d<=dMax; d+=dStep) {
            if (c-d<zMin || c+d>zMax) continue;
            double diff=BzAxis(c+d)-BzAxis(c-d); s+=diff*diff; ++n;
        }
        if (n>0 && s<sBest) { sBest=s; cSym=c; }
    }
    printf("On-axis centre (Bz symmetry, |d|<=%.0f mm):     z = %.1f mm\n", dMax, cSym);

    // (3) Parabola vertex: least-squares Bz = a0 + a1 z + a2 z^2; centre = -a1/(2 a2)
    {
        const double W=600.0, step=20.0;
        double S0=0,S1=0,S2=0,S3=0,S4=0, T0=0,T1=0,T2=0;
        for (double z=cSym-W; z<=cSym+W; z+=step) {
            if (z<zMin || z>zMax) continue;
            double y=BzAxis(z), z2=z*z;
            S0+=1; S1+=z; S2+=z2; S3+=z2*z; S4+=z2*z2; T0+=y; T1+=y*z; T2+=y*z2;
        }
        double m[3][4]={{S0,S1,S2,T0},{S1,S2,S3,T1},{S2,S3,S4,T2}};
        for (int col=0; col<3; ++col) {
            int piv=col;
            for (int r=col+1;r<3;++r) if (std::abs(m[r][col])>std::abs(m[piv][col])) piv=r;
            for (int j=0;j<4;++j) std::swap(m[col][j],m[piv][j]);
            double d=m[col][col]; for (int j=col;j<4;++j) m[col][j]/=d;
            for (int r=0;r<3;++r) if (r!=col){ double f=m[r][col]; for (int j=col;j<4;++j) m[r][j]-=f*m[col][j]; }
        }
        double a1=m[1][3], a2=m[2][3], cPar=(a2!=0.0)?-a1/(2.0*a2):0.0;
        printf("On-axis centre (parabola vertex, |z-c|<=%.0f mm): z = %.1f mm\n", W, cPar);
    }

    // (4) Off-axis Br (=bx) sign change at r=100 mm, central region
    printf("\nBr zero crossing at r=100 mm:\n");
    const double rTest=100.0;
    double prev=BrAt(rTest,zMin), prevz=zMin;
    for (double z=zMin+dz; z<=zMax+1e-6; z+=dz) {
        double cur=BrAt(rTest,z);
        if (prev*cur<=0.0 && (prev!=0.0||cur!=0.0) && std::abs(z)<400.0) {
            double zc=prevz - prev*(z-prevz)/(cur-prev);
            printf("  z=%.1f..%.1f mm  Br=%+.6f..%+.6f T  -> zero at z=%.2f mm\n",
                   prevz, z, prev, cur, zc);
        }
        prev=cur; prevz=z;
    }
}
