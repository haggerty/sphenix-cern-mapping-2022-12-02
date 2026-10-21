# sPHENIX solenoid cryostat surveys

Position and orientation of the solenoid cryostat from the survey reports in this
directory, expressed in sPHENIX coordinates.

## Conventions

- **sPHENIX coordinates:** x = west, y = up, z = north (beam axis). Units mm and mrad.
- **Survey coordinates** (all reports): x = west, y = south, z = up. Conversion to sPHENIX:
  `x = x_s`, `y = z_s`, `z = −y_s`.
- **Yaw = dx/dz** of the solenoid axis; + means the north end points west. For a
  right-handed rotation about +y (up), this equals the rotation angle.
- **Pitch = dy/dz** of the solenoid axis; + means the north end is up. A right-handed
  rotation about +x tips the north end down, so the right-handed pitch angle is −dy/dz.
  The slope is quoted because its sign reads directly.

## Results

| | 2022-08-16 | 2026-10-07 | change (2026 − 2022) |
|---|---|---|---|
| Report | [`2022-08-16_magnet_location.pdf`](2022-08-16_magnet_location.pdf) | [`2026-10-07_cryostat_location.pdf`](2026-10-07_cryostat_location.pdf) | |
| Reference frame | RHIC CONTROL::Local_08 oclockBH | RHIC CONTROL::Local_08 oclockBH (same) | |
| What was recorded | position and rotation of the "S-Phenix Magnet" frame, plus its 3×3 axis matrix (inches, degrees) | fitted cryostat cylinder (ends, centre, radius 1382.4 mm, length 3846 mm), end circles and end planes (mm) | |
| Centre x (mm) | +2.03 | +3.12 | +1.09 |
| Centre y (mm) | +1.46 | +0.94 | −0.52 |
| Centre z (mm) | −4.52 | −2.93 | +1.59 |
| **Yaw** dx/dz (mrad) | **+0.85** | **+1.42** | +0.57 |
| **Pitch** dy/dz (mrad) | **+0.17** | **+0.69** | +0.52 |

The centre moved by ~2 mm and both angles grew by ~0.5 mrad. The two surveys do not
define the axis the same way (a surveyor-defined magnet frame in 2022, a cylinder fit
in 2026), so part or all of that difference may come from the method rather than
movement. Within the 2026 report alone, the end-plane normals differ from the cylinder
axis by up to 0.4 mrad (below).

## How the numbers were obtained

**2022-08-16.** The rows of the report's axis matrix are the magnet frame's axes in the
working frame. The magnet Y axis (pointing south) is (−0.000847, 1.000000, −0.000172);
reversing it to point north and converting to sPHENIX gives yaw +0.847, pitch +0.172
mrad. The centre is the frame's translation, (0.079889, 0.177981, 0.057425) in, converted
to mm (×25.4) and to sPHENIX. *Assumption:* the frame origin is the cryostat centre —
how the surveyors defined the frame is not recorded in the report.

**2026-10-07.** The axis runs from the fitted cylinder's south end (0.392, 1925.924,
−0.387) mm to its north end (5.838, −1920.066, 2.258) mm; converted to sPHENIX this gives
yaw +1.416, pitch +0.688 mrad. The centre is the fitted cylinder centre (3.115, 2.929,
0.936) mm. Other features in the same report:

| 2026 feature | yaw (mrad) | pitch (mrad) | centre z (mm) |
|---|---:|---:|---:|
| cylinder fit (adopted) | +1.42 | +0.69 | −2.93 |
| end-circle centres | +1.41 | +0.61 | −3.41 |
| south end-plane normal | +1.27 | +0.41 | — |
| north end-plane normal | +1.00 | +0.34 | — |

**Sign check.** Each report also lists "projected angles" computed by the survey
software. The axis vectors above reproduce them exactly (2022 Y axis: 90.0485° and
−0.0098°; 2026 cylinder: −89.9189° and 179.9606°), which confirms the reports were read
the right way round.

## Other surveys

- **2022-08-12** ([`2022-08-12_magnet_location.pdf`](2022-08-12_magnet_location.pdf)) —
  the magnet at roll-in, in the *sPhenix Assembly Ref-2022* frame (not the RHIC frame;
  the magnet was ~1.53 m south of its final position). Yaw +1.40, pitch −0.06 mrad. Not
  directly comparable with the two surveys above.
- **2022-12-06** — 18 raw laser scans of the magnet, four days after the 2022-12-02 field
  mapping (`sPhenix-Mag-12-06-22_07-12-2022_10-59-50.zip`, not in this repo). No fitted
  results; they would need to be processed by the survey group.
