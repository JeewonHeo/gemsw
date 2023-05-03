import os

run_eras = {
        8149: "Run2022D",
        8220: "Run2022E",
        8456: "Run2022G",
        }

runs = {
        8149: [357802, 357803, 357804, 357805, 357806, 357807, 357808, 357809, 357812, 357813, 357814, 357815],
        8220: [359691, 359693, 359694],
        8456: [362433, 362435, 362437, 362438, 362439],
        }

for fill in run_eras:
    # if fill == 8220: continue
    if fill == 8149: continue
    if fill == 8456: continue
    for run_number in runs[fill]:
        run_era = run_eras[fill]
        command = f"crab submit -c crabConfig_background.py "
        command += f"General.requestName={run_era}_{run_number}_05Apr2023_BKG_analysis "
        command += f"Data.inputDataset=/ZeroBias/{run_era}-v1/RAW "
        command += f"Data.runRange={run_number} "
        command += f"Data.outputDatasetTag={fill}_{run_number}_200 "
        print(command)
        os.system(command)
