# limCode2020

From Clara Chung's SURP 2020 project with R. Bond and Patrick Breysse.
Original code before any changes made by Clara can be obtained from Patrick.

For examples of how to use lim and limlam, see Example_Calcs.ipynb by Patrick. 

You can find George's COMAP simulation catalogue at https://www.cita.utoronto.ca/~gstein/data/CO/COMAP_fullvolume/midz/z4.77-6.21_700Mpc/npz/ or from the CITA directory /mnt/raid-cita/echung/surp2020/lim_Clara/limlam_mocker/catalogues/COMAP_z4.77-6.21_700Mpc.
You can get the .h5 catalogue from Martine, or find it on CITA server in the directory /mnt/raid-cita/echung/surp2020/lim_Clara/limlam_mocker/catalogues/galaxy_catalogue.h5.

When loading peakpatch catalogue for the first time (load_halos.py -> load_peakpatch_catalogue()), it may be preferable to set saveHalos=True and save the halos to limCode2020/outputs/. If nothing is saved to /outputs/ and saveHalos=False is run, it will throw an error saying you have to save halos before loading them from /outputs/. Nothing is saved to /outputs/ initially from github because the files are too large. Same with the .h5 catalogue obtained from Martine. You can access the already saved /outputs/ on CITA server at /mnt/raid-cita/echung/surp2020/lim_Clara/outputs/.
