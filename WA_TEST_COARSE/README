# West Australia TC test case 
Here you can find all needed files to run this **TC Veronica WA** example.

To make vizualization easy I put 2 python scripts 
(1) to plot stations (plot_station.py) and 
(2) for spatial variables (plot_maps.py) for different variables 
    (i.e. create figs like this animation)

## Flowchart:
(1) Run the model for 7 days (2019/03/18 - 2019/03/25) and save history.nc, 
station.nc, hotfile_out.nc (use wwminput_forward.nml, just copy it into wwminput.nml)
(2) make plots for a few stations (i.e. 0, 1 and 5) 
    using: python3 plot_station.py station.nc --station=0 --out_path=./ (see --help). 
    You can make spatial maps for HS if want.
Don't forget to copy history.nc and station.nc in separate folder (./forward_save) for later reference.

(3) run WWM again but now **for only 4 days** and starting **3 days later** (2019/03/21) 
    using previously created hotstart_out.nc (step 1) from record 4.
    (keep boundary and wind forcing sections the same, model should find correct record in time). 
    This is setup in wwminput_restart.nml, so just copy that over old wwminput.nml 

(4) make plots again and compare results for period when they overlap. Hotstart should make that the same.


Ivica
