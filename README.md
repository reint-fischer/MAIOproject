# MAIOproject
MAIO Ocean Turbulence project with Giovanni Negro

To study ocean turbulence with drifters we are going to make, analyse and interpret observations.
To facilitate the collaborative process, we will use this repository to structure the different tasks.

The repository consists of a project description folder with background information, a data folder in which all intermediate data is stored, a figures folder in which results are stored and the code files used for the analysis.

The code consists of a number of files which perform different tasks performed by functions.
1. selectmedsea: code to select trajectories in the Mediterrenean and discard all other data
    a. OpenDrifterdata(filelocation): function to open NOAA netcdf file with drifter trajectories
    b. selectMedsea(ID,lat,lon,time): function to select trajectories between 30N and 45N and 5W and 36E
    c. Conversion(D): function to convert drifter locations from lon-lat to x-y in km.
    
2. findpairs: code to find and save the chance pairs in the Mediterrenean
    a. SelectPair(D): function that takes the Mediterrenean trajectories and computes pairs of IDs that come within 10 km in a 6 hour window
    - Unpair: remove the double pairs
    
3. computedistance: code to compute distance timeseries for the chance pairs
    a. ComputeDistance(ID1,ID2,Data_Mediterrenean): function to compute the distance between the pairs, find the minimum and return forward and backward timeseries from that minimum.
    
4. plottraj: code that plots chance pair trajectories and the dispersion distance.

