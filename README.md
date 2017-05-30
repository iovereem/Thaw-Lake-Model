# Thaw-Lake-Model
1-D numerical model of permafrost and subsidence processes
Originally written by Nora Matell
Now maintained and to be translated into Python by Irina Overeem

ThawLake1D Model couples a permafrost thermal model, a lake ice model and a subsidence model. It models heat conductio through a lake-permafrost system, evaluating temperature with depth.

There are three components to this model: 
- a permafrost thermal model, it calculates thermal condutivity and specific heat capacity
- lake ice model calculated conduction of heat within water depending on incoming solar radiation 
- a subsidence model: melted permafrost cells with excess ice are subsided by the excess ice content. 
The combinded processes allow for a talik to form under a lake. 

The model is 1D, the uppermost 10m has a high resolution (0.05m), below 10 m to 500m depth gridcells are set to be 1m.
Simulations run at daily timestep. Total simulation times are typicaly over 100's of years.

Model is designed and calibrated for Alaska Coastal Plain.
We calibrated the model against temperature data in the subsurface from the Drew point, AK, USGS meteorological station.
