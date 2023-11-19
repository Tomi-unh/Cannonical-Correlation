# Cannonical Correlation

  # Energetic Neutral Atoms Temperature Maps.

 Understanding the dynamics of the magnetotail both spatially and temporally is one of the biggest areas of research in space physics. How the energy transports from the tail side connection of the magnetic field all the way into the ionoshpere is crucial to understand this dynamics. Many challenges stand in the way of understand these various phenomena, not the least of which is the over reliance on in situ observations. Much of the progress made over the decades to resolve these issues revolves around having the various satellites being in the right place at the right time. This limitation slows down how much of the tail can be analyzed at a given time becuase almost all observations are localized and we rarely ever get a global context of these phenomena. Multi-point obeservations can offer some global view of the magnetotail but the coordination required to get a global context of an event is a challenge. Simulations have been used to bridge this gap, but that also relies on in situ measurements both for inputs into the algorithm and also for validation of the results.       
 Having a global view of the magnetotial can give another avenue of observing and analyzing the magnetotail and putting in situ measurements into a global context. This is where a temeprature map of the magnetotail comes into play [Keesee et. al, 2008 (10.1029/2008JA013130)]. This technique makes use of energetic neutral atoms (ENA) to image the magnetotial. ENAs are created from charge exchange collisions of energetic ions with cold neutrals. The collision causes the eneretic ions to gain an electron from the cold neutral and become a neutral itself. Little loss of information is lost in this process due to the energy difference between the cold neutral and the energetic ion. The little loss of information ensures that the newly formed energetic neutral still keeps the properties (majority of it at least) of the parent ion population. Being neutrals, they can be measured remotely from distant satellites like the Two Wide-Angle Imaging Neutral-Atom Spectrometers (TWINS) and used as a proxy for getting information about the parent ion population. The measured ENAs can then be projected back down into the equatorial plane as a function of several properties of the parent ions, in this case, the temperature. This process creates an ENA temeperature map of the magnetotail. The full details of the techniques can be seen in the paper linked above. 

# Cannonical Correlation Analysis.

Having the temperautre maps of the magnetotail, it is important to compare the dynamincs seen in the maps with ionoshperic phenomena to fully understand the energy transfer between these two coupled regions. How does the temporal and spatial features seen in the maps correlate with what is obersed in the ionosphere? This is the main question this project seek to answer, and one of the ways to do this is through to use of cannonical correlation analysis (CCA). Cannonical correlation is a powerful technque used to understand the relationships between multivaraite datasets with a  latent variable shared between them. This process finds the linear combination of the two sets of variables that miaximizes the correlation between them. With CCA, the correlation between two different magnetometer stations and their field data can be determined, with goemagnetic field activities being the latent variable that drives the variations that both stations observe. Following the technique detailed by Dods et. al, 2015 (doi:10.1002/2015JA021456), a network of the SuperMAG stations can be generated. The network essentially shows a map of the time variation correlations between the stations, where the stations are the nodes of the network, and the edges are the correlation. An edge is formed if and only if there exists a correlation between two station that exceeds a minimum threshold. Correlation that falls below the defined threshold represents the lack of connection between the nodes (station). 

# How to Run:
  
  1.) First thing is to get the actual correlation coefficient of all the stations with one another.     
  Run the main function in Cannonical_CrossCorrelation.py script twice.   
  Once for the whole month surrounding the time of interest (set the Duration arg to 28 in the main function)  
  Another for the day surrounding the time of interest (set the Duration to 1). This gets the correlation coefficients.    
  2.) Run the month_analysis function from the NetworkAnalyzer class in the Network_Analysis.py file. This gets the station threshold matrix  
  3.) Run the main function from the NetworkAnalyzer class to get the network diagrams, along with binned degree diagrams and Degree Centrality (Top 3 Stations).

  




# Sample Results
The network created helps to understand the spatial and temporal correlations between the temperature maps in the magnetotail and the ionospheric dynamics. 
