# intensitynet 1.2

* Added new functionality. Now the package can work with marked events (numerical, categorical or both). Each edge in the graph provided by the intensitynet object contains all the information of its associated events.  

* Changed internal function name 'AllEdgeIntensities()' to 'EdgeIntensitiesAndProportions()'  In order to adapt its significate to the new functionalities.

* Changed funcion name 'CalculateEventIntensities()' to 'RelateEventsToNetwork()'. In order to adapt its significate to the new functionalities.

* Greatly reduced the computational time to calculate event-related intensities with the function 'RelateEventsToNetwork()'



# intensitynet 1.1

* Modified DESCRIPTION name references 'IntensityNet' to 'intensitynet' to match the package name.

* Added option to select event the error distance. Now the user can determine the maximum distance from an event to an edge to be considered part of that edge.

* Added an option to the plots ('plot' and 'PlotHeatmap' functions) to also plot the events. The events will be displayed as orange squares.

* Fixed bug when using the 'plot()' function with intensitynet objects which showed properly the axis in a .pdf but not in a dynamic plot. 
