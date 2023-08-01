## Stopover Site Segmentation

MoveApps

Github repository: *github.com/movestore/stopovers*

## Description
Segmentation alogrithm that efficiently extracts or annotates segments of each individual track that fulfill the criterion to stay within a disk of given radius for a minimum duration of time. This algorithm was often used to extract stopover sites of migrating animals.

## Documentation
The App is based on a segementation algorithm that has been developed and published in collaboration with Maike Buchin and Stef Sijben from the Univeristy of Bochum. (https://link.springer.com/chapter/10.1007/978-3-642-32316-4_2).

The algorithm passes through each individual movement track and calculates segments that fullfil the criteria to lie in a given maximum radius for a given minimum time duration. The calculations are based on great-circle distances of lat/long locations and relate each location to a potential smallest disk midpoint. The smallest disk calculations are implemented acoding to an algorithm described in Section 4.7 of De Berg, Cheong, Van Kreveld, Overmars: Computational Geometry, Algorithms and Applications, 3rd ed. Springer, 2008. The several computational optimisation steps make this App highly efficient.

The algorithm is based on great circle distance calculations if the loctions are WGS84 (longitute/latitude), else it assumes carthesian coordinates (planar). Note the option to annotate your tracks only or filter for locations in stopovers only with generating one track per stopover (see Setting 'annot').


### Input data
move2 location object

### Output data
move2 location object, either the input with 'stopover' indicated or only the locations within any stopover

### Artefacts
`stopover_sites.csv`: .csv-file with Table of all individuals' stopover sites and their properties (arrival, departure, mid location, stopover duration and radius)

### Settings 
**Minimum stopover duration (`duration`):** Defined minimum duration that the animal has to stay in a given radius for it to be considered a stopover site. Unit: `seconds`.

**Maximum stopover radius (`radius`):** Defined maximum radius the animal has to stay in for a given duration of time for it to be considered stopover site. Unit: `metres`.

**Annotate full data set with 'stopover' (`annot`):** Select if you want the full input data set returned, annotated with an attribute 'stopover' (entries 'move' or 'stopover) or a move2 location object returned that contains one track per stopover and all locations of the stopover included (default). Movement location will be dropped in that case.

### Null or error handling:
**Setting `duration`:** If no duration AND no radius are given, the input data set is returned with a warning. If no duraiton is given (NULL), but a radius is defined then a default duration of 86400 seconds = 24 hours = 1 day is set. 

**Setting `radius`:** If no radius AND no duration are given, the input data set is returned with a warning. If no radius is given (NULL), but a duration is defined then a default radius of 10000 m = 10 km is set. 

**Data:** If there are no resting locations retained after all analyses, NULL is returned, likely leading to an error.
