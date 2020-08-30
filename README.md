# MMS SLAMS detection and magnetospheric classification

This repository holds the matlab implementation of the SLAMS detection- and magnetospheric classification algorithm of my master thesis. Consider reading it to get a better understanding of this repsitory.

## Prerequisites
* Matlab R2019b
* The IRFU-Matlab package on your Matlab path, available at https://github.com/irfu/irfu-matlab (Cloned on 21 Feb 2020)
* A local MMS database with

Newer versions of Matlab and the IRFU package will probably work too. **Please note** that you do not need any of this if you are just interested in the raw SLAMS database. This is only if you want to continue development, create a new SLAMS database with other settings or use the database loader, [load_SLAMS.m](load_SLAMS.m)

## File structure

1. [**SLAMS_database/**](SLAMS_database) Holds the SLAMS database and is the default directory for new databases.

2. [**cache/**](cache) Some intermediate data is cached here to avoid time consuming recalculations.

3. [**events/**](events) Contains csv-files with ``sitl``- and ``burst_segment``-data about SLAMS and the quasi-parallel bowshock from the MMS science data center. It was used to find suitable intervals for the validation dataset.

4. [**labeled/**](labeled) Holds the manually created validation dataset in two folders: ``data/`` and ``plot/``. ``data/`` contains txt-files (one file for each validation interval) with alternating SLAMS start- and stop time data. Each file is named according to this naming convention: ``[id]_[interval start time]_[interval stop time].txt`` (time is specified as nanoseconds since 2000). ``plot/`` contains images of the intervals where SLAMS are marked for reference.

5. [**tools/**](tools) This is essentially a toolbox that holds many useful functions.

6. [**SLAMS_analysis.m**](SLAMS_analysis.m) The function used to investigate SLAMS propteries (Section 3.6 and 4.2 of the thesis). Requires a SLAMS-database struct as input, which can be produced by ``load_SLAMS.m``.

7. [**SLAMS_finder.m**](SLAMS_finder.m) The main file which holds the ``SLAMS_finder``-class. The class implements both the SLAMS finder algorithm and the magnetospheric classifier. To find SLAMS, first initialize the class, then run the ``evaluate`` method and provide a time interval. ``evaluate`` then outputs a struct of SLAMS. Usage example can be found in ``find_SLAMS.m``.

8. [**SLAMS_marker.m**](SLAMS_marker.m) Contains the code used to manually mark SLAMS to create the validation dataset.

9. [**SLAMS_results.m**](SLAMS_results.m) The function used to produce the results in section 4.2 of the thesis. Requires a SLAMS-database struct as input, which can be produced by ``load_SLAMS.m``.

10. [**find_SLAMS.m**](find_SLAMS.m) This code can produce new SLAMS-databases by searcing multiple intervals. **Before running it**, ensure that the ``SLAMS_db_path`` variable is the location of the SLAMS databse on your machine.

11. [**load_SLAMS.m**](load_SLAMS.m) A function that takes the name of a SLAMS-databse (the name of the folder where the database is located) and loads it into a struct. **Before running it**, ensure that the ``SLAMS_db_path`` variable is the location of the SLAMS databse on your machine.

12. [**validate_SLAMS_finder.m**](validate_SLAMS_finder.m) The code used to validate the SLAMS detection algorithm against the validation set (section 3.2 of the thesis).

## SLAMS database structure

...

## Authors

* **Carl Foghammar Nömtak** - [cfognom](https://github.com/cfognom)

## License

This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details
