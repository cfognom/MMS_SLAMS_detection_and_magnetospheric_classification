# MMS SLAMS detection and magnetospheric classification

This repository holds the matlab implementation of the SLAMS detection- and magnetospheric classification algorithm of my master thesis. Consider reading it to get a better understanding of this repsitory.

# Prerequisites
* Matlab R2019b
* The IRFU-Matlab package on your Matlab path, available at https://github.com/irfu/irfu-matlab (Cloned on 21 Feb 2020)

Newer versions of Matlab and the IRFU package will probably work too. **Please note** that you do not need any of this if you are just interested in the raw SLAMS database.

# File structure

1. **/SLAMS_database:** holds the SLAMS database and is the default directory for new databases.

2. **/cache:** some intermediate data is cached here to avoid time consuming recalculations.

3. **/events:** contains csv-files with 'sitl'- and 'burst_segment'-data about SLAMS and the quasi-parallel bowshock from the MMS science data center. It was used to find suitable intervals for the validation dataset.

4. **/labeled:** holds the manually created validation dataset in two folders: **/data** and **/plot**. **/data** contains txt-files (one file for each validation interval) with alternating SLAMS start- and stop time data. Each file is named according to this naming convention: [id]\_[interval start time]\_[interval stop time].txt (time is specified as nanoseconds since 2000). **/plot** contains images of the intervals where SLAMS are marked for reference.

5. **/tools:** is essentially a toolbox that holds many useful functions that are used frequently.

6. **SLAMS_analysis:** 

7. **SLAMS_finder:**

8. **SLAMS_marker**

9. **SLAMS_results**

10. **find_SLAMS**

11. **load_SLAMS**

12. **validate_SLAMS_finder**
