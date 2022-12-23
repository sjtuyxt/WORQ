# WORQ Tree

## Dataset ##
| Name      |    \# of records | Download  |
| :-------- | :--------:| :-- |
| Gowalla  | 6.4M |  https://snap.stanford.edu/data/loc-gowalla.htmlb   |
| BerlinMOD     |   56.1M |  https://secondo-database.github.io/BerlinMOD/BerlinMOD.html  |
| Chicagotaxi      |    167M | https://data.cityofchicago.org/Transportation/Taxi-Trips/wrvz-psew  |
| LondonCourier     |    9.9M | https://researchdata.edu.au/courier-trajectories-aposecourier-datasetapos/939716  |

## Code ##
``` bash
g++ -g -Wall -std=c++11 {Bkd|Alstrup|LSMR|PH|WorqNoLazy|Worq}.cpp -o run
./run
```
The program will read "DB.txt" as input and read query range in "query.txt". 

The output will be saved in "result.txt".
