# crystalfdps

A working  sample code using FDPS, writen entirely in 
Crystal language.

To cmpile, after you download and place the files (either from
zip file or by git clone), 
```
   shards install
   make fdpscr
```
should create an executable. Try
```
  ./fdpscr -h
```
to get the list of command-line options.

Edit Makefile to use OpenMP (enabled by default) and MPI

```
  ./fdpscr -n 100 > /dev/null
```
should give

```
FDPS on Crystal test code
     //==================================\\
     ||                                  ||
     || ::::::: ::::::. ::::::. .::::::. ||
     || ::      ::    : ::    : ::       ||
     || ::::::  ::    : ::::::'  `:::::. ||
     || ::      ::::::' ::      `......' ||
     ||     Framework for Developing     ||
     ||        Particle Simulator        ||
     ||     Version 5.0g (2019/09)       ||
     \\==================================//

       Home   : https://github.com/fdps/fdps 
       E-mail : fdps-support@mail.jmlab.jp
       Licence: MIT (see, https://github.com/FDPS/FDPS/blob/master/LICENSE)
       Note   : Please cite the following papers.
                - Iwasawa et al. (2016, Publications of the Astronomical Society of Japan, 68, 54)
                - Namekata et al. (2018, Publications of the Astronomical Society of Japan, 70, 70)

       Copyright (C) 2015 
         Masaki Iwasawa, Ataru Tanikawa, Natsuki Hosono,
         Keigo Nitadori, Takayuki Muranushi, Daisuke Namekata,
         Kentaro Nomura, Junichiro Makino and many others
******** FDPS has successfully begun. ********
time: 0.0, energy error: -0.0
time: 0.25, energy error: -6.676560796472655e-6
time: 0.5, energy error: -7.65611452148078e-6
time: 0.75, energy error: 3.913752066804322e-5
time: 1.0, energy error: 1.745688979141589e-5
```

