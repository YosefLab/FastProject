Multi platform bh_tsne

Van der Maaten's bh_tsne for various platforms (Windows, Linux and MAC) and architecture (32 or 64 bit) is statically linked against OpenBLAS(http://xianyi.github.io/OpenBLAS) library. Even though Mac OS X has its own cblas (from Xcode), OpenBLAS seems to be faster. The MatLab function fast_tsne automatically detects the operating system as well as the architecture and runs the corresponding executable. If you have trouble with the 64bit versions then just modify the file fast_tsne.m to run only the 32 bit versions. In both Windows (tested in Win 7 and 8) and Linux (tested in Suse, Ubuntu and Linux Mint) the 64 bit version is about 2 seconds faster than the corresponding 32 bit executable for 6000 samples. All the executable are compiled using gcc v4.x.x. The source.rar archive contains the source as well as the OpenBLAS static libraries for all the supported operating systems. It also has dos batch files to compile the executables under Windows and bash shell script to compile the program under non-windows operating systems like Linux, Mac or BSD.


Lagnajeet Pradhan
University of Texas at Dallas
Email : lxp090020@utdallas.edu
