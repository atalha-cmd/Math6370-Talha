###################################################################
#  Makefile for Lab 3
#
#  Md Abu Talha
#  SMU Mathematics
#  Math 6370
###################################################################

all : driver1.exe driver2.exe

driver1.exe : driver1.cpp vectors.cpp Makefile
	g++ -O2 driver1.cpp vectors.cpp -p -o driver1.exe

driver2.exe : driver2.cpp tridiag_matvec.cpp Makefile
	g++ -O0 driver2.cpp tridiag_matvec.cpp -g -o driver2.exe

####### End of Makefile #######
