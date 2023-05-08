c***************************************************************************
       subroutine getch (tt)

c	NOV 98	MALV 		First Version. From the work:
c                               Rodrigo et al, Atmosfera, 3, pp.31-43, 19990
c       apr 02  FGG             New reactions
c       jul 03  MAC             Updated to agree with lower atmosphere
c***************************************************************************
	
      implicit none

      include 'param.h'
      include 'param_v3.h'

c Arguments 	
      real*8 		tt		! Temperature 


c**************************************************************************

ch2: h + o2 + co2 --> ho2 + co2
      ch2 = 1.3d-31 * exp(tt / 300.) ** (-1.6)

ch3: o + ho2 --> oh + o2
      ch3 = 3.0d-11 * exp(200.d0 / tt)

ch4: co + oh --> co2 + h
      ch4 = 1.57d-13 + 3.54d-33

ch5: ho2 + ho2 --> h2o2 + o2
      ch5 = 1.5d-12 * exp(19.d0 / tt)

ch7: oh + ho2 --> h2o + o2
      ch7 = 4.8d-11 * exp(250.d0 / tt)

ch9: o(1d) + h2o --> 2oh
      ch9 = 2.2d-10

ch10: o + o + co2 --> o2 + co2
      ch10 = 2.5 * 5.2d-35 * exp(900. / tt) 

ch11: o + oh --> o2 + h
      ch11 = 2.2d-11 * exp(120.d0 / tt)

ch13: h + ho2 --> h2 + o2
      ch13 = 6.48d-12

ch14: o(1d) + h2 --> h + oh
      ch14 = 1.1d-10

ch15: oh + h2 --> h + h2o
      ch15 = 5.5d-12 * exp (-2000.d0 / tt)

ch18: oh + h2o2 --> h2o + ho2
      ch18 = 2.9d-12 * exp (-160.d0 / tt)

ch19: o(1d) + co2 --> o + co2
      ch19 = 7.4d-11 * exp(120.d0 / tt)

ch20: o(1d) + o2 --> o + o2
      ch20 = 3.2d-11 * exp (70.d0 / tt)

ch21: o + o2 + co2 --> o3 + co2
      ch21 = 2.5 * 6.0d-34 * (tt / 300.) ** (-2.4) 

      return 
      end









