c-----------------------------------------------------------------------
c INCLUDE 'aerdust.h'
c Spectral  Radiative properties of the Martian dust
c SPECIAL "bright dust" (Clancy and Lee, 1991)


c visible (Clancy and Lee, 1991)
      integer ndustvis
      parameter(ndustvis=2)
      REAL longdustvis(ndustvis),epdustvis(ndustvis)
     &    ,omegdustvis(ndustvis),gdustvis(ndustvis)


c IR (composite: Toon 1977 (5-15micrometers) + Forget, 1996 (>15micrometers) 
      integer ndustir
      parameter(ndustir=124)
      REAL longdustir(ndustir),epdustir(ndustir)
     &    ,omegdustir(ndustir),gdustir(ndustir)

c Ratio solsir = Tau(vis) / Tau (9 microns)
      REAL solsir
      PARAMETER (solsir=2.E+0)

c-----------------------------------------------------------------------
       

      data longdustvis/
     & 0.210000E-06,0.414800E-05/
      data epdustvis/2.5,2.5/
      data omegdustvis/0.92,0.92/
      data gdustvis/0.55,0.55/
      data longdustir/
     & 0.556000E-05,0.625000E-05,0.714000E-05,0.725000E-05,0.735000E-05
     &,0.746000E-05,0.758000E-05,0.769000E-05,0.781000E-05,0.794000E-05
     &,0.800000E-05,0.807000E-05,0.813000E-05,0.820000E-05,0.826000E-05
     &,0.833000E-05,0.840000E-05,0.848000E-05,0.855000E-05,0.862000E-05
     &,0.870000E-05,0.877000E-05,0.880000E-05,0.884956E-05,0.892857E-05
     &,0.900901E-05,0.909091E-05,0.917431E-05,0.925926E-05,0.934579E-05
     &,0.943396E-05,0.952381E-05,0.961539E-05,0.970874E-05,0.980392E-05
     &,0.990099E-05,0.100000E-04,0.101010E-04,0.102041E-04,0.103093E-04
     &,0.104167E-04,0.105263E-04,0.106383E-04,0.107527E-04,0.108696E-04
     &,0.109890E-04,0.111111E-04,0.112360E-04,0.113636E-04,0.114942E-04
     &,0.116279E-04,0.117647E-04,0.119048E-04,0.120482E-04,0.121951E-04
     &,0.123457E-04,0.125000E-04,0.126582E-04,0.128205E-04,0.129870E-04
     &,0.131579E-04,0.133333E-04,0.135135E-04,0.136986E-04,0.138889E-04
     &,0.140845E-04,0.142857E-04,0.144928E-04,0.147059E-04,0.149254E-04
     &,0.151515E-04,0.153846E-04,0.156250E-04,0.158730E-04,0.161290E-04
     &,0.163934E-04,0.166667E-04,0.169492E-04,0.172414E-04,0.175439E-04
     &,0.178571E-04,0.181818E-04,0.185185E-04,0.188679E-04,0.192308E-04
     &,0.196078E-04,0.200000E-04,0.204082E-04,0.208333E-04,0.212766E-04
     &,0.217391E-04,0.222222E-04,0.227273E-04,0.232558E-04,0.238095E-04
     &,0.243902E-04,0.250000E-04,0.256410E-04,0.263158E-04,0.270270E-04
     &,0.277778E-04,0.285714E-04,0.294118E-04,0.303030E-04,0.312500E-04
     &,0.322581E-04,0.333333E-04,0.344828E-04,0.357143E-04,0.370370E-04
     &,0.384615E-04,0.400000E-04,0.500000E-04,0.600000E-04,0.700000E-04
     &,0.800000E-04,0.100000E-03,0.120000E-03,0.140000E-03,0.170000E-03
     &,0.200000E-03,0.230000E-03,0.260000E-03,0.300000E-03/
      data epdustir/
     & 1.8223799,1.3646801,0.6918110,0.6400360,0.6000190
     &,0.5436750,0.5075750,0.5349990,0.5741300,0.6528150
     &,0.7026240,0.8045070,0.9152510,0.9964510,1.1199900
     &,1.2508200,1.4647000,1.6723000,1.7668300,1.8770100
     &,2.0622301,2.2681601,2.3134601,2.3844199,2.4568000
     &,2.4623699,2.5139699,2.7634001,3.1062000,3.2699001
     &,3.0982399,2.8657801,2.7559099,2.6911600,2.6657801
     &,2.6558800,2.6744900,2.6645300,2.6016099,2.5456300
     &,2.4316900,2.3356500,2.2314799,2.1130800,2.0792401
     &,2.0617199,2.0427401,1.9252200,1.9074700,1.9007500
     &,1.7614900,1.5749300,1.5687200,1.5357200,1.3301001
     &,1.1643699,1.0552500,1.3750900,1.3406500,1.4579200
     &,1.2703400,1.2499200,0.9752970,1.0887500,1.0551100
     &,0.8680390,0.9080590,0.9270510,0.8930280,0.7252980
     &,0.7761930,0.7722600,0.6197060,0.8046810,0.9906260
     &,1.1675600,1.3232200,1.4563700,1.5742900,1.6683700
     &,1.6386300,1.6082799,1.5773100,1.5452300,1.5131000
     &,1.4803500,1.4470100,1.4130900,1.3785900,1.3430200
     &,1.3074800,1.2714100,1.2348300,1.1977800,1.1602700
     &,1.1218801,1.0836500,1.0450701,1.0062100,0.9670900
     &,0.9273350,0.8879470,0.8485210,0.8090940,0.7697440
     &,0.7305380,0.6912350,0.6527100,0.6146080,0.5770320
     &,0.5400900,0.5038810,0.3378440,0.2453480,0.1899010
     &,0.1543380,0.1124400,0.0889283,0.0739620,0.0593130
     &,0.0496814,0.0428006,0.0376627,0.0324836/
      data omegdustir/
     & 0.9455460,0.9255780,0.6110340,0.5893190,0.4896250
     &,0.3736060,0.2719810,0.1933450,0.1591350,0.1584020
     &,0.1842620,0.2511510,0.2616600,0.2910640,0.3243100
     &,0.3532210,0.3744530,0.3675350,0.3614250,0.3678870
     &,0.3820890,0.3866010,0.3842250,0.3821830,0.3814160
     &,0.3779760,0.3837150,0.4103900,0.4473680,0.4952290
     &,0.5245740,0.5125160,0.4980300,0.4833010,0.4788940
     &,0.4769360,0.4898900,0.5542500,0.5697640,0.6304710
     &,0.6634610,0.6500820,0.6469200,0.6292330,0.5847980
     &,0.5847600,0.6166800,0.6622360,0.6334560,0.6793890
     &,0.7708860,0.7351760,0.6715080,0.7118780,0.8014710
     &,0.7203560,0.5123400,0.4070190,0.4575710,0.4804760
     &,0.5493510,0.5477680,0.5564710,0.5089200,0.5221230
     &,0.5291490,0.4292840,0.4024670,0.4110070,0.3907940
     &,0.3085770,0.2905000,0.2411480,0.3065880,0.3493080
     &,0.3836540,0.4065040,0.4211960,0.4325340,0.4396850
     &,0.4379220,0.4360720,0.4341270,0.4320340,0.4298880
     &,0.4276300,0.4252520,0.4227480,0.4201070,0.4172670
     &,0.4143320,0.4112310,0.4079530,0.4044860,0.4008160
     &,0.3968690,0.3927530,0.3883870,0.3837510,0.3788250
     &,0.3735170,0.3679440,0.3620120,0.3556900,0.3489460
     &,0.3417450,0.3339700,0.3257680,0.3170010,0.3076270
     &,0.2976060,0.2868930,0.2218360,0.1676840,0.1255320
     &,0.0942349,0.0544707,0.0330144,0.0214071,0.0119898
     &,0.0073572,0.0047672,0.0033577,0.0022074/
      data gdustir/
     & 0.7727880,0.7891680,0.8149010,0.8143680,0.8141810
     &,0.8135650,0.8105500,0.8033680,0.7936510,0.7786140
     &,0.7633720,0.7376310,0.7232250,0.7045770,0.6794170
     &,0.6532190,0.6216970,0.6061300,0.6006830,0.5882300
     &,0.5672470,0.5528980,0.5515480,0.5486130,0.5436640
     &,0.5431460,0.5333350,0.5041380,0.4678710,0.4283810
     &,0.4064880,0.4156190,0.4282650,0.4454850,0.4535160
     &,0.4693080,0.4752720,0.4885620,0.5029940,0.5025310
     &,0.5189140,0.5377080,0.5533670,0.5710530,0.5778970
     &,0.5762890,0.5706920,0.5760510,0.5783530,0.5683710
     &,0.5678870,0.5920300,0.5962880,0.5909690,0.5981750
     &,0.6187920,0.6339570,0.6020770,0.6024800,0.5877340
     &,0.5985190,0.5957170,0.6125820,0.6001100,0.5978010
     &,0.6065740,0.5983230,0.5910160,0.5887040,0.5944530
     &,0.5810070,0.5743440,0.5758550,0.5615340,0.5437820
     &,0.5243770,0.5039760,0.4829940,0.4609260,0.4389000
     &,0.4355460,0.4321260,0.4286380,0.4247460,0.4211400
     &,0.4174550,0.4136880,0.4098350,0.4058920,0.4015100
     &,0.3974090,0.3932030,0.3888870,0.3844560,0.3799040
     &,0.3748690,0.3700950,0.3651760,0.3601050,0.3548730
     &,0.3490900,0.3435280,0.3377980,0.3318590,0.3256990
     &,0.3192870,0.3122380,0.3053810,0.2982310,0.2907650
     &,0.2829610,0.2747820,0.2274640,0.1885150,0.1562300
     &,0.1306820,0.0917524,0.0638560,0.0486232,0.0311882
     &,0.0219897,0.0157278,0.0129793,0.0101294/
