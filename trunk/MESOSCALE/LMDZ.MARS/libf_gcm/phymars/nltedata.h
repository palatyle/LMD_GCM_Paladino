c INCLUDE nlte_escape.h
       !
       ! Escape functions and VMRs from tabulated values.
       ! Origin: nlte_escape.dat file form Miguel L.
       ! (Y. Wanherdrick, 09/2000)

       integer     np
       parameter ( np = 68 )                ! # data points in tabulation
       real        pnb(np)                  ! Pressure in tabulation
       real        ef1(np)                  ! Esc.funct.#1, tabulated
       real        ef2(np)                  ! Esc.funct.#2, tabulated
       real        co2vmr(np)               ! CO2 VMR tabulated
       real        o3pvmr(np)               ! CO2 VMR tabulated
       real        n2covmr(np)              ! N2+CO VMR tabulated
 
 
      DATA pnb/
     &    12.0000,    11.0000,    10.8000,
     &    10.6000,   10.40000,   10.20000,
     &   10.00000,    9.80000,    9.60000,
     &    9.40000,    9.20000,    9.00000,
     &    8.80000,    8.60000,    8.40000,
     &    8.20000,    8.00000,    7.80000,
     &    7.60000,    7.40000,    7.20000,
     &    7.00000,    6.80000,    6.60000,
     &    6.40000,    6.20000,    6.00000,
     &    5.80000,    5.60000,    5.40000,
     &    5.20000,    5.00000,    4.80000,
     &    4.60000,    4.40000,    4.20000,
     &    4.00000,    3.80000,    3.60000,
     &    3.40000,    3.20000,    3.00000,
     &    2.80000,    2.60000,    2.40000,
     &    2.20000,    2.00000,    1.80000,
     &    1.60000,    1.40000,    1.20000,
     &    1.00000,   0.800000,   0.599999,
     &   0.400000,   0.200000,  0.,
     &  -0.200000,  -0.400001,  -0.600000,
     &  -0.800000,   -1.00000,   -1.20000,
     &   -1.40000,   -1.60000,   -1.80000,
     &   -2.00000,   -3.00000/
 
      DATA ef1/
     &    4.58112E-04,    4.58112E-04,    4.58112E-04,
     &    4.58112E-04,    4.58112E-04,    4.58112E-04,
     &    4.58112E-04,    4.58112E-04,    4.58112E-04,
     &    4.58112E-04,    4.58112E-04,    4.58112E-04,
     &    4.58112E-04,    4.58112E-04,    4.58112E-04,
     &    4.58112E-04,    4.58112E-04,    4.58112E-04,
     &    4.58112E-04,    4.58112E-04,    4.58112E-04,
     &    4.58112E-04,    4.61707E-04,    4.76886E-04,
     &    4.95638E-04,    5.20935E-04,    5.55511E-04,
     &    6.01219E-04,    6.63734E-04,    7.50691E-04,
     &    8.63474E-04,    1.00900E-03,    1.19642E-03,
     &    1.42690E-03,    1.71398E-03,    2.06663E-03,
     &    2.48974E-03,    3.01578E-03,    3.64350E-03,
     &    4.40323E-03,    5.32066E-03,    6.40456E-03,
     &    7.72069E-03,    9.25684E-03,    1.10905E-02,
     &    1.32374E-02,    1.57643E-02,    1.87388E-02,
     &    2.22072E-02,    2.63099E-02,    3.10614E-02,
     &    3.66948E-02,    4.32373E-02,    5.15022E-02,
     &    6.21455E-02,    7.77212E-02,    9.92027E-02,
     &   0.131155,   0.179470,   0.258913,
     &   0.380549,   0.530450,   0.643180,
     &   0.741061,   0.826336,   0.922787,
     &   0.997203,    1.00000/
 
      DATA ef2/
     &    1.98992E-03,    1.98992E-03,    1.98992E-03,
     &    1.98992E-03,    1.98992E-03,    1.98992E-03,
     &    1.98992E-03,    1.98992E-03,    1.98992E-03,
     &    1.98992E-03,    1.98992E-03,    2.01376E-03,
     &    2.09450E-03,    2.22993E-03,    2.42056E-03,
     &    2.68018E-03,    3.04398E-03,    3.43896E-03,
     &    3.80282E-03,    4.20622E-03,    4.76121E-03,
     &    8.01698E-03,    1.19947E-02,    1.69149E-02,
     &    2.24497E-02,    2.85244E-02,    3.54813E-02,
     &    4.39264E-02,    5.46248E-02,    6.75367E-02,
     &    8.29931E-02,    1.01717E-01,   0.123422,
     &   0.148468,   0.177096,   0.208816,
     &   0.244003,   0.282013,   0.322559,
     &   0.365542,   0.410518,   0.457384,
     &   0.505358,   0.553627,   0.600472,
     &   0.644807,   0.687185,   0.727429,
     &   0.764734,   0.798562,   0.828699,
     &   0.854797,   0.877717,   0.897874,
     &   0.915258,   0.929904,   0.942381,
     &   0.952906,   0.962173,   0.970191,
     &   0.976437,   0.981501,   0.985406,
     &   0.988560,   0.991111,   0.993653,
     &   0.995561,    1.00000/
 
      DATA co2vmr/
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.949619,   0.947694,   0.945830,
     &   0.944016,   0.940557,   0.937068,
     &   0.932366,   0.893661/
 
      DATA o3pvmr/
     &    5.06756E-08,    9.16539E-07,    1.68217E-06,
     &    3.00843E-06,    5.03151E-06,    8.07489E-06,
     &    1.23137E-05,    1.79029E-05,    2.45308E-05,
     &    3.27431E-05,    4.26692E-05,    5.44396E-05,
     &    6.78865E-05,    8.33147E-05,    1.00148E-04,
     &    1.18846E-04,    1.39681E-04,    1.64909E-04,
     &    1.93617E-04,    2.25161E-04,    2.60834E-04,
     &    3.01501E-04,    3.44953E-04,    3.91011E-04,
     &    4.40377E-04,    4.90820E-04,    5.43200E-04,
     &    5.95335E-04,    6.45420E-04,    6.93166E-04,
     &    7.43729E-04,    7.93710E-04,    8.44394E-04,
     &    8.94318E-04,    9.44732E-04,    9.94964E-04,
     &    1.04901E-03,    1.10008E-03,    1.16302E-03,
     &    1.22989E-03,    1.30026E-03,    1.37131E-03,
     &    1.45556E-03,    1.55186E-03,    1.66328E-03,
     &    1.77802E-03,    1.91546E-03,    2.07503E-03,
     &    2.24903E-03,    2.47117E-03,    2.71728E-03,
     &    2.99739E-03,    3.33582E-03,    3.73507E-03,
     &    4.20819E-03,    4.76887E-03,    5.42558E-03,
     &    6.20815E-03,    7.14473E-03,    8.28545E-03,
     &    9.51779E-03,    1.08140E-02,    1.22359E-02,
     &    1.36870E-02,    1.51495E-02,    1.67196E-02,
     &    1.85485E-02,    3.36252E-02/
 
      DATA n2covmr/
     &    2.71412E-02,    2.71464E-02,    2.71490E-02,
     &    2.71523E-02,    2.71558E-02,    2.71617E-02,
     &    2.71672E-02,    2.71749E-02,    2.71837E-02,
     &    2.71943E-02,    2.72058E-02,    2.72189E-02,
     &    2.72326E-02,    2.72483E-02,    2.72661E-02,
     &    2.72848E-02,    2.73054E-02,    2.73279E-02,
     &    2.73514E-02,    2.73775E-02,    2.74048E-02,
     &    2.74345E-02,    2.74672E-02,    2.75021E-02,
     &    2.75404E-02,    2.75826E-02,    2.76340E-02,
     &    2.77013E-02,    2.78220E-02,    2.79707E-02,
     &    2.81759E-02,    2.84339E-02,    2.87587E-02,
     &    2.91600E-02,    2.96561E-02,    3.02558E-02,
     &    3.09922E-02,    3.18062E-02,    3.27010E-02,
     &    3.35635E-02,    3.44388E-02,    3.53327E-02,
     &    3.62143E-02,    3.70941E-02,    3.79315E-02,
     &    3.87626E-02,    3.95524E-02,    4.03071E-02,
     &    4.10071E-02,    4.16229E-02,    4.21231E-02,
     &    4.25167E-02,    4.27964E-02,    4.29773E-02,
     &    4.30488E-02,    4.29638E-02,    4.28049E-02,
     &    4.26788E-02,    4.26822E-02,    4.29426E-02,
     &    4.34634E-02,    4.42559E-02,    4.53038E-02,
     &    4.65879E-02,    4.80262E-02,    4.96303E-02,
     &    5.14885E-02,    6.91651E-02/
 
