
      integer npres                ! Number of pressures in NIR correction
      parameter (npres=42)         ! table

      common /NIRdata/ pres1d,corgcm,oco21d,alfa,p1999
      real    pres1d(npres)
      real    corgcm(npres)
      real    oco21d(npres),alfa(npres),p1999(npres)
