// profile for WRF initialization
// update: january 2008


//// ----------------------------------------------
//// Allow user to define input parameters
//// ----------------------------------------------
//lat=evstr(x_dialog('Latitude (deg)?','20'));
//lon=evstr(x_dialog('East Longitude (deg)?','0'));
//l_s=evstr(x_dialog('Aerocentric longitude (deg)?','90'));
//lt=evstr(x_dialog('Local time ?','12'));
//
//alt=evstr(x_dialog('Altitude top (km)?','50'));
//resol=evstr(x_dialog('Vertical resolution (points)?','100'));


[finput]=file('open', 'input_coord', 'old');
[vinput]=read(finput,6,1);

lon   = vinput(1);
lat   = vinput(2);
l_s   = vinput(3);
lt    = vinput(4);
alt   = vinput(5);
resol = vinput(6);

disp(lon)
disp(lat)
disp(l_s)
disp(lt)
disp(alt)
disp(resol)

//vert_coord=evstr(x_dialog('Vertical coordinate type (1-4)?','3'));
vert_coord = 3;
//dust_scenario=evstr(x_dialog('Dust scenario (1-8)?','3'));
dust_scenario = 3;
//hireskey=evstr(x_dialog('High resolution (0-1)?','1'));
hireskey=1;


// create atltitude array at the required resolution
alt_profile = [0,1,5,10,20,50,100,[linspace(200,alt*1000,resol)]];
resol=resol+7;



// MCD profiles
//// - profils_mcd(Ls,Loct,Lat,Lon,dust,zkey,hireskey,xz)
//getf('./profils_mcd.sci');
exec('./profils_mcd.sci');//,7);
[p,m] = profils_mcd(l_s,lt,lat,lon,dust_scenario,vert_coord,hireskey,alt_profile);


// In Scilab arrays have to be initialized 
utab=ones(max(size(alt_profile)));
vtab=ones(max(size(alt_profile)));
ttab=ones(max(size(alt_profile)));


// % "-espace occupe- . -precision apres la virgule- -type de nombre-"
// \n pour passer a la ligne

fid = mtlb_fopen("input_sounding","w")
fid2= mtlb_fopen("input_more","w")
fid3= mtlb_fopen("input_therm","w")
  // pression de surface (1ere couche)
     ref_p = p(1,1,1)  // equiv. a m(1,1,19)
     ref_p_mbar = ref_p/100;
     ref_dummy = 0;
  // temp_pot de surface (1ere couche)
     disp(m(1,1,15)); // equiv. a p(1,1,3)
     //ref_temppot = m(1,1,15)*(610./ref_p)^(1.0/4.4);
     ref_temppot = m(1,1,15)*(610./ref_p)^(1.0/3.9);
  // ecriture  	
     mfprintf(fid,"%10.2f",ref_p_mbar);
     mfprintf(fid,"%12.2f",ref_temppot);
     mfprintf(fid,"%12.2f\n",ref_dummy);

  // ecriture du profil
for i=1:resol
  dummy=0;
  pression=p(1,i,1);
  altitude=m(1,i,2);    // altitude zero MOLA
  temperature=p(1,i,3);
  density=p(1,i,2);
  cp = m(1,i,8);
  rr = m(1,i,49);
  rcp = rr / cp ;
  //temppot=temperature*(610./pression)^rcp;
  //temppot=temperature*(610./pression)^(1.0/4.4);
  temppot=temperature*(610./pression)^(1.0/3.9);
  u=p(1,i,4);
  v=p(1,i,5);

      mfprintf(fid,"%10.2f",altitude);
      mfprintf(fid,"%12.2f",temppot);
      ttab(i)=temppot;
      mfprintf(fid,"%12.2f",dummy);
      mfprintf(fid,"%12.2f",u);
      utab(i)=u;
      mfprintf(fid,"%12.2f\n",v);
      vtab(i)=v;
      mfprintf(fid3,"%12.2f",rr);
      mfprintf(fid3,"%12.2f",cp);
      mfprintf(fid3,"%18.6e",pression);
      mfprintf(fid3,"%18.6e",density);
      mfprintf(fid3,"%12.2f\n",temperature);
end
altimetry=m(1,1,2)
//altimetry=0. //OK
surftemp=m(1,1,15)
mfprintf(fid2,"%10.2f",altimetry);
mfprintf(fid2,"%10.2f",surftemp);

exit

// Plot 
alt_profile=alt_profile/1000;
xset("window",1)
style = [1,5];
strf = "111";
leg = "Zonal@Meridional";
rect = [0,min([utab,vtab])-5,max(alt_profile),max([utab,vtab])+5]; //BEWARE : xmin ymin xmax ymax
plot2d([alt_profile',alt_profile'],[utab,vtab],style,strf,leg,rect);
	// BEWARE : column vector in plot2d
xtitle("MCD 4.2 wind profiles at longitude "+string(lon)+" deg, latitude "+string(max(lat))+" deg, localtime "+string(lt)+" h, ls "+string(l_s)+" deg","Altitude (km)","Wind (m/s)");

xset("window",2)
style = [1];
strf = "111";
leg = "";
rect = [0,min(ttab)-5,max(alt_profile),max(ttab)+5]; //BEWARE : xmin ymin xmax ymax
plot2d(alt_profile',ttab,style,strf,leg,rect);
	// BEWARE : column vector in plot2d
xtitle("MCD 4.2 temperature at longitude "+string(lon)+" deg, latitude "+string(max(lat))+" deg, localtime "+string(lt)+" h, ls "+string(l_s)+" deg","Altitude (km)","Temperature (K)");

//exit



