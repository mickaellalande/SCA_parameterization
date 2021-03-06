!
! $Id: ini_histday_seri.h 2399 2015-11-20 16:23:28Z emillour $
!
!ym Ne fonctionnera pas en mode parallele
      IF (is_sequential) THEN
      
      IF (type_run.EQ."AMIP") THEN
!
       zstophy = pdtphys
       zout = ecrit_day
!
         idayref = day_ref
         CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
!
         CALL gr_fi_ecrit(1,klon,nbp_lon,nbp_lat,longitude_deg,zx_lon)
         DO i = 1, nbp_lon
            zx_lon(i,1) = longitude_deg(i+1)
            zx_lon(i,nbp_lat) = longitude_deg(i+1)
         ENDDO
         DO ll=1,klev
            znivsig(ll)=REAL(ll)
         ENDDO
         CALL gr_fi_ecrit(1,klon,nbp_lon,nbp_lat,latitude_deg,zx_lat)
!
         imin_debut=1 
         nbpti=1
         jmin_debut=1 
         nbptj=1
!
         CALL histbeg("histday_seri.nc",  &
                       nbp_lon,zx_lon(:,1),nbp_lat,zx_lat(1,:), &
                       imin_debut,nbpti,jmin_debut,nbptj, &
                       itau_phy, zjulian, dtime, &
                       nhori, nid_day_seri)
!
         CALL histvert(nid_day_seri, "presnivs",  &
                      "Vertical levels","mb", &
                       klev, presnivs/100., nvert)
!
         CALL histdef(nid_day_seri, "bilTOA",  &
                      "Net radiation at model top", "W/m2", &
                      nbp_lon,nbp_lat,nhori, 1,1,1, -99, 32,  &
                      "ave(X)", zstophy,zout)
!
         CALL histdef(nid_day_seri, "bils",  &
                      "Net downward energy flux at surface","W/m2", &
                      nbp_lon,nbp_lat,nhori, 1,1,1, -99, 32,  &
                      "ave(X)", zstophy,zout)
!
         CALL histdef(nid_day_seri, "ecin",  &
                      "Total kinetic energy (per unit area)","J/m2", &
                      nbp_lon,nbp_lat,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
!
!IM 151004 BEG
         IF(1.EQ.0) THEN
!
         CALL histdef(nid_day_seri, "momang",  &
                     "Total relative angular momentum (per unit area)", &
                     "kg/s", &
                      nbp_lon,nbp_lat,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
!
         CALL histdef(nid_day_seri, "frictor",  &
                     "Friction torque (per unit area)", "N/m", &
                      nbp_lon,nbp_lat,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
!
         CALL histdef(nid_day_seri, "mountor",  &
                     "Mountain torque (per unit area)", "N/m", &
                      nbp_lon,nbp_lat,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
!
         ENDIF !(1.EQ.0) THEN
!
         CALL histdef(nid_day_seri, "momang",  &
                     "Axial angular momentum (per unit area)", &
                     "kg/s", &
                      nbp_lon,nbp_lat,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
!
         CALL histdef(nid_day_seri, "torsfc",  &
              "Total surface torque (including mountain torque)", "N/m", &
                      nbp_lon,nbp_lat,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
!
!IM 151004 END        
!
         CALL histdef(nid_day_seri, "tamv",  &
                      "Temperature (mass-weighted vert. ave)", "K", &
                      nbp_lon,nbp_lat,nhori, 1,1,1, -99, 32, &
                      "ave(X)", zstophy,zout)
!
         CALL histdef(nid_day_seri, "psol",  &
                      "Surface pressure", "Pa", &
                      nbp_lon,nbp_lat,nhori, 1,1,1, -99, 32,  &
                      "ave(X)", zstophy,zout)
!
         CALL histdef(nid_day_seri, "evap",  &
                      "Evaporation and sublimation (per unit area)",  &
                      "kg/(m2*s)", &
                      nbp_lon,nbp_lat,nhori, 1,1,1, -99, 32,  &
                      "ave(X)", zstophy,zout)
!
!          call histdef(nid_day_seri, 
!    .         "SnowFrac", 
!    .         "Snow-covered area ", "%",  
!    .         nbp_lon,nbp_lat,nhori, 1,1,1, -99, 32,
!    .         "ave(X)", zstophy,zout)
!
!        CALL histdef(nid_day_seri, "snow_depth", 
!IM 080904  .                "Snow Depth (water equivalent)", "m",
!IM 191104  .                "Snow Depth (water equivalent)", "kg/m2",
!    .                "Snow Mass", "kg/m2",
!    .                nbp_lon,nbp_lat,nhori, 1,1,1, -99, 32, 
!    .               "ave(X)", zstophy,zout)
!
           call histdef(nid_day_seri,  &
               "tsol_"//clnsurf(is_oce),  &
               "SST over open (ice-free) ocean ", "K",   &
               nbp_lon,nbp_lat,nhori, 1,1,1, -99, 32, &
               "ave(X)", zstophy,zout)
!
!=================================================================
!
         CALL histend(nid_day_seri)
!
!=================================================================
      ENDIF ! fin de test sur type_run.EQ.AMIP
      
      ENDIF ! is_sequential
