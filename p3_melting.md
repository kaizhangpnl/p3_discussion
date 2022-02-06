## Melting 

The links below point to the NGD Phase II baseline version (with p3). 

## ice_melting before sedimentation  

https://github.com/E3SM-Project/E3SM/blob/7792b21d910705e1924c0014c6913eaefc721005/components/eam/src/physics/cam/micro_p3.F90#L781-L786

```
      !............................................................
      ! melting
      call ice_melting(rho(k),t_atm(k),pres(k),rhofaci(k),&
           table_val_qi2qr_melting,table_val_qi2qr_vent_melt,latent_heat_vapor(k),latent_heat_fusion(k),dv,sc,mu,kap,&
           qv(k),qi_incld(k),ni_incld(k),&
           qi2qr_melt_tend,ni2nr_melt_tend)
```

https://github.com/E3SM-Project/E3SM/blob/7792b21d910705e1924c0014c6913eaefc721005/components/eam/src/physics/cam/micro_p3.F90#L2354-L2399

```
subroutine ice_melting(rho,t_atm,pres,rhofaci,    &
table_val_qi2qr_melting,table_val_qi2qr_vent_melt,latent_heat_vapor,latent_heat_fusion,dv,sc,mu,kap,qv,qi_incld,ni_incld,    &
           qi2qr_melt_tend,ni2nr_melt_tend)
   ! melting
   ! need to add back accelerated melting due to collection of ice mass by rain (pracsw1)
   ! note 'f1pr' values are normalized, so we need to multiply by N
   ! currently enhanced melting from collision is neglected
   ! include RH dependence

   implicit none

   real(rtype), intent(in) :: rho
   real(rtype), intent(in) :: t_atm
   real(rtype), intent(in) :: pres
   real(rtype), intent(in) :: rhofaci
   real(rtype), intent(in) :: table_val_qi2qr_melting ! melting
   real(rtype), intent(in) :: table_val_qi2qr_vent_melt ! melting (ventilation term)
   real(rtype), intent(in) :: latent_heat_vapor
   real(rtype), intent(in) :: latent_heat_fusion
   real(rtype), intent(in) :: dv
   real(rtype), intent(in) :: sc
   real(rtype), intent(in) :: mu
   real(rtype), intent(in) :: kap
   real(rtype), intent(in) :: qv
   real(rtype), intent(in) :: qi_incld
   real(rtype), intent(in) :: ni_incld

   real(rtype), intent(out) :: qi2qr_melt_tend
   real(rtype), intent(out) :: ni2nr_melt_tend

   real(rtype) :: qsat0

   if (qi_incld .ge.qsmall .and. t_atm.gt.T_zerodegc) then
      qsat0 = qv_sat( T_zerodegc,pres,0 )

      qi2qr_melt_tend = ((table_val_qi2qr_melting+table_val_qi2qr_vent_melt*bfb_cbrt(sc)*bfb_sqrt(rhofaci*rho/mu))*((t_atm-   &
      T_zerodegc)*kap-rho*latent_heat_vapor*dv*(qsat0-qv))*2._rtype*pi/latent_heat_fusion)*ni_incld

      qi2qr_melt_tend = max(qi2qr_melt_tend,0._rtype)
      ni2nr_melt_tend = qi2qr_melt_tend*(ni_incld/qi_incld)

   endif

   return

end subroutine ice_melting
```


## ice_complete_melting after sedimentation  


https://github.com/E3SM-Project/E3SM/blob/7792b21d910705e1924c0014c6913eaefc721005/components/eam/src/physics/cam/micro_p3.F90#L1591-L1594


```
       !.........................................................
       ! Instantenous melting of ice/snow at T = t_snow_melt = 2c    
       call ice_complete_melting(kts,kte,ktop,kbot,kdir,qi(i,:),ni(i,:),qm(i,:),latent_heat_fusion(i,:),exner(i,:),th_atm(i,:), & 
            qr(i,:),nr(i,:),qc(i,:),nc(i,:))
```


https://github.com/E3SM-Project/E3SM/blob/7792b21d910705e1924c0014c6913eaefc721005/components/eam/src/physics/cam/micro_p3.F90#L4452-L4495


```
subroutine ice_complete_melting(kts,kte,ktop,kbot,kdir,qi,ni,qm,latent_heat_fusion,exner,th_atm, & 
            qr,nr,qc,nc)

   implicit none
   
   integer, intent(in) :: kts, kte
   integer, intent(in) :: ktop, kbot, kdir
   real(rtype), intent(in), dimension(kts:kte) :: exner,latent_heat_fusion
   real(rtype), intent(inout), dimension(kts:kte) :: qi, ni, qc, nc, qr, nr, qm, th_atm
   
   real(rtype) :: t_snow_melt,del_mass,del_num,rv_tmp,rv 
   integer :: k

   t_snow_melt = 273.15 + 2.0_rtype         
   k_loop_mlt:  do k = kbot,ktop,kdir
      if (qi(k).ge.qsmall .and. (th_atm(k)/exner(k)) > t_snow_melt) then
         del_mass = qi(k)
         del_num = ni(k)
         rv_tmp = 3.0_rtype*qi(k)/ni(k)/4.0_rtype/pi/900.0_rtype  ! density of pure ice [kg/m3]
         rv = 1.0e6_rtype*bfb_cbrt(rv_tmp)                        ! in [um]
         if((qm(k)/qi(k)) < 0.1_rtype)then
            if(rv < 100.0_rtype)then
               ! ... Ice crystas
               qc(k) = qc(k) + del_mass
               nc(k) = nc(k) + del_num
            else
               ! ... Unrimed snow
               qr(k) = qr(k) + del_mass
               nr(k) = nr(k) + del_num
            endif
         else
            ! ... Medium rimed snow
            qr(k) = qr(k) + del_mass
            nr(k) = nr(k) + del_num
         endif
         qi(k) = 0.0_rtype
         ni(k) = 0.0_rtype
         qm(k) = 0.0_rtype
         th_atm(k) = th_atm(k) - del_mass*latent_heat_fusion(k)/cp
      endif
   enddo k_loop_mlt
   
return
end subroutine ice_complete_melting  
```



