      program mammary
c-------------------------------------------------------------------------c
c     program to simulate the dynamics of epithelial turnover in the mouse
c     mammary gland epithelium
c
c     at each periodic cycle the algorithm involves choosing random, but 
c     separated domains, which then enter into a phase of growth and
c     regression
c
c     for details of the dynamics, see main text, supplementary theory and
c     comments below. In short, the algorithm replaces a region of size 
c     ndomain/2 by duplicating a domain of ndomain/4 cells on either side. 
c     Then further rounds of local stochastic loss and replacement remove
c     pathological commensurability effects
c
c     the resulting clone size distributions and associated averages are
c     collated and output
c-------------------------------------------------------------------------c
      implicit double precision (a-h,o-z)
c-------------------------------------------------------------------------c
c     nsamp = number of samples for simulation
c     dq    = probability that a given cell enters cycle per estrous cycle
c-------------------------------------------------------------------------c
      parameter (nsamp=1000,dcysc=1.0d0,dq=0.1d0*dcysc)
c-------------------------------------------------------------------------c
c     nx      = total linear dimension of one-dimensional lattice
c     ndomain = total size of the domain that enters into cycle
c     ntots   = total number of domains that become active during one cycle
c     nmx     = maximum clone size for arrays
c-------------------------------------------------------------------------c
      parameter (nx=1000000,ndomain=1000,itime=6,
     +     ntots=int(dq*dble(nx)/dble(ndomain)),nmx=10000)
c-------------------------------------------------------------------------c
c     latt  = array of x-coordinates indexing centre of activated domains
c     nlat  = barcode of clone at position nx
c     nclsz = array for clone size distribution
c-------------------------------------------------------------------------c
      integer latt(ntots),icycle(0:itime),nlat(nx),
     +    nclsz(nx,itime)
c-------------------------------------------------------------------------c
c     cldist = clone size distribution
c     dcum   = cumulative clone size dsitribution
c     dav    = store for averages
c     davln  = store for log averages
c-------------------------------------------------------------------------c
      real*8 cldist(0:nmx,itime),dcum(itime),dav(itime),
     +   davln(itime),davlns(itime)
c---- seed random number generator
      call random_seed
c---- cycle count for outputs (in units of estrous cycle number)
      icycle(0)=0
      icycle(1)=4
      icycle(2)=8
      icycle(3)=16
      icycle(4)=32
      icycle(5)=64
      icycle(6)=128
c---- set counter for maximum clone size
      max1=0
      max2=0
c---- zero clone size distribution
      do 30 l1=0,nmx
         do 20 l2=1,itime
            cldist(l1,l2)=0.0d0
 20      continue
 30   continue
c-------------------------------------------------------------------------c
c     loop through samples
c-------------------------------------------------------------------------c
      do 650 ns=1,nsamp
         if (dble(ns/10).eq.dble(ns)/dble(10))
     +      write(*,*) 'sample number=',ns
c---- initialise unique barcode for stem cell nlat()
         do 50 lx=1,nx
            nlat(lx)=lx
 50      continue
c---- loop over cycle time points
         do 500 it=1,itime
c---- zero clone size distribution for individual barcodes
            do 55 l1=1,nx
               nclsz(l1,it)=0 
 55         continue
c---- loop over estrous cycles
            do 400 ncy=icycle(it-1),icycle(it)
c-------------------------------------------------------------------------c
c     get set of lattice sites latt(.) for activation
c-------------------------------------------------------------------------c
               npt=0
 60            call random_number(random)
               ix=int(random*dble(nx))+1
               if (npt.eq.0) goto 80
c---- check that chosen sites are separated by a distance >ndomain
               iflag=0
               do 70 l1=1,npt
                  if (min(iabs(latt(l1)-ix),nx-iabs(latt(l1)-ix)).lt.
     +               ndomain) iflag=1
 70            continue
               if (iflag.eq.1) goto 60
 80            npt=npt+1
               latt(npt)=ix
               if (npt.lt.ntots) goto 60
c-------------------------------------------------------------------------c
c     effect stem cell loss/replacement by duplication from wings of domain
c        i.e., ndomain/4 cells on either side of domain expand by duplication 
c        to replace the ndomain/2 cells in the middle of the domain
c-------------------------------------------------------------------------c
               do 100 l1=1,ntots ! loop over all domains
                  do 88 l2=1,ndomain/4
                     nxi1=mod(latt(l1)+ndomain/4-1+l2-1+nx-1,nx)+1
                     nxf1=mod(latt(l1)+2*l2-2+nx-1,nx)+1
                     nxf2=mod(latt(l1)+2*l2-1+nx-1,nx)+1
                     nlat(nxf1)=nlat(nxi1) ! copy one stem cell
                     nlat(nxf2)=nlat(nxi1) ! and duplicate other
                     nxi1=mod(latt(l1)-ndomain/4-l2+nx-1,nx)+1
                     nxf1=mod(latt(l1)-2*l2+nx-1,nx)+1
                     nxf2=mod(latt(l1)-2*l2+nx-1-1,nx)+1
                     nlat(nxf1)=nlat(nxi1) ! copy one stem cell
                     nlat(nxf2)=nlat(nxi1) ! and duplicate other
 88               continue
c-------------------------------------------------------------------------c     
c     then effect further rounds of loss/replacement inside domain to avoid
c     spurious commensurability effects
c-------------------------------------------------------------------------c
                  do 91 l3=1,ndomain/4 ! adjustable
                     call random_number(random)
                     ixi=int(random*dble(ndomain))+1
                     nxi=mod(latt(l1)-ndomain/2+ixi+nx-1,nx)+1
                     call random_number(random)
                     if (random.gt.0.5d0) then ! replacement from left/right
                        nxf=mod(nxi+1+nx-1,nx)+1
                     else
                        nxf=mod(nxi-1+nx-1,nx)+1
                     endif
                     nlat(nxf)=nlat(nxi)
 91               continue
c---- output sample domain for diagnostics -------------------------------c
c                  do 99 l2=1,2*ndomain
c                     nxi=mod(latt(l1)-ndomain+l2+nx-1,nx)+1
c                     write(16,*) l2,nlat(nxi)
c 99               continue
c                  stop
 100           continue
 400        continue
c-------------------------------------------------------------------------c
c     update clone size arrays
c-------------------------------------------------------------------------c
            do 420 l1=1,nx
               nclsz(nlat(l1),it)=nclsz(nlat(l1),it)+1
 420        continue
c-------------------------------------------------------------------------c
c     determine sample averaged clone size distribution
c-------------------------------------------------------------------------c
            do 430 l1=1,nx
               nstem=nclsz(l1,it)
               if (nstem.gt.max1) max1=nstem
               cldist(nstem,it)=cldist(nstem,it)+
     +            1.0d0/dble(nx)/dble(nsamp)
 430        continue
c---- determine total surviving clone number
            nclone1=0 ! stem cell
            do 450 l1=1,nx
               if (nclsz(l1,it).gt.0) nclone1=nclone1+1
 450        continue
c---- checksum on normalisation
c            dnorm1=0.0d0
c            do 470 l1=0,nmx
c               dnorm1=dnorm1+cldist(l1,it)*dble(nsamp)/dble(ns)
c 470        continue
c            write(*,*) 'cycle no.=',it,'stem clone no.=',nclone1,
c     +         nx,'normalisations=',dnorm1
 500     continue
 650  continue
c-------------------------------------------------------------------------c
c     output distributions and averages
c-------------------------------------------------------------------------c
      open(unit=21,file='cldist.dat',status='replace')
      do 680 it=1,itime
         dav(it)=0.0d0 ! average clone size
         davln(it)=0.0d0 ! average log total clone size
         davlns(it)=0.0d0 ! average square log total clone size
         dcum(it)=1.0d0 ! cumulative clone size distribution
 680  continue
      do 700 l1=1,max1
         write(21,*) l1,(cldist(l1,it)/(1.0d0-cldist(0,it)),
     +      it=1,itime), ! surviving stem cell clone size distribution
     1       (dmax1(dcum(it)-cldist(l1,it)/(1.0d0-cldist(0,it)),0.0d0),
     2      it=1,itime) ! cumulative stem cell clone size distribution
         do 690 it=1,itime
            dcum(it)=dcum(it)-
     +         cldist(l1,it)/(1.0d0-cldist(0,it))
            dav(it)=dav(it)+dble(l1)*
     +         cldist(l1,it)/(1.0d0-cldist(0,it))
            davln(it)=davln(it)+dlog(dble(l1))*
     +         cldist(l1,it)/(1.0d0-cldist(0,it))
            davlns(it)=davlns(it)+dlog(dble(l1))**2.0d0*
     +         cldist(l1,it)/(1.0d0-cldist(0,it))
 690     continue
 700  continue
      close(unit=21)
c-------------------------------------------------------------------------c
c     output averages for stem and total:
c        cycle number
c        surviving clone fraction
c        average clone size
c        average log clone size
c        variance of log clone size
c-------------------------------------------------------------------------c
      open(unit=22,file='clsize.dat',status='replace')
      do 750 it=1,itime
         write(22,*) icycle(it),
     +      1.0d0-cldist(0,it),dav(it),davln(it),
     1      davlns(it)-davln(it)**2.0d0
 750  continue
      close(unit=22)
c-------------------------------------------------------------------------c
c     output sample clones
c-------------------------------------------------------------------------c
      open(unit=16,file='sample.dat',status='replace')
      do 800 l1=1,10*ndomain
         write (16,*) l1,nlat(l1)
 800  continue
      close(unit=16)
      stop
      end
            
          
