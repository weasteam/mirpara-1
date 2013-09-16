c     Program ct2out - Creates an ASCII display of a nucleic acid
c                      secondary structure from a ct file.
c
c     Input - a ct file containing one or more structures
c             (only standard input is accepted)
c     Output - text output (standard output)
c     Objective - This program creates a simple structure
c                 plot in text format.
c
      implicit integer (a-z)
      parameter (maxn=20000)
      character*1 lorc
      character*1 seq(2*maxn)
      character*80 ctrec,ctlabel,oldlabel
      real alog10
      integer basepr(maxn),hstnum(2*maxn)
      logical flag
      data lorc/'l'/

c
c     Read foldings one by one.
c
      numseq = 0
      nold = -1
      do i = 1,80
         oldlabel(i:i) = ' '
         ctlabel(i:i) = ' '
      enddo
 10   read(5,1035,end=999) ctrec
      read(ctrec,*) n
      k=int(alog10(float(n-1) + 0.1)) + 1
      i=1
      do while (ctrec(i:i).eq.char(9).or.ctrec(i:i).eq.' ')
         i = i + 1
      enddo
      do j = i+k,80
         ctlabel(j-i-k+1:j-i-k+1) = ctrec(j:j)
      enddo
      lablen = 80
      do while (ctlabel(lablen:lablen).eq.' ') 
         lablen = lablen - 1
      enddo
      if (n.ne.nold.or.ctlabel.ne.oldlabel) then
         numseq = numseq + 1
         iter = 0
         oldlabel = ctlabel
         nold = n
      endif
      do i = 1,n
         basepr(i) = 0
      enddo
      do i = 1,n
         read(5,1035,end=998,err=997) ctrec
 1035    format(a80)
         read(ctrec,*,err=997) k,seq(i)
         kb = index(ctrec,seq(i))
         read(ctrec(kb+1:80),*,err=997) iup,idown,basepr(i),hstnum(i)
         if (i.ne.k) go to 997
         if (basepr(i).gt.0) basepr(basepr(i)) = i
c     Check for circular sequence
         if (i.eq.1.and.iup.eq.n) lorc = 'c'
         if (i.eq.n) then
            if(lorc.eq.'c'.and.idown.ne.1.or.lorc.eq.'l'.and.idown.eq.1)
     .           then
               write(6,*) 'Sequence is declared both linear and ',
     .              'circular. Error.'
               call exit(1)
            endif
         endif
      enddo
      iter = iter + 1
c
c     Check for knots
c
      do i = 1,n
         if (basepr(i).ne.0) then
            j = max(i,basepr(i))
            ip = min(i,basepr(i))
            if (basepr(ip).ne.j.or.basepr(j).ne.ip) then
               write(6,1041) hstnum(ip),hstnum(j)
               call exit(1)
            endif
            do k = ip+1,j-1
               if (basepr(k).ne.0) then
                  l = basepr(k)
                  if (l.le.ip.or.l.ge.j) then
                     write(6,1042) hstnum(ip),hstnum(j),hstnum(k),
     .                    hstnum(l)
                     call exit(1)
                  endif
               endif
            enddo
         endif
      enddo
 1041 format(' Base pair ',i5,'.',i5,' has 5'' or 3'' end in ',
     .     'another base pair.')
 1042 format(' Base pair ',i5,'.',i5,' conflicts with ',i5,'.',i5)
      
c
c     Double up if circular
c
      if (lorc.eq.'c') then
         do i = 1,n
            basepr(i+n) = 0
            seq(i+n) = seq(i)
            hstnum(i+n) = hstnum(i)
         enddo
c
c     Find a base pair that closes a hairpin loop.
c     First find a base pair.
c
         i = 1
         do while (basepr(i).eq.0.and.i.le.n)
            i = i + 1
         enddo
         if (i.gt.n) then
c
c     No structure !
c
            write(6,1050) numseq,iter,ctlabel
 1050       format('Sequence ',i3,' Structure ',i5,':',/,a80,/,
     .           'No structure found.')
            go to 10
         endif
c
         flag = .false.
         do while (.not.flag)
            j = i + 1
            do while (basepr(j).eq.0.and.j.lt.basepr(i))
               j = j + 1
            enddo
            if(j.eq.basepr(i)) then
               flag = .true.
            else
               i = j
            endif
         enddo
c
c     i.basepr(i) close a hairpin loop
c
c     Circular permutation of structure : start at basepr(i)
c
         do j = 1,i
            if (basepr(j).gt.0) then
               if (basepr(j).lt.i) then
                  basepr(n+j) = basepr(j) + n
                  basepr(n+basepr(j)) = n + j
               else
                  basepr(basepr(j)) = n + j
                  basepr(n+j) = basepr(j)
               endif
            endif
         enddo
c
c     Restore folding in circular case
c
         do j = n + 1,n + i
            if(basepr(j).le.n) basepr(basepr(j)) = j - n
         enddo
         
      endif

c
      write(6,2010) numseq, iter
 2010 format('Sequence ',i3,' Structure ',i6)

c
c     Text output.
c
      call txtout(1,n,seq,hstnum,basepr,ctlabel,lablen,err)
c
c     Read next structure
c
      write(6,2015)
 2015 format('__________________________________________________',/)
      go to 10
c
c     Error in CT file.
c
 997  write(6,9080)
 9080 format('STOP: Error in CT file.')
      call exit(1)
c     
c     Error - incomplete CT file.
c     
 998  write(6,9090)
 9090 format('STOP:  Premature end of CT file.')
      call exit(1)
c     
 999  stop
      end
 
c     Initialize the stack.
      subroutine initst
      implicit integer (a-z)
      dimension stk(500,4)
      common /stk/ stk,sp
 
      sp = 0
      return
      end
c     Add A,B,C,D to the bottom of the stack.
      subroutine push(a,b,c,d)
      implicit integer (a-z)
      dimension stk(500,4)
      common /stk/ stk,sp
      
      sp = sp + 1
      if (sp.gt.500) then
         write (6,*) 'STOP: STACK OVERFLOW'
         call exit(1)
      endif
      stk(sp,1) = a
      stk(sp,2) = b
      stk(sp,3) = c
      stk(sp,4) = d
      return
      end
c     Retrieve A,B,C,D from the bottom of the stack and decrease the
c     stack size by one.
      function pull(a,b,c,d)
      implicit integer (a-z)
      dimension stk(500,4)
      common /stk/ stk,sp
 
      if (sp.eq.0) then
         pull = 1
         return
      endif
      a = stk(sp,1)
      b = stk(sp,2)
      c = stk(sp,3)
      d = stk(sp,4)
      sp = sp - 1
      pull = 0
      return
      end

c     Text output of a secondary structure.
      subroutine txtout(n1,n2,seq,hstnum,basepr,ctlabel,lablen,error)
c
      implicit integer (a-z)
      character*1 array(6,5000),seq(1)
      integer hstnum(1),basepr(1)
      character*80 ctlabel
c
      data amax/5000/
c
c     Write sequence label and computed energy 
c

      hstn1 = hstnum(n1)
      hstn2 = hstnum(n2)
      write(6,102) (ctlabel(i:i),i=1,lablen)
 102  format(80a1)
      write(6,103) hstn1,hstn2
 103  format('Bases',i5,' to',i5,/)

c     
c     Initialize traceback
c     
      call initst
      call push(n1,n2,0,0)
c     
c     Fill in output array
c
      do while (pull(i,j,countr,xx).eq.0)
c
c     Look for dangling ends
c
         ip = i
         jp = j
c        Chew off dangling 5' bases
         do while (basepr(ip).eq.0.and.ip.lt.j)
            ip = ip+1
         enddo
         if (ip.eq.jp) then
c        Hairpin loop found - dump array
            size = j - i + 1
            hsize = (size-1)/2
            if (hsize.gt.0) then
               do k = 1,hsize
                  array(1,countr+k) = ' '
                  array(2,countr+k) = seq(hstnum(i+k-1))
                  if(10*(hstnum(i+k-1)/10).eq.hstnum(i+k-1))
     .                 call digit(1,countr+k,hstnum(i+k-1),amax,array)
                  array(3,countr+k) = ' '
                  array(4,countr+k) = ' '
                  array(5,countr+k) = seq(hstnum(j-k+1))
                  array(6,countr+k) = ' '
                  if (10*(hstnum(j-k+1)/10).eq.hstnum(j-k+1))
     .                 call digit(6,countr+k,hstnum(j-k+1),amax,array)
               enddo
               countr = countr + hsize + 1
               if (2*hsize.eq.size-1) then
                  array(1,countr) = ' '
                  array(2,countr) = ' '
                  array(3,countr) = '\\'
                  array(4,countr) = seq(hstnum(i+hsize))
                  array(5,countr) = ' '
                  array(6,countr) = ' '
               else
                  array(1,countr) = ' '
                  array(2,countr) = ' '
                  array(3,countr) = seq(hstnum(i+hsize))
                  array(5,countr) = ' '
                  array(6,countr) = ' '
                  if (10*(hstnum(i+hsize)/10).eq.hstnum(i+hsize))
     .                 call digit(1,countr,hstnum(i+hsize),amax,array)
                  array(4,countr) = seq(hstnum(i+hsize+1))
                  if (10*(hstnum(i+hsize+1)/10).eq.hstnum(i+hsize+1))
     .                 call digit(6,countr,hstnum(i+hsize+1),amax,array)
               endif
            endif
            call dump_array(amax,array,countr)
            do k = 1,countr
               do k1 = 1,6
                  array(k1,k) = ' '
               enddo
            enddo
         else
c           Chew off dangling 3' bases
            do while (basepr(jp).eq.0)
               jp = jp-1
            enddo
c
c           Test for bifurcation
c
            if (basepr(ip).lt.basepr(jp)) then
               do k1 = 1,6
                  do k2 = 1,2
                     array(k1,countr+k2) = ' '
                  enddo
               enddo
               if (basepr(i).gt.0) then
                  array(3,countr+1) = '-'
                  array(3,countr+2) = '-'
               else
                  array(2,countr+1) = '.'
                  array(2,countr+2) = '-'
               endif
               array(5,countr+1) = '\\'
               countr = countr + 2
               call push(basepr(ip)+1,j,countr,0) 
               call push(i,basepr(ip),countr,0)
            else
               k = max0(ip-i,j-jp)
               if (k.gt.0) then
c              Unpaired bases are in an interior or bulge loop.
                  do kk = 1,k
                     array(1,countr+kk) = ' '
                     array(3,countr+kk) = ' '
                     array(4,countr+kk) = ' '
                     array(6,countr+kk) = ' '
                     if (i+kk-1.lt.ip) then
                        array(2,countr+kk) = seq(hstnum(i+kk-1))
                        if (10*(hstnum(i+kk-1)/10).eq.hstnum(i+kk-1))
     .                       call digit(1,countr+kk,hstnum(i+kk-1),
     .                       amax,array)
                     else
                        array(2,countr+kk) = '-'
                     endif
                     if (j-kk+1.gt.jp) then
                        array(5,countr+kk) = seq(hstnum(j-kk+1))
                        if (10*(hstnum(j-kk+1)/10).eq.hstnum(j-kk+1))
     .                       call digit(6,countr+kk,hstnum(j-kk+1),
     .                       amax,array)
                     else
                        array(5,countr+kk) = '-'
                     endif
                  enddo
                  countr = countr + k
               endif
c      
c              Stacking must occur
c      
               i = ip
               j = jp
               do while (basepr(i).eq.j)
c              Base pair case
                  countr = countr + 1
                  array(1,countr) = ' '
                  array(2,countr) = ' '
                  array(5,countr) = ' '
                  array(6,countr) = ' '
                  array(3,countr) = seq(hstnum(i))
                  if (10*(hstnum(i)/10).eq.hstnum(i))
     .                 call digit(1,countr,hstnum(i),amax,array)
                  array(4,countr) = seq(hstnum(j))
                  if (10*(hstnum(j)/10).eq.hstnum(j))
     .                 call digit(6,countr,hstnum(j),amax,array)
                  i = i + 1
                  j = j - 1
               enddo
               call push(i,j,countr,0)
            endif
         endif
      enddo
c     Normal exit.
      return
      end
 
      subroutine dump_array(amax,array,countr)
      integer amax,countr,k1,k2
      character*1 array(6,amax)
      do k1 = 1,6
         write(6,100) (array(k1,k2),k2=1,countr)
 100     format(5000a1)
      enddo
      write(6,*) ''
      return
      end

c     Puts the number base_num in row 'row' and column 'col' of the array.
c     The least significant digit ends up in column 'col'. If the number
c     is too large to fit, a period is put in column 'col'.
      subroutine digit(row,col,base_num,amax,array)
      implicit integer (a-z)
      integer hstnum(1)
      character*1 array(6,amax)

      position = col
      p = base_num/10
      q = base_num - 10*p
      do while (q.gt.0.or.p.gt.0)
         write(array(row,position),100) q
 100     format(i1)
         if (position.eq.1.and.p.gt.0) then
c     Number overflows to the left. Put a dot in column 'col' and return.
            do k = position,col
               array(row,k) = ' '
            enddo
            array(row,col) = '.'
            return
         endif
         position = position - 1
         r = p
         p = r/10
         q = r - 10*p
      enddo
      return
      end
